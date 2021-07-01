import os
import pysam
import uuid
import contextlib
from shutil import which


@contextlib.contextmanager
def sorted_bam_file(
        write_path,
        origin_bam=None,
        header=None,
        read_groups=None,
        local_temp_sort=True,
        input_is_sorted=False,
        mode='wb',
        fast_compression=False,  # Use fast compression for merge (-1 flag)
        temp_prefix = 'concall-TMP',
        **kwargs
        ):
    """ Get writing handle to a sorted bam file

    Args:
        write_path (str) : write to a  bam file at this path

        origin_bam (pysam.AlignmentFile ) : bam file to copy
         header for or
        header (dict) : header for the bam file to write

        read_groups(set/dict) : set or dictionary which contains read groups. The dictionary should have the format { read_group_id (str)
                { 'ID': ID, 'LB':library,
                'PL':platform,
                'SM':sampleLib,
                'PU':readGroup }

        local_temp_sort(bool) : create temporary files in current directory

        input_is_sorted(bool) : Assume the input is sorted, no sorting will be applied

        mode (str) : Output mode, use wbu for uncompressed writing.

        **kwargs : arguments to pass to the new pysam.AlignmentFile output handle


    Write some pysam reads to a sorted bam file
    Example:
        >>> import pysam
        >>> from tagging.bamfunctions import sorted_bam_file
        >>> test_sam = pysam.AlignmentFile('test.sam','w',reference_names=['chr1','chr2'],reference_lengths=[1000,1000])
        >>> read_A = pysam.AlignedSegment(test_sam.header)
        >>> read_A.reference_name = 'chr2'
        >>> read_A.reference_start = 100
        >>> read_A.query_sequence = 'TTGCA'
        >>> read_A.query_name= 'READ_A'
        >>> read_A.cigarstring = '5M'
        >>> read_A.query_qualities = [30] * len(read_A.query_sequence)
        >>> read_A.set_tag('RG','HVKCCBGXB.4.MYLIBRARY_1')
        >>> read_B = pysam.AlignedSegment(test_sam.header)
        >>> read_B.reference_name = 'chr1'
        >>> read_B.reference_start = 100
        >>> read_B.query_sequence = 'ATCGGG'
        >>> read_B.cigarstring = '6M'
        >>> read_B.query_name= 'READ_B'
        >>> read_B.query_qualities = [30] * len(read_B.query_sequence)
        >>> read_B.set_tag('RG','HVKCCBGXB.4.MYLIBRARY_2')
        >>> read_groups = set(( 'HVKCCBGXB.4.MYLIBRARY_2','HVKCCBGXB.4.MYLIBRARY_1'))
        >>> with sorted_bam_file('out.bam', header=test_sam.header,read_groups=read_groups) as out:
        >>>     out.write(read_A)
        >>>     out.write(read_B)
        Results in the bam file:
        @HD	VN:1.6	SO:coordinate
        @SQ	SN:chr1	LN:1000
        @SQ	SN:chr2	LN:1000
        @RG	ID:HVKCCBGXB.4.MYLIBRARY_2	SM:MYLIBRARY_2	LB:MYLIBRARY	PU:HVKCCBGXB.4.MYLIBRARY_2	PL:ILLUMINA
        @RG	ID:HVKCCBGXB.4.MYLIBRARY_1	SM:MYLIBRARY_1	LB:MYLIBRARY	PU:HVKCCBGXB.4.MYLIBRARY_1	PL:ILLUMINA
        READ_B	0	chr1	101	0	6M	*	0	0	ATCGGG	??????	RG:Z:HVKCCBGXB.4.MYLIBRARY_2
        READ_A	0	chr2	101	0	5M	*	0	0	TTGCA	?????	RG:Z:HVKCCBGXB.4.MYLIBRARY_1


    """
    unsorted_path = None
    unsorted_alignments = None
    target_dir = None

    unsorted_path = f'{write_path}.unsorted'
    # Create output folder if it does not exists
    target_dir = os.path.dirname(unsorted_path)
    if not os.path.exists(target_dir) and len(target_dir)>0 and target_dir!='.':
        try:
            os.makedirs(target_dir, exist_ok=True)
        except Exception as e:
            pass

    if header is not None:
        pass
    elif type(origin_bam) is str:
        with pysam.AlignmentFile(origin_bam) as ain:
            header = ain.header.copy()
    elif origin_bam is not None:
        header = origin_bam.header.copy()
    else:
        raise ValueError("Supply a header or origin_bam object")


    try: # Try to open the temp unsorted file:
        unsorted_alignments = pysam.AlignmentFile(
            unsorted_path, mode, header=header, **kwargs)
    except Exception:
        # Raise when this fails
        raise

    # Yield a handle to the alignments,
    # this handle will be released when the handle runs out of scope
    yield unsorted_alignments
    unsorted_alignments.close()

    if read_groups is not None:
        add_readgroups_to_header(unsorted_path, read_groups)

    # Write, sort and index
    if input_is_sorted is False:
        sort_and_index(
            unsorted_path,
            write_path,
            remove_unsorted=True,
            local_temp_sort=local_temp_sort,
            fast_compression=fast_compression,
            prefix=temp_prefix
            )
    else:
        os.rename(unsorted_path, write_path)


def sort_and_index(
        unsorted_path,
        sorted_path,
        remove_unsorted=False,
        local_temp_sort=True,
        fast_compression=False,
        prefix='TMP'
        ):
    """ Sort and index a bam file
    Args:
        unsorted_path (str) : path to unsorted bam file

        sorted_path (str) : write sorted file here

        remove_unsorted (bool) : remove the unsorted file

        local_temp_sort(bool): create temporary files in target directory

        fast_compression(bool): fast compression

        prefix(str): prefix for temp file. default: 'TMP'
    Raises:
        SamtoolsError when sorting or indexing fails
    """

    if prefix is None:
        raise Exception('prefix cannot be None')

    if local_temp_sort:
        base_directory = os.path.abspath(os.path.dirname(sorted_path))
        prefix = f'{prefix}.{uuid.uuid4()}'
        temp_path_first = f'{base_directory}/{prefix}'
        if temp_path_first.startswith('/TMP'):
            # Perform sort in current directory
            temp_path_first = f'./{prefix}'

        # Try to sort at multiple locations, if sorting fails try the next until all locations have been tried
        temp_paths = [temp_path_first, f'/tmp/{prefix}', f'./{prefix}']
        for i, temp_path in enumerate(temp_paths):
            failed = False
            try:
                pysam.sort(
                    '-o',
                    sorted_path,
                    '-T',
                    f'{temp_path}',  # Current directory with a random prefix
                    unsorted_path, ('-l 1' if fast_compression else '-l 3')
                )
            except Exception as e:
                print(f'Sorting failed at temp: {temp_path}. {e}')
                failed = True
                if i == len(temp_paths)-1:
                    raise

            if not failed:
                break
    else:
        pysam.sort("-o", sorted_path, unsorted_path, ('-l 1' if fast_compression else '-l 3'))
    pysam.index(sorted_path)
    if remove_unsorted:
        os.remove(unsorted_path)


def add_readgroups_to_header(
        origin_bam_path,
        readgroups_in,
        target_bam_path=None,
        header_write_mode='auto'):
    """ Add the readgroups in the set readgroups to the header of origin_bam_path.

    This function first loads the header of the origin to memory.
    The supplied readgroups are added to this header.
    The new header is then exported to a SAM file. The SAM file is then
    concatenated to the original bam file.

    Args:
        origin_bam_path(str) : path to bam file to which to add readgroups to header

        readgroups_in(set/dict) : set or dictionary which contains read groups. The dictionary should have the format { read_group_id (str)
                { 'ID': ID, 'LB':library,
                'PL':platform,
                'SM':sampleLib,
                'PU':readGroup }

        target_bam_path(str) : path to write bam file including the readgrouped header to. When not supplied the output is written to the input bam file

    """

    # Create a read group dictionary
    if isinstance(readgroups_in, set):
        readGroupsDict = {}
        for readGroup in readgroups_in:
            flowCell, lane, sampleLib = readGroup.split('.')
            try:
                library, _ = sampleLib.rsplit('_', 1)
            except Exception as e:
                # the library is not part of the sample name:
                library = 'undefinedLibrary'
            readGroupsDict[readGroup] = {
                'ID': readGroup,
                'LB': library,
                'PL': 'ILLUMINA',
                'SM': sampleLib,
                'PU': readGroup}
    elif isinstance(readgroups_in, dict):
        readGroupsDict = readgroups_in
    else:
        raise ValueError("supply a set or dict for readgroups_in")

    with pysam.AlignmentFile(origin_bam_path, "rb") as origin:
        header = origin.header.copy()
        hCopy = header.to_dict()
        hCopy['RG'] = list(readGroupsDict.values())
        replace_bam_header(origin_bam_path, hCopy,
                            target_bam_path=target_bam_path,
                            header_write_mode=header_write_mode)



def replace_bam_header(origin_bam_path, header, target_bam_path=None, header_write_mode='auto'):

    if target_bam_path is None:
        target_bam_path = origin_bam_path

    # Write the re-headered bam file to this path
    complete_temp_path = origin_bam_path.replace('.bam', '') + '.rehead.bam'

    # When header_write_mode is auto, when samtools is available, samtools
    # will be used, otherwise pysam
    if header_write_mode == 'auto':
        if which('samtools') is None:
            header_write_mode = 'pysam'
        else:
            header_write_mode = 'samtools'

    if header_write_mode == 'pysam':

        with pysam.AlignmentFile(complete_temp_path, "wb", header=header) as out, \
             pysam.AlignmentFile(origin_bam_path) as origin:
            for read in origin:
                out.write(read)

        os.rename(complete_temp_path, target_bam_path)

    elif header_write_mode == 'samtools':

        # Write the new header to this sam file:
        headerSamFilePath = origin_bam_path.replace(
            '.bam', '') + '.header.sam'

        # Write the sam file with the complete header:
        with pysam.AlignmentFile(headerSamFilePath, 'w', header=header):
            pass

        # Concatenate and remove origin
        rehead_cmd = f"""{{ cat {headerSamFilePath}; samtools view {origin_bam_path}; }} | samtools view -b > {complete_temp_path} ;
                mv {complete_temp_path } {target_bam_path};rm {headerSamFilePath};
                """
        os.system(rehead_cmd)
    else:
        raise ValueError(
            'header_write_mode should be either, auto, pysam or samtools')

