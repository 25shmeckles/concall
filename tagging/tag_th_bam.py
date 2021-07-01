import pysam
import pandas as pd
import os
import argparse
from time import time
import progressbar
from bamfunctions import sorted_bam_file


def convert_num_to_base(value):
    if value == 4:
        out = "T"
    elif value == 1:
        out = "C"
    elif value == -1:
        out = "A"
    elif value == -4:
        out = "G"
    else:
        out = "X"
    return out


class Bamprocessing:
    def __init__(self,
                 mapped_bam_path,
                 bb_single_df,
                 timestamp_df,
                 keep_original_name=False,
                 export_df=True,
                 ):
        self.mapped_bam_path = mapped_bam_path
        self.bb_single_df = pd.read_pickle(bb_single_df).applymap(convert_num_to_base) # Pandas DataFrame
        self.timestamp_df = pd.read_pickle(timestamp_df) # OrderedDict
        self.keep_original_name = keep_original_name
        self.export_df = export_df
        self.prefix = "".join(self.mapped_bam_path.split(".")[:-2])
        # TODO: check if indexed
        pysam.index(self.mapped_bam_path)
        rc = int(pysam.view("-c", self.mapped_bam_path).strip())
        self.bar = progressbar.ProgressBar(max_value=rc)

    def tag_th_bam(self):
        with pysam.AlignmentFile(self.mapped_bam_path) as inputfile:
            with sorted_bam_file(f"{self.prefix}.tagged.sorted.bam", header=inputfile.header) as out:
                for i, read in enumerate(inputfile):
                    self.bar.update(i)
                    params = read.query_name.rsplit(
                        '_', 10)
                    read.set_tag('RN', params[0])  # read name
                    read.set_tag('UC', params[1])  # unique consensus id "rep0"
                    read.set_tag('RC', params[6])  # repeat count
                    read.set_tag('FL', params[2])  # full read length
                    read.set_tag('MS', params[7])  # matched score
                    # Time stamp
                    try:
                        read.set_tag('TI', self.timestamp_df[read.get_tag('RN')])  # timestamp
                    except KeyError:
                        # logging
                        read.set_tag('TI', "NotFound")
                    # Barcode
                    try:
                        read.set_tag('BX', "".join(self.bb_single_df.loc[read.query_name].values))
                    except KeyError as e:
                        # logging
                         read.set_tag('BX',"NotFound")
                    if self.keep_original_name is True:
                        pass
                    else:
                        # how to get read name?
                        if read.is_unmapped:
                            read.query_name = f"unmapped_{i}"
                        else:
                            # TODO: Test set has different query name format to DER4498_P260
                            read.query_name = f"{read.reference_name.split('-')[0]}:{read.pos}_{i}"
                    out.write(read)



def tag_bam(dic, in_consensus, filename_prefix, out_consensus=None, min_repeats=0):

    if out_consensus == None:
        out_consensus = in_consensus

    g = pysam.AlignmentFile(f"{in_consensus}/{filename_prefix}.sorted.bam")
    with pysam.AlignmentFile(f"{out_consensus}/{filename_prefix}_{min_repeats}.tagged.bam", 'wb',
                             header=g.header) as sub_file:
        for read in g:
            read_name = read.query_name.split("_")[0]
            try:
                read.set_tag('BC', dic[read_name]['bb_all'])  # bb repeat count
                read.set_tag('IC', dic[read_name]['ins_all'])  # ins repeat count
                read.set_tag('BM', dic[read_name]['bb_len_median'])  # length
                read.set_tag('IM', dic[read_name]['ins_len_median'])  # length
                read.set_tag('RV', dic[read_name]['bb_rev'])  # rev
                read.set_tag("FL", dic[read_name]['raw_read_length'])  # full read length
                read.set_tag('TI', dic[read_name]['timestamp'])  # timestamp
                read.set_tag('RX', dic[read_name]['RX'])  # Backbone 4 bp sequence

            except Exception as e:
                read.set_tag('BB', "None")  # bb repeat count
            #                 print(e, "cannot find read in count dict.")
            try:
                dic[read_name]['chr'] = read.reference_name
                dic[read_name]['positions'] = (read.get_reference_positions()[0], read.get_reference_positions()[-1])
                seq = read.get_reference_sequence()
                dic[read_name]['sequence'] = seq
                dic[read_name]['mapped_read_length'] = len(seq)
            except Exception as e:
                # print(e) # list index out of range
                # upmapped reads
                pass
            if dic[read_name]['bb_all'] >= min_repeats:
                sub_file.write(read)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapped_bam_path", help="bam file", type=str)
    parser.add_argument("--bb_single", help="bb_single dataframe", type=str)
    parser.add_argument("--timestamp_df", help="timestamp for tagging",
                        type=str, default="bb_single dataframe for tagging")
    args = parser.parse_args()
    a = time()
    p = Bamprocessing(args.mapped_bam_path, args.bb_single, args.timestamp_df)
    p.tag_th_bam()
    b = time()
    print(f"\nTagging bam file finished. Took {(b-a)/60} minutes.")