import pysam
import collections
import pandas as pd
import os
import pickle
from tagging.get_bb import get_barcode
# https://towardsdatascience.com/pairwise-sequence-alignment-using-biopython-d1a9d0ba861f


def tag_bam(dic, in_consensus, filename_prefix, out_consensus=None, min_repeats=0):
    """
    Function to tag bam files with stats file.
    :param dic: dictoinary containing stats file information from the old conbow pipeline
    :param in_consensus:
    :param filename_prefix:
    :param out_consensus:
    :param min_repeats:
    :return:
    """
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

    #     os.system(f"samtools view -h -F 256 {out_consensus}/{filename_exc}_{min_repeats}.tagged.bam > {out_consensus}/{filename_exc}_{min_repeats}.tagged.bam")
    os.system(
        f"samtools sort -l 9 {out_consensus}/{filename_prefix}_{min_repeats}.tagged.bam -o \
        {out_consensus}/{filename_prefix}_{min_repeats}.tagged.sorted.bam")
    os.system(f"samtools index {out_consensus}/{filename_prefix}_{min_repeats}.tagged.sorted.bam")


def get_stats_tag_bam(in_consensus,
                      filename_prefix,
                      stats_file_path,
                      time_stamp,
                      bb_bam,
                      barcode_pos=[32, 62, 79, 97],
                      barcode_name="BBCR",
                      out_consensus=None,
                      ver="old",
                      write_file=True,
                      min_repeats=0):
    # if output path is not defined, use input path
    if out_consensus == None:
        out_consensus = in_consensus

    # read timestamp file, medaka statistics file
    with open(time_stamp, 'rb') as handle:
        time = pickle.load(handle)
    with open(stats_file_path, 'r') as f:
        count = [line.rstrip() for line in f]

    # get barcode from bb bam file
    barcode_df = get_barcode(bb_bam, barcode_pos=barcode_pos)
    print("time", len(time))
    print("count", len(count))
    print("bar_df", barcode_df.shape)

    count_dic = collections.defaultdict(dict)
    bb_cannot_find = []
    for j, i in enumerate(count):

        readname = i.split(",")[0]
        try:
            count_dic[readname]["bb_all"] = int(i.split(",")[1])
            count_dic[readname]["bb_rev"] = int(i.split(",")[2])
            count_dic[readname]["bb_len_median"] = float(i.split(",")[3])
            count_dic[readname]["ins_all"] = int(float(i.split(",")[4]))
            count_dic[readname]["ins_len_median"] = float(i.split(",")[5])
        except Exception as e:
            print(e, "something is wrong.")
        try:
            count_dic[readname]['timestamp'] = time[readname]
        except Exception as e:
            count_dic[readname]['timestamp'] = "none"
            # print(e, "cannot find time_stamp")
            pass
        try:
            # tag as 'RX' to replace "UMI" in single cell sequencing
            count_dic[readname]['RX'] = barcode_df.loc[readname]['BB']
        except Exception as e:
            #
            count_dic[readname]['RX'] = None
            #            print(e, "cannot find read in barcode_df")
            bb_cannot_find.append(readname)
        try:
            if ver == 'new':
                count_dic[readname]["raw_read_length"] = int(i.split(",")[6])
                count_dic[readname]["RCA_length"] = int(float(i.split(",")[7]))
        except Exception as e:
            print(e, "CSV and Data frame length is not equal. Could be version problem")

    # write sam files
    if write_file is True:
        tag_bam(dic=count_dic, in_consensus=in_consensus, filename_prefix=filename_prefix, out_consensus=out_consensus,
                min_repeats=min_repeats)

    count_df = pd.DataFrame.from_dict(count_dic).T

    return count_df, barcode_df, bb_cannot_find


if __name__ == '__main__':
    in_consensus = '../DER4386'
    out_consensus = "../DER4386/"
    filename_exc = "all-ins"
    stats_file_path = "../DER4386/DER4386-new-stats.csv"
    wg_df = get_stats_tag_bam(in_consensus, filename_exc, stats_file_path, out_consensus=None)

#    DER4622_bc_df_BB25 = get_barcode(bb_bam, barcode_pos=[32, 62, 79, 97], barcode_name="BB25", barcode_fa=None)
