import re
import collections
import pysam
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns
import time
import progressbar
import subprocess
import numpy as np


def convert_base_to_num(value):
    if value == "T":
        out = 4
    elif value == "C":
        out = 1
    elif value == "A":
        out = -1
    elif value == "G":
        out = -4
    else:
        out = -100
    return out


def convert_num_to_base(value):
    if value == 4:
        out = "T"
    elif value == 1:
        out = "C"
    elif value == -1:
        out = "A"
    elif value == -4:
        out = "G"
    elif value == -100:
        out = "X"
    else:
        out = value
    return out


class Barcode:
    def __init__(self,
                 bb_name=None,
                 bb_only_bam=None,
                 out_path=None,
                 bb_fa_path="../../concall/data/ref/BBCR_ref.fa",
                 pos=None,
                 max_depth=100000000):
        self.bb_name = bb_name
        self.sample_bam = bb_only_bam
        assert bb_only_bam is not None
        assert os.path.exists(self.sample_bam), f"{self.sample_bam} does not exist."
        # TODO: if not sorted, sort and index
        self.bb_fa_path = bb_fa_path
        assert os.path.exists(self.bb_fa_path), f"{self.bb_fa_path} does not exist."
        self.pos = pos
        self.out_bb_pickle = None
        if self.pos is None:
            self.get_ref_barcode_positions()
            print(self.pos)
        try:
            self.mpileup_range = [self.pos[0], self.pos[-1]]
        except Exception as e:
            print(e)
        # out path -- a directory to hold 3 files
        self.out_path = out_path
        # pickle file
        self.out_bb_pickle = None
        self.bb_single = None
        self.read_count = None
        self.real_read_count = None
        self.max_depth = max_depth
        self.get_read_count()
        self.read_type = collections.Counter()

    def get_read_count(self):
        p = subprocess.check_output(["samtools", "view", "-c", "-F", "260", self.sample_bam])
        p = p.decode("utf-8").strip()
        self.real_read_count = int(p)

    def get_ref_barcode_positions(self):
        with open(self.bb_fa_path, 'r') as f:
            bb = f.readlines()
        bb = bb[1].strip()
        bb = bb.upper()
        bb_pos = [match.start() for match in re.finditer(re.escape("N"), bb)]
        self.pos = bb_pos

    def get_barcode_from_consensus(self):
        barcode = collections.defaultdict(lambda: collections.defaultdict(list))
        if self.max_depth > self.real_read_count:
            max_value = self.real_read_count
        else:
            max_value = self.max_depth
        bar = progressbar.ProgressBar(max_value=max_value)
        i = 0
        with pysam.AlignmentFile(self.sample_bam, "rb") as samfile:
            for pileupcolumn in samfile.pileup(self.bb_name,
                                               self.mpileup_range[0],
                                               self.mpileup_range[1],
                                               max_depth=self.max_depth):
                if (pileupcolumn.pos == self.pos[0]) or \
                        (pileupcolumn.pos == self.pos[1]) or \
                        (pileupcolumn.pos == self.pos[2]) or \
                        (pileupcolumn.pos == self.pos[3]):
                    for pileupread in pileupcolumn.pileups:

                        # This is an estimate progress bar. It's not precise!
                        i += 0.25
                        if i > self.max_depth:
                            i = self.max_depth
                        bar.update(i)
                        if not pileupread.is_del and not pileupread.is_refskip:
                            self.read_type['base'] += 1
                            read_name_per_rep = pileupread.alignment.query_name
                            consensus_name = pileupread.alignment.query_name.split("_")[0]
                            base_rep = convert_base_to_num(
                                pileupread.alignment.query_sequence[pileupread.query_position])
                            barcode[read_name_per_rep][f"b{pileupcolumn.pos}"].append(base_rep)
                            barcode[read_name_per_rep]["read_name"] = consensus_name
                            # query position is None if is_del or is_refskip is set.
                        elif pileupread.is_del or pileupread.is_refskip:
                            self.read_type['indel'] += 1
                            barcode[read_name_per_rep][f"b{pileupcolumn.pos}"].append(0)
                            barcode[read_name_per_rep]["read_name"] = consensus_name

                        # if it doesn't map to this base
                        else:
                            self.read_type['nan'] += 1
                            barcode[read_name_per_rep][f"b{pileupcolumn.pos}"].append(0)
                            barcode[read_name_per_rep]["read_name"] = consensus_name
        single_barcode = collections.defaultdict(lambda: collections.defaultdict(int))
        for k, v in barcode.items():
            for k2, v2 in v.items():
                if type(v2) == list:
                    if len(v2) > 1:
                        v3 = [x for x in v2 if x != 0]
                        if len(set(v3)) > 1:
                            # chimeric
                            self.read_type['chimeric'] += 1
                            single_barcode[k][k2] = -15
                        else:
                            # list exclude 0 has only len == 1
                            self.read_type['multiple-alignment'] += 1
                            single_barcode[k][k2] = v3[0]
                    else:
                        self.read_type['1-mapped-base'] += 1
                        single_barcode[k][k2] = v2[0]
        self.bb_single = pd.DataFrame.from_dict(single_barcode).T
        self.bb_single.fillna(-20, inplace=True)

        assert self.bb_single.isna().value_counts().shape == (1,)
        self.out_bb_pickle = pd.DataFrame.from_dict(barcode).T
        self.out_bb_pickle = self.out_bb_pickle[[f"b{self.pos[0]}",
                                                 f"b{self.pos[1]}",
                                                 f"b{self.pos[2]}",
                                                 f"b{self.pos[3]}",
                                                 "read_name"]]
        self.out_bb_pickle.to_pickle(f"{self.out_path}.pickle.gz")
        pd.to_pickle(self.bb_single.value_counts(), f"{self.out_path}_value_counts.pickle.gz")

        self.read_count = self.bb_single.shape[0]
        print(f"\nSamtools pileup: {self.real_read_count} reads. Pysam end results: {self.read_count} reads.")
        print(self.read_type)
        self.bb_single.to_pickle(f"{self.out_path}_bb_single.pickle.gz")

    def view_barcode(self):
        if self.bb_single is None:
            if os.path.exists(f"{self.out_path}_bb_single.pickle.gz"):
                self.bb_single = pd.read_pickle(f"{self.out_path}_bb_single.pickle.gz")
                assert self.bb_single.isna().value_counts().shape == (1,)
            else:
                print("The barcode is not extracted yet. Please run get_barcode_from_consensus()")
                exit()
        # view barcode distribution
        # nan: -20, chimera: -15, bases: -4:G, -1:A, 1:C, 4:T
        # mask < -20 ==> not masked
        cg = sns.heatmap(self.bb_single.mask(self.bb_single < -20),
                         cmap="RdYlGn",
                         cbar_kws={"ticks": [-20, -10, -4, -1, 1, 4]},
                         yticklabels=False
                         )
        plt.savefig(f"{self.out_path}.png", dpi=150)


def get_barcode(bb_bam, barcode_pos, barcode_name, barcode_fa=None):
    barcode = collections.defaultdict(lambda: collections.defaultdict(list))
    # random positions, 0-based: [32, 62, 79, 97]
    start = barcode_pos[0]
    end = barcode_pos[-1]
    if os.path.isdir(bb_bam):
        for file in os.listdir(bb_bam):
            if file.endswith("sorted.bam"):
                # file path + file name
                samfile = pysam.AlignmentFile( bb_bam + file, "rb")
                for pileupcolumn in samfile.pileup(barcode_name, start, end):
                    if (pileupcolumn.pos == 32) or (pileupcolumn.pos == 62) or (pileupcolumn.pos == 79) or \
                    (pileupcolumn.pos == 97):
                        for pileupread in pileupcolumn.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                read_name = (pileupread.alignment.query_name).split("_")[0]
                                base_rep = convert_base_to_num(pileupread.alignment.query_sequence
                                [pileupread.query_position])
                                barcode[read_name][f"b{pileupcolumn.pos}"].append(base_rep)
                                # query position is None if is_del or is_refskip is set.
                            elif pileupread.is_del or pileupread.is_refskip:
                                barcode[read_name][f"b{pileupcolumn.pos}"].append(0)
                            else:
                                barcode[read_name][f"b{pileupcolumn.pos}"].append(0)
                                # barcode[read_name][f"b{pileupcolumn.pos}_range7"] =
                                # pileupread.alignment.query_sequence[pileupread.query_position -
                                # 3:pileupread.query_position + 3]

                samfile.close()
    bb_combi = pd.DataFrame.from_dict(barcode).T
    bb_combi.fillna(-100, inplace=True)
    bb_info = collections.defaultdict(dict)

    for i, row in bb_combi.iterrows():
        a = []
        most_freq = ""
        for x in [row['b32'], row['b62'], row['b79'], row['b97']]:
            if x == -100:
                x = [0, 0]
            else:
                pass
            try:
                most_common_2 = collections.Counter(x).most_common(2)
                most_common = most_common_2[0][0]
                if most_common == 0:
                    if len(most_common_2) == 1:
                        most_freq += "X"
                    else:
                        # pick the second common base rather than X
                        most_freq += convert_num_to_base(most_common_2[1][0])
                else:
                    most_freq += convert_num_to_base(most_common)
            except Exception as e:
                print(e)
                print(x)
            a.append(np.asarray(x))
        bb_info[i]['bc'] = most_freq
        max_len = max([len(x) for x in a])
        signal = [np.pad(x, (0, max_len - len(x)), 'constant') for x in a]
        new_signal = np.vstack(signal)
        break_points = get_chimeric_reads(i, new_signal.T, plot=False)
        if len(break_points) > 0:
            bb_info[i]['break_points'] = break_points
        else:
            bb_info[i]['break_points'] = False
    bb_info_df = pd.DataFrame.from_dict(bb_info).T
    return bb_info_df


def get_chimeric_reads(read_name, signal, fig_out="../data/bc_breakpoint", plot = False):
    algo = rpt.Pelt(model="rbf").fit(signal)
    result = algo.predict(pen=10)
    rpt.display(signal, [0,1], result)
    if plot:
        if not os.path.exists(fig_out):
            os.makedirs(fig_out)
        plt.savefig(f"{fig_out}/{read_name[:6]}.png")
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--bb_name", help="Backbone type. This determines where the barcode is located",
                        type=str, default="BBCR")
    parser.add_argument("--bb_only_bam", help="bam file with the backbone sequences consensus.", type=str)
    parser.add_argument("--bb_bam", help="bam file with the backbone sequences mapped with bowtie. directory 01_bowtie",
                        type=str)
    parser.add_argument("--bb_fa_path", help="directory where bb_ref.fa file located")
    parser.add_argument("--out_path", help="out readname to barcode pickle file, has to end with .pickle.gz")
    parser.add_argument("--testing", action='store_true')
    parser.add_argument("--view", action='store_true')
    parser.add_argument("--dont_get_bb", action='store_true')
    args = parser.parse_args()
    # read name is the contig name. Use fasta to map to barcode again. Construct reference with only barcode
    # if testing: subset reads via pysam.
    if args.testing is True:
        set_max_depth = 1000
    else:
        set_max_depth = 100000000
    # current version: get bb from consensus
    p2 = Barcode(args.bb_name,
                 args.bb_only_bam,
                 args.out_path,
                 args.bb_fa_path,
                 max_depth=set_max_depth)
    a = time.time()
    if not args.dont_get_bb:
        p2.get_barcode_from_consensus()
    if args.view:
        p2.view_barcode()
    b = time.time()
    print(f"\nProcess finished without error. Total: {p2.read_count} reads. Took {round((b-a)/60, 2)} mins.")
    # previous version: get bb from individual maps: get_barcode
    # bb_info_df = get_barcode(args.bb_bam, barcode_pos=[32, 62, 79, 97], barcode_name="BB24", barcode_fa=None)
    # bb_info_df.to_pickle("../data/bb_info.pickle")
