import gzip
import pickle
import collections
from dateutil import parser
import os
import argparse
import re
import pysam

def get_timestamp(in_path, out_path):
    time_string = collections.OrderedDict()
    time_stamp = collections.OrderedDict()
    for filename in os.listdir(in_path):
        with pysam.FastxFile(f"{in_path}/{filename}") as fh:
            for entry in fh:
                readname = entry.name
                if not readname:
                    continue
                time_string[readname] = entry.comment.split(" ")[3].split("=")[1][:-1] # tstamp[readname] = time (hh:mm:ss)
                t = parser.isoparse(entry.comment.split(" ")[3].split("=")[1])
                time_stamp[readname] = t.timestamp()
    if out_path is None:
        out_path = f"{in_path}_timestamp.pickle.gz"

    with gzip.open(out_path, 'wb') as handle:
        pickle.dump(time_string, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    pa = argparse.ArgumentParser(description='Get timestamp of nanopore reads from raw fastq files.')
    pa.add_argument('-i', "--in_path", type=str,
                        help='provide exact (not relative) path to fq/fa folder of reads of interests.')
    pa.add_argument('-o', "--out_path", type=str,
                        help='provide path to store timestamp info as pickle.')

    args = pa.parse_args()
    get_timestamp(args.in_path, args.out_path)

    with gzip.open(args.out_path, 'rb') as f:
        file_content = pickle.load(f)
        print("timestamp.pickle.gz is a OrderedDict containing", len(file_content), "readname: timestamp pairs.")