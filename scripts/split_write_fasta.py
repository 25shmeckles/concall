import sys
sys.path.append("/hpc/local/CentOS7/cog/lib/python2.7/site-packages/")
import pysam
import collections
import numpy as np
import pandas as pd
import argparse
import portalocker
from multiprocessing import Pool
# this code is copied from sam.py

parser = argparse.ArgumentParser(description='Split original long read by bowtie mapped backbone (or insert in case of pjet).')
parser.add_argument("--sam", "-s", description="path to samfile containing mapped short reads (backbones) generated by bowtie.")
parser.add_argument("--fasta","-f", description="path to original fasta long read file")
parser.add_argument("--bb", "-b", description="output path to backbone fasta files, containing all subreads of all reads. Format is compatible with medaka smolecule module <read-name>_<sub-read-ID>.")
parser.add_argument("--ins", "-i", description="output path to insert fasta files, containing all subreads of all reads. Format is compatible with medaka smolecule module <read-name>_<sub-read-ID>.")
parser.add_argument("--stats","-st", description="path to store overall stats for long reads.")
parser.add_argument("--st_sub","ss", description="path to stats per subreads.")
parser.add_argument("--length","-l", description="extract only the reads with certain length between the matches. Given as a tuple (mix_gap, max_gap). This option is deprecated since it is covered by medaka pipeline downstream.")
args = parser.parse_args()

samFile = pysam.AlignmentFile(args.sam, "rb")
inFile = args.fasta
bb_outFile = args.bb
ins_outFile = args.ins
statFile = args.stats
sub_statFile = args.st_sub



# Create a dict with lists of SimpleRead using referenceName as key

myDict = collections.defaultdict(list)
for read in samFile:
    # Remove matches that suck
    alnErr = sum(x[1] for x in read.cigartuples if x[0] in [1,2])
    if alnErr/float(read.reference_end) < 0.2:
        #read.reference_start + read.infer_query_length(always=True)
        myDict[read.reference_name].append((read.reference_start,read.reference_end,read.flag&(1<<4)>0,read.query_name))

# Combine matches per read
sortedKeys =list( myDict.keys())
sortedKeys.sort()
for key in sortedKeys:
    myDict[key].sort()
    # print(key,myDict[key])

# exit()
min_gap = 50
max_gap = 500
# Work through the fastq file

read_stats = collections.namedtuple("read_name","read_seq","read_plus","read_phred","read_id")
ins_subread_stats = collections.namedtuple("read_name","subread_repeat","direction","ref_name","start_position","end_position") 
bb_subread_stats = collections.namedtuple("read_name","subread_repeat","direction","ref_name","start_position","end_position")
bb_fastafile = ""
ins_fastafile = ""
with portalocker.Lock("", timeout=1) as statFile:

#while True:
#    readName=fqFile.next().rstrip()
#    readSeq=fqFile.next().rstrip()
#    readPlus=fqFile.next().rstrip()
#    readPhred=fqFile.next().rstrip()
#    readId=readName[1:].split()[0]

with open(inFile,'r') as fqFile, open(ins_outFile,'w') as dumpFile, \
        open(bb_outFile,'w') as dump_bb_File, open(bb_outFile+'_stat.csv','w') as statFile:
    try:
        while True:
            readName=fqFile.next().rstrip()
            readSeq=fqFile.next().rstrip()
            readPlus=fqFile.next().rstrip()
            readPhred=fqFile.next().rstrip()
            readId=readName[1:].split()[0]


            cutId = 1
            if readId not in myDict:
                statFile.write(', '.join(
                    [str(y) for y in
                        [readId,
                        0,
                        0,
                        0,
                        0,
                        0]])+'\n')
                continue

            cur_bbs = myDict[readId]
            # Split sequence, dump backbone
            for i,x in enumerate(cur_bbs):
                # print(x[1]-x[0])
                # Extend the identifier
                dump_bb_File.write(readName + ' c_start=' + str(x[0]) + ' c_ref=' + x[3] + ' c_str=' + ('-' if x[2] else '+') + ' c_idx=' + str(i) + '\n')#+':'+str(x[0])+'-'+str(x[1])

                # Dump the actual sub sequence with primers
                dump_bb_File.write(readSeq[x[0]:x[1]] + '\n')
                dump_bb_File.write(readPlus+'\n')

                # Add perfect phred scores for forced primers
                dump_bb_File.write(readPhred[x[0]:x[1]] + '\n')

            # Mind the gaps, should be inserts
            gaps = []
            if min_gap < cur_bbs[0][0] < max_gap:
                gaps.append((0,cur_bbs[0][0]))

            for i,x in enumerate(cur_bbs[1:]):
                if min_gap < x[0] - cur_bbs[i][1] < max_gap:
                    gaps.append((cur_bbs[i][1],x[0]))

            if min_gap < len(readSeq) - cur_bbs[-1][1] < max_gap:
                gaps.append((cur_bbs[-1][1],len(readSeq)))

            for i,x in enumerate(gaps):
                # Extend the identifier
                dumpFile.write(readName + ' c_start=' + str(x[0]) + ' c_ref=gap' + ' c_idx=' + str(i) + '\n')#+':'+str(x[0])+'-'+str(x[1])

                # Dump the actual sub sequence with primers
                dumpFile.write(readSeq[x[0]:x[1]] + '\n')
                dumpFile.write(readPlus+'\n')

                # Add perfect phred scores for forced primers
                dumpFile.write(readPhred[x[0]:x[1]] + '\n')

            statFile.write(', '.join(
                [str(y) for y in
                    [readId,
                    len(cur_bbs),
                    sum([x[2] for x in cur_bbs]),
                    np.median([x[1]-x[0] for x in cur_bbs]),
                    len(gaps),
                    (np.median([x[1]-x[0] for x in gaps]) if len(gaps)>0 else 0)]])+'\n')

    except StopIteration:
        pass

print('Splitting data finished.')
