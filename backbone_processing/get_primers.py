import sys
import argparse
complement = str.maketrans('ATCGN', 'TAGCN')

def check_N(seq):
    return seq.count("N")

def reverse_complement(seq):
    k = list(seq.translate(complement))
    k.reverse()
    rev_comp = "".join(k)
    return rev_comp


def prepare_primer_sequences(bb_seq:str,
                             base_count_3prime:int=20,
                             base_count_5prime:int=20,
                             skip_count_3prime:int=0,
                             skip_count_5prime:int=0,
                             write_name=False):
    """

    :param bb_seq: backbone sequence
    :param base_count_3prime: how many base pairs do you want to match on the 3 prime end.
    :param base_count_5prime: how many base pairs do you want to match on the 5 prime end.
    :param skip_count_3prime: how many base pairs do you want to skip on 3 prime end (useful when the 3' and 5' end has
    complement sequences.
    :param skip_count_5prime: how many base pairs do you want to skip on 5 prime end (useful when the 3' and 5' end has
    complement sequences.
    :param write: write 5' and 3' primer to ../data/seg/{write_name}.fa if name is provided; if name is not provided,
    the file will not be written.
    :return: (3 prime seq, 5 prime seq) as a tuple.

    >>> prepare_primer_sequences("GGGCATGCACAGATGTACACGTGACGCAACGANTGATGTTAGCTATTTGTTCAATGACATATNCTGGTATGATCAATACN\
    AGATCTGATATTGATATNCTGATACTCATATATGTAGAATATCACATTATATTTATTATAATACATCGTCGAACATATACACAATGCATCTTATCTATACGTATCGGGATA\
    GCGTTGGCATAGCACTGGATGGCATGACCCTCATTAGATGCTGCATGACATAGCCC", 20, 20, 0 ,0, "bb22")
    ('GATGCTGCATGACATAGCCC', 'GTGTACATCTGTGCATGCCC')
    """

    total_len = len(bb_seq)
    seq_3prime = bb_seq[total_len-(base_count_3prime+skip_count_3prime): total_len-skip_count_5prime]
    seq_5prime = reverse_complement(bb_seq[skip_count_5prime:(skip_count_5prime+ base_count_5prime)])
    fail = False
    if check_N(seq_3prime) != 0:
        print("Warning: 3' primer sequence contains N. May cause problem in downstream analysis.")
        fail = True
    if check_N(seq_5prime) != 0:
        print("Warning: 5' primer sequence contains N. May cause problem in downstream analysis."
              " Try a shorter sequence.")
        fail = True
    if fail:
        sys.exit(1)

    if write_name:
        with open(f"./data/seg/3_prime_{write_name}.fa", "w") as f:
            f.write("\n".join([f">3_prime_{write_name}",seq_3prime]))
        with open(f"./data/seg/5_prime_{write_name}.fa", "w") as f:
            f.write("\n".join([f">5_prime_{write_name}",seq_5prime]))
    return (seq_3prime, seq_5prime)


def get_default_seq(bb_name):
    '''

    :param bb_name:
    :return:
    >>> get_default_seq("BBCR")[:34]
    'GGGCGGTATGTCATGCACACGAATCCCGAAGANT'
    '''
    with open(f"./data/backbones/backbones.fa", "r") as f:
        flist = [x.strip() for x in f.readlines()]
        i = flist.index(f">{bb_name}")
        return flist[i+1]




if __name__ == "__main__":
    # import doctest
    # doctest.testmod()
    # Example usage: python get_primers.py --bb_type BB36 --prime3 25 --prime5 20 --prime5_skip 4 --prime3_skip 4
    parser = argparse.ArgumentParser(description='Get primer sequences.')
    parser.add_argument("--bb_type", help="backbone name (according to the backbones.fa file.", default=None)
    parser.add_argument("--bb_seq", help="backbone sequence. Only provide this when it is not listed in backbones.fa file", default=None)
    parser.add_argument("--bb_name", help="provide bb_name when bb_seq is provided. Will be used to write new primer files.", default=None)
    parser.add_argument("--prime3", help="length of primer to produce", type=int, default=20)
    parser.add_argument("--prime3_skip", help="how many base pairs do you want to skip on 5 prime end (useful "
                                             "when the 3' and 5' end has complement sequences", type=int, default=0)
    parser.add_argument("--prime5", help="length of primer to produce,", type=int, default=20)
    parser.add_argument("--prime5_skip", help="how many base pairs do you want to skip on 5 prime end (useful "
                                             "when the 3' and 5' end has complement sequences", type=int, default=0)
    args = parser.parse_args()
    if args.bb_type:
        bb_seq = get_default_seq(args.bb_type)
        bb_name = args.bb_type
    else:
        bb_seq = args.bb_seq
        bb_name = args.bb_name
    prepare_primer_sequences(bb_seq=bb_seq,
                         base_count_3prime=args.prime3,
                         base_count_5prime=args.prime5,
                         skip_count_3prime=args.prime3_skip,
                         skip_count_5prime=args.prime5_skip,
                         write_name=bb_name)


