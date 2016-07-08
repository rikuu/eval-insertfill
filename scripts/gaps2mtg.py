import sys

# Gap2Seq defaults
fuz = 10
k = 31

if len(sys.argv) != 2:
    print('Usage:', sys.argv[0], '<gaps.fa>')
    sys.exit(1)

# >chr17:432929-433022
# ttatgcctataatcctagaacttggggaggccaagatgggaNNNNNNNNNNNggtcagaagtttaagatcagcctaagcaacatagcaaaatc

# >bkpt1_chr17_pos_9575_fuzzy_0_HET left_kmer
# CTGAGACGGGCTGCCACCGATTTATTGACGC
# >bkpt1_chr17_pos_9575_fuzzy_0_HET right_kmer
# CGCTGGGACAGTGGCACCGAAACGGAAGGTG

def parse_gap(seq, i):
    split1 = seq[0].split(':')
    contig = split1[0][1:]

    split2 = split1[1].split('-')
    pos = int(split2[0]) + k + fuz

    gap = ''.join([s.rstrip() for s in seq[1:]])
    left = gap[fuz:k+fuz]
    right = gap[-(k+fuz):-fuz]

    assert(len(left) == k)
    assert(len(right) == k)

    # TODO: Figure out fuzzy and HET
    comment = 'bkpt%i_%s_pos_%i_fuzzy_0_HET' % (i, contig, pos)
    print('>%s left_kmer\n%s\n>%s right_kmer\n%s' % (comment, left, comment, right))

with open(sys.argv[1], 'r') as f:
    seq, i = [], 1
    for line in f:
        if line[0] == '>' and len(seq) > 1:
            parse_gap(seq, i)
            seq, i = [], i+1
        seq += [line]
    parse_gap(seq, i)
