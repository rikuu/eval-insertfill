import sys
import numpy as np
from math import log

if len(sys.argv) != 8:
    print('Usage:', sys.argv[0], '<genome.fa>',
        '<num_gaps> <min_length> <max_length>',
        '<inserts.bed> <gaps.bed> <breakpoints.bed>')
    sys.exit(1)

genome_file = sys.argv[1]
num_gaps = int(sys.argv[2])
length_start = float(sys.argv[3])
length_end = float(sys.argv[4])
inserts_file = sys.argv[5]
gaps_file = sys.argv[6]
breakpoints_file = sys.argv[7]

# Default k+fuz from Gap2Seq
flank_length = 41

genome_length = 0
contig = ''

# Parse genome length and contig identifier
# NOTE: Assumes genome consists of a single contig
with open(genome_file, 'r') as f:
    for line in f:
        if line[0] == '>':
            contig = line[1:-1]
        else:
            genome_length += len(line) - 1

# Generate gap lengths in logarithmic space
lengths = np.logspace(log(length_start, 10), log(length_end, 10), num_gaps)

# Generate random starting positions for gaps
def gen(length):
    start = np.random.randint(flank_length, genome_length - flank_length)
    end = start + int(length)
    return (start, end)

# Check for overlapping gaps
def overlap(gaps, gap):
    s1, e1 = gap[0], gap[1]
    for g in gaps:
        s2, e2 = g[0], g[1]
        if not ((s1 < s2 and e1 < e2) or (s1 > s2 and e1 > e2)):
            return True
    return False

gaps = []
for length in reversed(lengths):
    gap = gen(length + 2*flank_length)
    while overlap(gaps, gap):
        gap = gen(length)
    gaps += [gap]
    # print(gaps)

gaps = sorted(gaps, key = lambda g: g[0])

# Write BED file
with open(inserts_file, 'w') as f:
    for start, end in gaps:
        f.write('%s\t%i\t%i\n' % (contig, start+flank_length, end-flank_length))

# Write BED file with flanks
with open(gaps_file, 'w') as f:
    for start, end in gaps:
        f.write('%s\t%i\t%i\n' % (contig, start, end))

# Write insertion breakpoints in BED format
with open(breakpoints_file, 'w') as f:
    accum = 0
    for start, end in gaps:
        f.write('%s\t%i\t%i\n' % (contig, start-accum,
            start-accum+(2*flank_length)))
        accum += (end - start) - 2*flank_length
