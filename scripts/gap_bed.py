import sys
import numpy as np
from math import log

if len(sys.argv) != 3:
    print('Usage:', sys.argv[0], '<inserts.bed> <gaps.bed>')
    sys.exit(1)

# TODO: get these from input
genome_length = 83257441
contig = 'chr17'

# k+fuz
flank_length = 41

num_gaps = 20
length_start = 11
length_end = 5000

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

# Write BED file
with open(sys.argv[1], 'w') as f:
  for start, end in gaps:
    f.write('%s\t%i\t%i\n' % (contig, start+flank_length, end-flank_length))

# Write BED file with flanks
with open(sys.argv[2], 'w') as f:
  for start, end in gaps:
    f.write('%s\t%i\t%i\n' % (contig, start, end))
