import sys, itertools
from random import randint

fasta = sys.argv[1]

gap_lengths = {}
for c, l in zip(sys.argv[2].split(','), sys.argv[3].split(',')):
    gap_lengths[c] = int(l)

with open(fasta, 'r') as f:
  for h, l in itertools.izip_longest(*[f]*2):
    h = h.rstrip()
    l = l.rstrip()

    contig = h.split('|')[0][1:]
    gap = gap_lengths[contig]

    print h
    if h[0] == '>' and contig in gap_lengths:
      print len(l) - 200 - gap
      start = randint(200, len(l) - 200 - gap)
      print l[:start] + 'n'*gap + l[start+gap+1:]
    else:
      print l
