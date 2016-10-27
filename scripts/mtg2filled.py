import sys

k = 31
fuz = 10

if len(sys.argv) != 3:
    print('Usage:', sys.argv[0], '<gaps.fa> <mtg.insertions.fa>')
    sys.exit(1)

def meta_parse(file, parser):
    with open(file, 'r') as f:
        seq = []
        for line in f:
            if line[0] == '>' and len(seq) > 1:
                parser(seq)
                seq = []
            seq += [line.rstrip()]
        parser(seq)

left_flank, right_flank = {}, {}

def parse_flanks(seq):
    split = seq[0].split(':')[1].split('-')
    pos = int(split[0]) + k + fuz

    gap = ''.join(seq[1:])

    left_flank[pos] = gap[:k+fuz]
    right_flank[pos] = gap[-(k+fuz):]

def parse_fills(seq):
    split = seq[0].split('_')
    contig = split[1]
    pos = int(split[split.index('pos')+1])

    fill = ''.join(seq[1:])
    fill = left_flank[pos] + fill + right_flank[pos]

    start = int(pos) - (k+fuz)
    end = start + len(fill)
    print('>%s:%i-%i\n%s' % (contig, start, end, fill))

meta_parse(sys.argv[1], parse_flanks)
meta_parse(sys.argv[2], parse_fills)
