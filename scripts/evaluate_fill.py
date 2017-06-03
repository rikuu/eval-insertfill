import sys, re
import numpy as np
from collections import defaultdict

DEBUG = False

# Defaults from Gap2Seq
k = 31
fuz = 10

if len(sys.argv) < 3:
    print('Usage:', sys.argv[0], '<known.fasta> <evaluate.fasta> [<evaluate2.fasta> ...]')
    sys.exit(1)

# Computes Levenshtein distance
def edit_distance(a, b):
    if len(a) < len(b): return edit_distance(b, a)
    if len(b) == 0: return len(a)

    a, b = np.array(tuple(a)), np.array(tuple(b))

    v1 = np.arange(b.size + 1)
    for s in a:
        v2 = v1 + 1
        v2[1:] = np.minimum(np.minimum(v2[1:], v2[0:-1] + 1),
            np.add(v1[:-1], b != s))
        v1 = v2

    return v1[-1]

# Parses a fasta file of gaps into a dictionary
def parse(file):
    lines, identifier, buf = defaultdict(list), '', ''
    with open(file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line[0] != '>':
                buf += line.upper()
            else:
                if len(buf) != 0:
                    lines[identifier].append(buf)
                    buf = ''
                identifier = line[1:]
    lines[identifier].append(buf)
    return lines

known = parse(sys.argv[1])
for v in known.values(): assert len(v) == 1

results = defaultdict(list)
for f in sys.argv[2:]:
    filled = parse(f)

    if DEBUG and sorted(known.keys()) != sorted(filled.keys()):
        print(f, '\n', sorted(known.keys()), '\n', sorted(filled.keys()))
        sys.exit(1)

    for i in known.keys():
        if not i in filled:
            results[i].append(str(len(known[i][0])) + '*')
        else:
            distances = []
            for fill in filled[i]:
                dist = ''
                if fill.upper().count('N') == len(fill) - (2*(k+fuz)):
                    dist = str(len(known[i][0])) + '*'
                else:
                    dist = str(edit_distance(known[i][0], fill))
                distances.append(dist)
            results[i].append(','.join(distances))

# Print plottable lengths
for i in sorted(known.keys(), key = lambda key: len(known[key][0])):
    print(len(known[i][0]), ' '.join(results[i]))
