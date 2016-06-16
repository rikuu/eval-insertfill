import sys
import numpy as np

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
    lines = {}
    identifier = ''
    buf = ''

    with open(file, 'r') as f:
        for line in f:
            if line[0] != '>':
                buf += line.rstrip().upper()
            else:
                if len(buf) != 0:
                    lines[identifier] = buf
                    buf = ''
                identifier = line.rstrip()[1:]

                # Parse masking format - Might not work at some point
                s = identifier[6:].split('-')
                identifier = s[1] - s[0]
    lines[identifier] = buf

    return lines

known = parse(sys.argv[1])

results = {}
for i in known.keys(): results[i] = []

for f in sys.argv[2:]:
    filled = parse(f)

    if sorted(known.keys()) != sorted(filled.keys()):
        print(f, '\n', sorted(known.keys()), '\n', sorted(filled.keys()))
        sys.exit(1)

    for i in known.keys():
        results[i] += [str(edit_distance(known[i], filled[i]))]

for i in sorted(known.keys()):
    print(i, ' '.join(results[i]))
