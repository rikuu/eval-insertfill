import sys

k = 31
fuz = 10

def parse(file):
    buffer = []

    with open(file, 'r') as f:
        for line in f:
            if line[0] != '>':
                buffer += list(line.rstrip())

    return buffer[k+fuz : -(k+fuz)]

known = parse(sys.argv[1])
filled = parse(sys.argv[2])

wrong = 0
for i, j in zip(known, filled):
    if i != j:
        wrong += 1

print 1. - (float(wrong) / len(known))
