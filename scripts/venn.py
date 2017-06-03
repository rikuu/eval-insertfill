import sys, re

if len(sys.argv) < 2:
    print('Usage: %s <vcf> <filled> [filled2] ..' % sys.argv[0])
    sys.exit(1)

def parse_vcf(file):
    lengths = {}
    with open(file, 'r') as f:
        for line in f:
            if line[0] == '#': continue
            fields = line.rstrip().split('\t')

            lengths[int(fields[1]) - 1] = len(fields[4]) - 1
    return lengths

def parse(file):
    fills = {}
    success, fail = [], []
    identifier, buf = '', ''

    with open(file, 'r') as f:
        for line in f:
            if line[0] != '>':
                buf += line.rstrip().upper()
            else:
                if len(buf) != 0:
                    if 'N' in buf or 'n' in buf:
                        fail.append(identifier)
                    else:
                        success.append(identifier)
                    fills[identifier] = buf
                    buf = ''

                # Identify insertion sites by breakpoint positions
                identifier = int(re.split(':|-', line.rstrip()[1:])[1])

    if 'N' in buf or 'n' in buf:
        fail.append(identifier)
    else:
        success.append(identifier)

    return set(success), set(fail), fills

estimated_lengths = parse_vcf(sys.argv[1])
print('%i insertions, length: (%i, %i, %2.f)' % (len(estimated_lengths.values()),
    min(estimated_lengths.values()), max(estimated_lengths.values()),
    sum(estimated_lengths.values()) / len(estimated_lengths.values())))

sets = [set() for i in range(len(sys.argv[2:]))]
for i, f in enumerate(sys.argv[2:]):
    success, fail, fills = parse(f)
    sets[i] = sets[i].union(success)

    diff = [len(fills[key]) - estimated_lengths[key] for key in fills.keys()]
    print('%s, S: %i F: %i D: (%i, %i, %2.f)' % \
        (f, len(success), len(fail), min(diff), max(diff), sum(diff) / len(diff))

# from matplotlib import pyplot as plt
# from matplotlib_venn import venn3, venn2
# assert(len(sets) == 3)
# # venn3((sets[0], sets[1], sets[2]), ('MindTheGap', 'Gap2Seq', 'Gap2Seq (with filter)'))
# venn2((sets[0], sets[2]), ('MindTheGap', 'Gap2Seq (with filter)'))
# plt.savefig('../../cpm/venn.pgf')
