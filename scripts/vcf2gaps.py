import sys

# Gap2Seq defaults
fuz = 10
k = 31

if len(sys.argv) != 3:
    print('Usage:', sys.argv[0], '<reference.fa> <variants.vcf>')
    sys.exit(1)

# Parse chromosomes from the reference into a dictionary
reference = {}
with open(sys.argv[1], 'r') as f:
    comment, seq = '', ''
    for line in f:
        if line[0] == '>':
            if seq != '':
                reference[comment] = seq
                seq = ''
            comment = line.rstrip()[1:]
        else:
            seq += line.rstrip()
    reference[comment] = seq

with open(sys.argv[2], 'r') as f:
    for line in f:
        if line[0] == '#': continue
        fields = line.rstrip().split('\t')

        # Parse VCF fields
        insert = fields[4][1:]
        comment, start, end = fields[0], int(fields[1]) - 1, \
            int(fields[1]) + len(insert) - 1

        # Extract kmers from reference seqeuences
        left = reference[comment][start:min(start+k+fuz, len(reference[comment]))]
        right = reference[comment][max(end-(k+fuz), 0):end]
        gap = left + 'N' * len(insert) + right

        print('>%s:%i-%i\n%s' % (comment, start, end, gap))
