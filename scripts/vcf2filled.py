import sys

# Defaults from Gap2Seq
fuz = 10
d_err = 500

if len(sys.argv) != 3:
    print('Usage:', sys.argv[0], '<inserts.bed> <variants.vcf>')
    sys.exit(1)

def approx_search(query, inserts):
    seq, start, end = query
    min_start, max_start = start - fuz, start
    min_end, max_end = end - d_err, end + d_err + fuz

    # insert = (sequence, start, end)
    for i in inserts:
        if (i[0] == seq) and \
                (i[1] >= min_start and i[1] <= max_start) and \
                (i[2] >= min_end and i[2] <= max_end):
            return '>%s:%i-%i' % (i[0], i[1], i[2])
    return ''

# Parse the BED file
inserts = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        fields = line.rstrip().split('\t')
        inserts += [(fields[0], int(fields[1]), int(fields[2]))]

# Find inserts in the VCF approximately matching the BED file
with open(sys.argv[2], 'r') as f:
    for line in f:
        if line[0] == '#': continue
        fields = line.rstrip().split('\t')

        insert = fields[4]
        seq, start, end = fields[0], int(fields[1]), int(fields[1]) + len(insert)

        header = approx_search((seq, start, end), inserts)
        if header != '':
            print(header+'\n'+insert)
