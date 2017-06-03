import sys

# Gap2Seq defaults
fuz = 10
k = 31

if len(sys.argv) not in [3, 4, 5, 6]:
    print('Usage:', sys.argv[0], '<reference.fa> <variants.vcf> [gaps.fa] [filled.fa] [gapped_reference.fa]')
    sys.exit(1)

REFERENCE = sys.argv[1]
VARIANTS = sys.argv[2]
GAPS = sys.argv[3] if len(sys.argv) >= 4 else None
FILLED = sys.argv[4] if len(sys.argv) >= 5 else None
GAPPED_REFERENCE = sys.argv[5] if len(sys.argv) == 6 else None

# Parse chromosomes from the reference into a dictionary
reference = {}
with open(REFERENCE, 'r') as f:
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

gaps = open(GAPS, 'w') if GAPS != None else None
filled = open(FILLED, 'w') if FILLED != None else None

gapped_reference = reference if GAPPED_REFERENCE != None else None

with open(VARIANTS, 'r') as f:
    for line in f:
        if line[0] == '#': continue
        fields = line.rstrip().split('\t')

        # Parse VCF fields
        insert = fields[4][1:]
        comment, start, end = fields[0], int(fields[1]) - 1, \
            (int(fields[1]) - 1) + (len(insert) - 1)

        # Validated inserts and reference genome use different names?
        if comment[:3] == 'chr': comment = comment[3:]
        comment = 'CHROMOSOME_' + comment

        # Extract kmers from reference seqeuences
        left_end = min(start+k+fuz, len(reference[comment]))
        right_start = max(end-(k+fuz), 0)
        left = reference[comment][start:left_end]
        right = reference[comment][right_start:end]
        gap = left + 'N' * len(insert) + right

        if gaps == None:
            print('>%s:%i-%i\n%s' % (comment, start, end, gap))
        else:
            gaps.write('>%s:%i-%i\n%s\n' % (comment, start, end, gap))

        if filled != None:
            fill = left + insert + right
            filled.write('>%s:%i-%i\n%s\n' % (comment, start, end, fill))

        if gapped_reference != None:
            gapped_reference[comment] = gapped_reference[comment][:start] + \
                'N' * len(insert) + gapped_reference[comment][start+1:]

if GAPPED_REFERENCE != None:
    with open(GAPPED_REFERENCE, 'w') as f:
        for comment in gapped_reference.keys():
            f.write('>%s\n%s\n' % (comment, gapped_reference[comment]))

if gaps != None: gaps.close()
if filled != None: filled.close()
