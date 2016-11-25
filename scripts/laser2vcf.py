# Chromsome Pos Length Type Alt Qual
# CHROMOSOME_I 1943 3 INS CAT 32.7908

# Chromosome Pos ID Ref Alt Qual Filter Info
# CHROMOSOME_I 1943 . _ _CAT 32.7908 PASS

import sys

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
        fields = line.rstrip().split(' ')
        chromosome = fields[0]
        pos = int(fields[1])

        # Skip deletions for now
        if fields[3] != 'INS': continue

        ref = reference[chromosome][pos].upper()
        alt = ref + fields[4]

        qual = float(fields[5])
        print('%s\t%i\t.\t%s\t%s\t%f\tPASS' % (chromosome, pos, ref, alt, qual))
