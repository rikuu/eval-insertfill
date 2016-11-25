import sys

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == '#': continue
        fields = line.rstrip().split('\t')

        id = 'chr' + fields[0]
        start = int(fields[3])

        attributes = {key: value for (key, value) in [attr.split('=') for attr in fields[8].split(';')]}

        # Parsing based on the format used by Vergara et al.
        assert('Note' in attributes)
        notes = {key: value for (key, value) in [attr.split(':') for attr in attributes['Note'].split(',')]}

        assert('length' in notes and 'seq' in notes)
        end = start + int(notes['length'])
        sequence = notes['seq']

        print('%s\t%i\t.\t%s\t%s\t.\tPASS\t' % (id, start, sequence[0], sequence))
