import sys

start = int(sys.argv[1])
end = int(sys.argv[2])
length = end - start

seq_length = 0
with open(sys.argv[3], 'r') as f:
  for line in f:
    if line[0] != '>':
      seq_length += len(line.rstrip())

print float(seq_length) / float(length)
