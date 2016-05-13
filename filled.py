import sys

with open(sys.argv[1], 'r') as f:
  for line in f:
    if line[0] != '>':
      if 'N' in line or 'n' in line:
        print '0'
        exit()
print 1
