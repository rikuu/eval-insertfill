import sys
from random import choice

alphabet = ['a', 'c', 'g', 't']

with open(sys.argv[2], 'w') as fo:
  with open(sys.argv[1], 'r') as f:
    for line in f:
      buffer = line
      if buffer[0] != '>':
        buffer = list(buffer)
        for i in range(len(buffer)):
          if buffer[i] == 'N' or buffer[i] == 'n':
            buffer[i] = choice(alphabet)
        buffer = ''.join(buffer)
      fo.write(buffer)
