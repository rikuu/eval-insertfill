import sys

with open(sys.argv[2], 'w') as fo:
  with open(sys.argv[1], 'r') as f:
    buffer = ''
    for line in f:
      if line[0] == '>':
        if len(buffer) > 0:
          fo.write(buffer + '\n')
          buffer = ''
        fo.write(line)

      if len(buffer) >= 60:
        if len(buffer) > 0:
          fo.write(buffer[:60] + '\n')
          buffer = buffer[60:]

      if line[0] != '>':
        buffer += line.replace('\n', '').replace('N', '').replace('n', '')

    fo.write(buffer)
