import sys

if len(sys.argv) < 4:
    print('Usage:', sys.argv[0], '<all_reads.fa> <known_reads.fa>',
        '<predicted_reads.fa> [<predicted_reads2.fa> ..]')

def parse_ids(file):
  ids = []
  with open(file, 'r') as f:
    for l in f:
      if l[0] == '>':
        if len(l) < 3: continue
        ids.append(l[1:])
  return set(ids)

recall = lambda tp, fn: tp / (tp + fn)
precision = lambda tp, fp: tp / (tp + fp)
fscore = lambda tp, fn, fp: (2*tp) / (2*tp + fn + fp)

all_ids = parse_ids(sys.argv[1])
known_ids = parse_ids(sys.argv[2])
known_false = all_ids - known_ids

scores = ''
for predictor in sys.argv[3:]:
    predicted_ids = parse_ids(predictor)
    predicted_false = all_ids - predicted_ids

    # true positives, reads correctly filtered in
    tp = float(len(predicted_ids.intersection(known_ids)))

    if tp == 0.:
        scores += '0 0 0 '
        continue

    # false positives, reads incorrectly filtered in
    fp = float(len(predicted_ids - known_ids))

    # true negatives, reads correcly filtered out
    tn = float(len(predicted_false.intersection(known_false)))

    # false negatives, reads incorrectly filtered out
    fn = float(len(predicted_false - known_false))

    # Recall, Precision, F-score
    scores += '%.18f %.18f %.18f ' % \
        (recall(tp, fn), precision(tp, fp), fscore(tp, fn, fp))

print(scores)
