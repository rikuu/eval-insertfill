import sys

def parse_ids(file):
  ids = []
  with open(file, 'r') as f:
    for l in f:
      if l[0] == '>':
        if l[-2] == '/':
          ids.append(l[1:-2])
        else:
          ids.append(l[1:])
  return set(ids)

def recall(tp, fn):
    if tp == 0:
        return 0.
    return tp / float(tp + fn)

def precision(tp, fp):
    if tp == 0:
        return 0.
    return tp / float(tp + fp)

def fscore(tp, fn, fp):
    if tp == 0:
        return 0.
    return (2*tp) / float(2*tp + fn + fp)

all_ids = parse_ids(sys.argv[1])
known_ids = parse_ids(sys.argv[2])
predicted_ids = parse_ids(sys.argv[3])

# true positives, reads correctly filtered in
tp = len(predicted_ids.intersection(known_ids))

# false positives, reads incorrectly filtered in
fp = len(predicted_ids - known_ids)

# true negatives, reads correcly filtered out
tn = len((all_ids - predicted_ids).intersection((all_ids - known_ids)))

# false negatives, reads incorrectly filtered out
fn = len((all_ids - predicted_ids) - (all_ids - known_ids))

# Recall, Precision, F-score
print recall(tp, fn), precision(tp, fp), fscore(tp, fn, fp)
