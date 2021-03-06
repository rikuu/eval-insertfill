import sys, copy
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt, log
from collections import defaultdict
import numpy as np

if len(sys.argv) < 3:
    print('Usage: %s <results> <"recall","precision","fscore"> <sampling rate> [.pgf]' % sys.argv[0])
    sys.exit(1)

SPINE_COLOR = 'gray'

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale Tableau20 colors to [0, 1] range
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

def latexify(fig_width=None, fig_height=None, columns=1, rows=1):
    if fig_width is None:
        fig_width = columns*3.39

    if fig_height is None:
        golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = rows * (fig_width / columns) * golden_mean # height in inches

    # MAX_HEIGHT_INCHES = 8.0
    # if fig_height > MAX_HEIGHT_INCHES:
    #     print("WARNING: fig_height too large:" + fig_height +
    #           "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
    #     fig_height = MAX_HEIGHT_INCHES

    params = {'backend': 'ps',
              'text.latex.preamble': ['\\usepackage{gensymb}'],
              'axes.labelsize': 8, # fontsize for x and y labels (was 10)
              'axes.titlesize': 8,
              'text.fontsize': 8, # was 10
              'legend.fontsize': 8, # was 10
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'figure.figsize': [fig_width, fig_height],
              'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)

def format_axes(ax):
    # Remove frame lines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Ticks only on bottom and left
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    for spine in ['left', 'bottom']:
        ax.spines[spine].set_color(SPINE_COLOR)
        ax.spines[spine].set_linewidth(0.5)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_tick_params(direction='out', color=SPINE_COLOR)

    ax.set_xscale('log')
    ax.set_xlabel("Insertion Length (log)")

    ax.set_ylim([0.0, 1.01])

    return ax

def avg(l):
    assert(len(l) != 0)
    return float(sum(l)) / len(l)

def median(l):
    if len(l) == 0: return 0
    if len(l) == 1: return l[0]

    s = sorted(l)
    m = int(len(l) / 2)
    if len(l) % 2 == 1:
        return s[m]
    else:
        return avg(s[m-1:m+1])

def plot_between(ax, plots, steps=50, legend=True):
    lengths = sorted(plots[0][0].keys())
    log_lengths = np.logspace(log(min(lengths), 10), log(max(lengths), 10), steps).tolist()
    smooth_lengths = [(i+j) / 2 for i, j in zip(log_lengths[:-1], log_lengths[1:])]

    between = lambda d, f, i, j: [f(d[x]) for x in lengths if x >= i and x < j]
    smooth = lambda d, f: [f(between(d, f, i, j)) for i, j in zip(log_lengths[:-1], log_lengths[1:])]

    handles = [ax.plot(smooth_lengths, smooth(plot[0], avg), plot[2], label=plot[1], color=tableau20[plot[3]*2]) for plot in plots]
    if legend:
        ax.legend(handles=[handle[0] for i, handle in enumerate(handles) if plots[i][4]], loc='best')

# Read scores
#          Recall  Precision  F-score
# unmapped   0        1         2
# overlap    3        4         5
# gapfiller  6        7         8
# gapfiller2 9        10        11
# filter     12       13        14
# filter2    15       16        17
dds = [defaultdict(lambda: defaultdict(list)) for i in range(18)]

with open(sys.argv[1], 'r') as f:
  for l in f:
    # GAPLENGTH MEAN STDDEV OVERLAP UNMAPPED FILTER
    d = l.rstrip().split()

    if len(d) != 18:
        print(len(d))
        continue

    length = int(d[0])
    mean = int(d[1])
    stddev = int(d[2])

    for i in range(15):
        dds[i][mean][length].append(float(d[i+3]))

def dictsum(*args):
    for arg in args:
        assert(args[0].keys() == arg.keys())

    s = {}
    for k in args[0].keys():
        s[k] = [v for v in arg[k] for arg in args]
    return s

# latexify(columns=1.5)
# fig = plt.figure()
# ax = fig.add_subplot(111)
#
# plots = []
# if sys.argv[2] == "recall":
#     ax.set_ylabel("Recall")
#     # plots = [(dds[0][150], 'Unmapped', ':', 0, True),
#     #     (dds[3][150], 'Overlap', '--', 1, True),
#     #     (dds[6][150], 'GapFiller', '-', 2, True),
#     #     (dds[9][150], 'GapFiller (insertion)', '-', 3, True),
#     #     (dds[12][150], 'Filter', '--', 4, False),
#     #     (dds[15][150], 'Filter', '-', 4, True)]
#     plots = [(dictsum(dds[0][150], dds[0][1500], dds[0][3000]), 'Unmapped', ':', 0, True),
#         #(dictsum(dds[3][150], dds[3][1500], dds[3][3000]), 'Overlap', '--', 1, True),
#         (dictsum(dds[6][150], dds[6][1500], dds[6][3000]), 'GapFiller', '-', 2, True),
#         (dds[9][150], 'Filter (150)', '--', 3, False),
#         (dds[12][150], 'Filter (150)', '-', 3, True),
#         (dds[9][1500], 'Filter (1500)', '--', 4, False),
#         (dds[12][1500], 'Filter (1500)', '-', 4, True),
#         (dds[9][3000], 'Filter (3000)', '--', 5, False),
#         (dds[12][3000], 'Filter (3000)', '-', 5, True)]
# if sys.argv[2] == "precision":
#     ax.set_ylabel("Precision")
#     plots = [(dictsum(dds[1][150], dds[1][1500], dds[1][3000]), 'Unmapped', ':', 0, True),
#         #(dictsum(dds[4][150], dds[4][1500], dds[4][3000]), 'Overlap', '--', 1, True),
#         (dictsum(dds[7][150], dds[7][1500], dds[7][3000]), 'GapFiller', '-', 2, True),
#         (dds[10][150], 'Filter (150)', '--', 3, False),
#         (dds[13][150], 'Filter (150)', '-', 3, True),
#         (dds[10][1500], 'Filter (1500)', '--', 4, False),
#         (dds[13][1500], 'Filter (1500)', '-', 4, True),
#         (dds[10][3000], 'Filter (3000)', '--', 5, False),
#         (dds[13][3000], 'Filter (3000)', '-', 5, True)]
# if sys.argv[2] == "fscore":
#     ax.set_ylabel("F-score")
#     plots = [(dictsum(dds[2][150], dds[2][1500], dds[2][3000]), 'Unmapped', ':', 0, True),
#         #(dictsum(dds[5][150], dds[5][1500], dds[5][3000]), 'Overlap', '--', 1, True),
#         (dictsum(dds[8][150], dds[8][1500], dds[8][3000]), 'GapFiller', '-', 2, True),
#         (dds[11][150], 'Filter (150)', '--', 3, False),
#         (dds[14][150], 'Filter (150)', '-', 3, True),
#         (dds[11][1500], 'Filter (1500)', '--', 4, False),
#         (dds[14][1500], 'Filter (1500)', '-', 4, True),
#         (dds[11][3000], 'Filter (3000)', '--', 5, False),
#         (dds[14][3000], 'Filter (3000)', '-', 5, True)]

latexify(columns=1.5, rows=2)
f, (ax1, ax2) = plt.subplots(2, sharex=True)

ax1.set_ylabel("Recall")
plots = [(dictsum(dds[0][150], dds[0][1500], dds[0][3000]), 'Unmapped', ':', 0, True),
    (dictsum(dds[6][150], dds[6][1500], dds[6][3000]), 'GapFiller', '-', 2, True),
    (dds[9][150], 'Filter (150)', '--', 3, False),
    (dds[12][150], 'Filter (150)', '-', 3, True),
    (dds[9][1500], 'Filter (1500)', '--', 4, False),
    (dds[12][1500], 'Filter (1500)', '-', 4, True),
    (dds[9][3000], 'Filter (3000)', '--', 5, False),
    (dds[12][3000], 'Filter (3000)', '-', 5, True)]
plot_between(format_axes(ax1), plots, sys.argv[2], legend=False)

ax2.set_ylabel("Precision")
plots = [(dictsum(dds[1][150], dds[1][1500], dds[1][3000]), 'Unmapped', ':', 0, True),
    (dictsum(dds[7][150], dds[7][1500], dds[7][3000]), 'GapFiller', '-', 2, True),
    (dds[10][150], 'Filter (150)', '--', 3, False),
    (dds[13][150], 'Filter (150)', '-', 3, True),
    (dds[10][1500], 'Filter (1500)', '--', 4, False),
    (dds[13][1500], 'Filter (1500)', '-', 4, True),
    (dds[10][3000], 'Filter (3000)', '--', 5, False),
    (dds[13][3000], 'Filter (3000)', '-', 5, True)]
plot_between(format_axes(ax2), plots, sys.argv[2], legend=True)

plt.tight_layout()
if len(sys.argv) == 4:
    plt.savefig(sys.argv[3], dpi=300)
else:
    plt.show()
