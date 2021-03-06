import sys, copy
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt, log
from collections import defaultdict
import numpy as np

if len(sys.argv) < 4:
    print('Usage: %s <results> <"table","tools","indel"> <sampling rate> [.pgf]' % sys.argv[0])
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
    # ax.set_yscale('log')
    ax.set_xlabel("Insertion Length (log)")
    ax.set_ylabel("Score")

    ax.set_ylim([0.0, 1.01])

    return ax

dds = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

with open(sys.argv[1], 'r') as f:
    for l in f:
        d = l.rstrip().split()
        length = int(d[0]) - 82

        def append(list, i, f=min, skip_unfilled=False):
            value = d[i]

            if '*' in value:
                if skip_unfilled:
                    return
                else:
                    value = value[:-1]

            if ',' in value:
                l = [float(j) for j in value.split(',')]
                value = f(l)
            else:
                value = float(value)

            list.append(value / (length + 82))

        if sys.argv[2] == 'indel':
            append(dds[0][150][length], 1)
            append(dds[1][150][length], 2)
            append(dds[0][1500][length], 3)
            append(dds[1][1500][length], 4)
            append(dds[0][3000][length], 5)
            append(dds[1][3000][length], 6)

            # NOTE: 9000 = all reads
            append(dds[0][9000][length], 7)
            append(dds[1][9000][length], 8)
        else:
            for i in range(len(d)-1):
                append(dds[i][9000][length], i+1)
def avg(l):
    assert(len(l) != 0)
    return sum(l) / len(l)

def median(l):
    if len(l) == 0: return 0
    if len(l) == 1: return l[0]

    s = sorted(l)
    m = int(len(l) / 2)
    if len(l) % 2 == 1:
        return s[m]
    else:
        return avg(s[m-1:m+1])

def dictsum(*args):
    for arg in args:
        assert(args[0].keys() == arg.keys())

    s = {}
    for k in args[0].keys():
        s[k] = [v for v in arg[k] for arg in args]
    return s

def plot_fills(ax, plots, steps=50, legend=True):
    lengths = sorted(plots[0][0].keys())
    log_lengths = np.logspace(log(min(lengths), 10), log(max(lengths), 10), steps).tolist()
    smooth_lengths = [(i+j) / 2 for i, j in zip(log_lengths[:-1], log_lengths[1:])]

    between = lambda d, f, i, j: [f(d[x]) for x in lengths if x >= i and x < j]
    smooth = lambda d, f: [f(between(d, f, i, j)) for i, j in zip(log_lengths[:-1], log_lengths[1:])]

    for i, plot in enumerate(plots):
        ax.plot(smooth_lengths, smooth(plot[0], avg), plot[2], label=plot[1], color=tableau20[2*i])

    if legend:
        ax.legend(loc='best')

##### Print table
if sys.argv[2] == 'table':
    avg = lambda data: sum(data) / len(data)
    part = lambda data, low, high: [x for x in data if x >= low and x <= high]
    count = lambda data, indices: avg([avg(data[index]) for index in indices])
    avg_in = lambda data, low, high: count(data, part(data.keys(), low, high))
    all_tools = lambda low, high: [round(avg_in(dds[i][9000], low, high), 3) for i in [2, 0, 1]]
    splits = [0, 100, 300, 500, 1000, 10000]
    for low, high in zip(splits[:-1], splits[1:]):
        print(low, '--', high, all_tools(low, high))
    #print([all_tools(low, high) for low, high in [(0, 100), (100, 500), (500, 1000), (1000, 5000), (5000, 10000)])
    sys.exit(0)

# latexify(fig_width=6.9*1.5, columns=1, rows=1)
latexify(rows=1.3, columns=1.5)
fig = plt.figure()
ax = fig.add_subplot(111)
plots = []

if sys.argv[2] == 'tools':
    plots = [(dds[2][9000], 'Pindel', '-'),
        (dds[4][9000], 'GapFiller', '-'),
        (dds[6][9000], 'Sealer', '-'),
        (dds[3][9000], 'MindTheGap', '-'),
        (dds[0][9000], 'Gap2Seq', '-'),
        (dds[5][9000], 'GapCloser', '-'),
        (dds[1][9000], 'Gap2Seq\n+ filter', '-')]
elif sys.argv[2] == 'indel':
    plots = [(dictsum(dds[0][150], dds[0][1500], dds[0][3000], dds[0][9000]), 'All reads', '--'),
        (dds[1][150], 'Filter (150)', '-'),
        (dds[1][1500], 'Filter (1500)', '-'),
        (dds[1][3000], 'Filter (3000)', '-')]
        #(dds[1][9000], 'Filter (All)', '-')]
elif sys.argv[2] == 'validated':
    plots = [(dds[6][9000], 'Pindel', '-'),
        (dds[3][9000], 'GapFiller', '-'),
        (dds[5][9000], 'Sealer', '-'),
        (dds[2][9000], 'MindTheGap', '-'),
        (dds[0][9000], 'Gap2Seq', '-'),
        (dds[4][9000], 'GapCloser', '-'),
        (dds[1][9000], 'Gap2Seq\n+ filter', '-')]
elif sys.argv[2] == 'biological':
    plots = [(dds[2][9000], 'GapFiller', '-'),
        (dds[3][9000], 'Sealer', '-')
        (dds[1][9000], 'MindTheGap', '-'),
        (dds[0][9000], 'Gap2Seq', '-')]

plot_fills(format_axes(ax), plots, sys.argv[3])
plt.tight_layout()
if len(sys.argv) == 5:
    plt.savefig(sys.argv[4], dpi=300)
else:
    plt.show()
