import sys, copy
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt, log
from collections import defaultdict
from scipy.interpolate import spline
import numpy as np

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
              'text.latex.preamble': ['\usepackage{gensymb}'],
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
    ax.set_xlabel("Gap Length (log)")
    ax.set_ylabel("Normalized edit distance")

    ax.set_ylim([0.0, 1.0])

    return ax

# Read scores
#        similarity
# normal    0
# filter    1
dds = [defaultdict(lambda: defaultdict(list)) for i in range(2)]

with open(sys.argv[1], 'r') as f:
  for l in f:
    # LENGTH MEAN STDDEV FILTER NORMAL
    d = l.rstrip().split()

    length = int(d[0])

    #mean = int(d[1])
    #stddev = int(d[2])

    norm = lambda i: 1. - (float(d[i]) / length)

    dds[0][150][length].append(norm(1))
    dds[1][150][length].append(norm(2))
    dds[0][1500][length].append(norm(3))
    dds[1][1500][length].append(norm(4))
    dds[0][3000][length].append(norm(5))
    dds[1][3000][length].append(norm(6))

    # NOTE: 9000 = all reads
    dds[0][9000][length].append(norm(7))
    dds[1][9000][length].append(norm(8))

def plot_fills(ax, normal, filter, legend=True):
    lengths = sorted(filter.keys())
    #smooth_lengths = np.logspace(log(min(lengths), 10), log(max(lengths), 10), steps)

    average = lambda l: sum(l) / float(len(l))
    stripe = lambda d, f: [f(d[i]) for i in lengths]
    #smooth = lambda d, f: spline(lengths, stripe(d, f), smooth_lengths)

    ax.plot(lengths, stripe(filter, average), '-',
            lengths, stripe(normal, average), '--')

    if legend:
        ax.legend(['Filter', 'All reads'])

latexify(fig_width=6.9*4, columns=3, rows=1)
fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
for ax, l in zip([ax1, ax2, ax3, ax4], [150, 1500, 3000, 9000]):
    plot_fills(format_axes(ax), dds[0][l], dds[1][l])
plt.tight_layout()
plt.show()

# latexify()
# for l in [150, 1500, 3000]:
#     fig = plt.figure()
#     plot_fills(format_axes(fig.add_subplot(111)), dds[0][l], dds[1][l], l == 150)
#     plt.tight_layout()
#     if len(sys.argv) >= 2:
#         plt.savefig(sys.argv[2] + "." + str(l) + ".pgf")
#     fig.clf()
