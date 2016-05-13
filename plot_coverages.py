import sys
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt
from collections import defaultdict

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

def latexify(fig_width=None, fig_height=None, columns=1):
    if fig_width is None:
        fig_width = 3.39 if columns==1 else 6.9 # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height +
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    params = {'backend': 'ps',
              'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize': 8, # fontsize for x and y labels (was 10)
              'axes.titlesize': 8,
              'text.fontsize': 8, # was 10
              'legend.fontsize': 8, # was 10
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)

def format_axes(ax):
    # Remove frame lines
    ax.spines["top"].set_visible(False)
    # ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # ax.spines["left"].set_visible(False)

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

    return ax

# latexify(fig_width=6.9/1.5)
fig = plt.figure()
ax = format_axes(fig.add_subplot(111))

overlap_over, unmapped_over, filter_over = [], [], []
with open(sys.argv[1], 'r') as f:
  for l in f:
    d = l.rstrip().split()

    mean = int(d[0])
    length = int(d[1])

    if float(d[2]) > 25:
        # ax.scatter(mean, length, color=tableau20[1], alpha=0.4)
        overlap_over.append((mean, length))

    if float(d[3]) > 25:
        # ax.scatter(mean, length, color=tableau20[4], alpha=0.4)
        unmapped_over.append((mean, length))

    if float(d[4]) > 25:
        ax.scatter(mean, length, color=tableau20[3], alpha=0.8)
        filter_over.append((mean, length))

ax.set_xlim([100,5000])
ax.set_ylim([100,5000])

# overlap_over = sorted(overlap_over, key=lambda x: x[1])
# unmapped_over = sorted(unmapped_over, key=lambda x: x[1])
# filter_over = sorted(filter_over, key=lambda x: x[1])
#
# means = lambda x: x[0]
# lengths = lambda x: x[1]
#
# ax.plot(map(lengths, overlap_over), map(means, overlap_over), color=tableau20[0])
# ax.plot(map(lengths, unmapped_over), map(means, unmapped_over), color=tableau20[4])
# ax.plot(map(lengths, filter_over), map(means, filter_over), color=tableau20[8])

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("Gap Length")
ax.set_ylabel("Mean Insert Size")

plt.tight_layout()
plt.show()
# plt.savefig(label + ".pdf")
# plt.savefig(label + ".pgf")
