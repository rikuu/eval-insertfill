import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

if len(sys.argv) < 2:
    print('Usage: %s <times> [.pgf]' % sys.argv[0])
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
              'axes.labelsize': 6,
              'axes.titlesize': 6,
              'legend.fontsize': 6,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'figure.figsize': [fig_width, fig_height],
              'font.family': 'serif',
              'font.size': 8
    }

    matplotlib.rcParams.update(params)

def format_axes(ax, log=False):
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

    if log: ax.set_yscale('log')

    return ax

def plot(ax, plots, log=False, ticks=True):
    width = 0.15
    for i, plot in enumerate(plots):
        ax.bar(i*width*2, plot[0], width, color=tableau20[2*i], log=log)

    if ticks:
        ax.set_xticks([width/1.9 + i*2*width for i in range(len(plots))])
        ax.set_xticklabels([plot[1] for plot in plots])

total = 0
fills = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        if 'insertions' in line:
            total = int(line.split()[0])

        if not 'S:' in line: continue
        line = line.rstrip().split()
        fills.append(int(line[line.index('S:')+1]))# + int(line[line.index('F:')+1]))

total = total / 100
plots = [(fills[6] / total, 'Pindel'),
    (fills[3] / total, 'GapFiller'),
    (fills[5] / total, 'Sealer'),
    (fills[2] / total, 'MindTheGap'),
    (fills[0] / total, 'Gap2Seq'),
    (fills[4] / total, 'GapCloser'),
    (fills[1] / total, 'Gap2Seq\n+ filter')]

# plots = [(0, 'Pindel'),
#     (300 / 3, 'GapFiller'),
#     (172 / 3, 'Sealer'),
#     (139 / 3, 'MindTheGap'),
#     (171 / 3, 'Gap2Seq'),
#     (281 / 3, 'GapCloser'),
#     (249 / 3, 'Gap2Seq\n+ filter')]

latexify(columns=1.5)
fig, ax = plt.subplots()
ax = format_axes(ax)

ax.set_ylabel('Insertions (\%)')
plot(ax, plots)
plt.tight_layout()
if len(sys.argv) == 3:
    plt.savefig(sys.argv[2], dpi=300)
else:
    plt.show()
