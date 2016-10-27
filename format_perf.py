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

    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.set_xlabel("Gap Length (log)")
    # ax.set_ylabel("1 - Normalized Edit Distance")

    # ax.set_ylim([0.0, 1.0])

    return ax

# assembly & tool & user time & wall clock time & peak memory (Gb)
text = \
"""ABySS & GAP2SEQ & 95.11 & 27.14 & 0.483876
ABySS2 & GAP2SEQ & 82.15 & 22.83 & 0.471852
Allpaths-LG & GAP2SEQ & 84.5 & 23.82 & 0.481816
Bambus2 & GAP2SEQ & 99.15 & 34.15 & 0.445292
MSR-CA & GAP2SEQ & 99.76 & 34.2 & 0.441872
SGA & GAP2SEQ & 347.27 & 167.56 & 0.461204
SOAPdenovo & GAP2SEQ & 75.22 & 19.96 & 0.4742
Velvet & GAP2SEQ & 117.63 & 41.4 & 0.474644
TOTAL & GAP2SEQ & 1000.8 & 371.1 & 3.7
AVERAGE & GAP2SEQ & 125.1 & 46.4 & 0.5
ABySS & GAP2SEQ21 & 0.23 & 0.44 & 0.011428
ABySS2 & GAP2SEQ21 & 0.19 & 0.33 & 0.01134
Allpaths-LG & GAP2SEQ21 & 0.21 & 0.41 & 0.011388
Bambus2 & GAP2SEQ21 & 0.22 & 0.45 & 0.01152
MSR-CA & GAP2SEQ21 & 0.23 & 0.42 & 0.011456
SGA & GAP2SEQ21 & 0.85 & 2.06 & 0.013016
SOAPdenovo & GAP2SEQ21 & 0.14 & 0.34 & 0.011276
Velvet & GAP2SEQ21 & 0.24 & 0.51 & 0.011588
TOTAL & GAP2SEQ21 & 2.3 & 5.0 & 0.1
AVERAGE & GAP2SEQ21 & 0.3 & 0.6 & 0.0"""

results = [{}, {}]
for line in text.split('\n'):
    data = line.split(' & ')
    if data[0] in ['TOTAL', 'AVERAGE']: continue

    if data[1] == 'GAP2SEQ':
        results[0][data[0]] = data[2:]
    else:
        results[1][data[0]] = data[2:]

print(results)

ind = np.arange(len(results[0].keys()))
width = 0.2

fig, ax = plt.subplots()
ax = format_axes(ax)

ax.set_ylabel('Runtime (min)')
ax.set_xticks(ind + width)
ax.set_xticklabels(results[0].keys())

gap2seq20 = ax.bar(ind, [float(results[0][key][1]) for key in results[0].keys()], width, color=tableau20[0], log=True)
gap2seq21 = ax.bar(ind+width, [float(results[1][key][1]) for key in results[1].keys()], width, color=tableau20[2], log=True)
ax.legend((gap2seq20[0], gap2seq21[0]), ('Gap2Seq 2.0', 'Gap2Seq 2.1'))

plt.show()
