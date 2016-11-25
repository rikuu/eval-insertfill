import sys, copy, re
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

def format_axes(ax, log=True):
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

    # ax.set_ylim([0.0, 1.0])

    return ax

def plot(ax, results, value, log=True, bwa=None):
    tools = [('GAPCLOSER', 'GapCloser'), ('GAPFILLER_BWA', 'GapFiller'),
        ('GAP2SEQ', 'Gap2Seq'), ('GAP2SEQ21', 'Gap2Seq (with filter)')]
    for i, tool in enumerate(tools):
        ax.bar(ind + i*width, [float(results[tool[0]][key][value]) for key in results[tool[0]].keys()],
            width, color=tableau20[2*i], log=log, label=tool[1])

    if bwa != None:
        ax.bar(ind + (len(tools)-1)*width, [float(bwa[key][value-1]) for key in results[tool[0]].keys()],
            width, color=tableau20[8], log=log, label='BWA')

    ax.legend()

# assembly & tool & user time & wall clock time & peak memory (Gb)
text = \
"""staph & ABySS & GAP2SEQ21 & 0.23 & 0.58 & 0.011436
staph & ABySS2 & GAP2SEQ21 & 0.18 & 0.35 & 0.011344
staph & Allpaths-LG & GAP2SEQ21 & 0.21 & 0.41 & 0.011388
staph & Bambus2 & GAP2SEQ21 & 0.22 & 0.45 & 0.01152
staph & MSR-CA & GAP2SEQ21 & 0.23 & 0.42 & 0.011456
staph & SGA & GAP2SEQ21 & 0.85 & 2.06 & 0.013016
staph & SOAPdenovo & GAP2SEQ21 & 0.14 & 0.34 & 0.011276
staph & Velvet & GAP2SEQ21 & 0.24 & 0.51 & 0.011588
staph & TOTAL & GAP2SEQ21 & 2.3 & 5.1 & 0.1
AVERAGE & TOTAL & GAP2SEQ21 & 0.3 & 0.6 & 0.0
staph & ABySS & GAP2SEQ & 95.11 & 27.14 & 0.483876
staph & ABySS2 & GAP2SEQ & 82.15 & 22.83 & 0.471852
staph & Allpaths-LG & GAP2SEQ & 84.5 & 23.82 & 0.481816
staph & Bambus2 & GAP2SEQ & 99.15 & 34.15 & 0.445292
staph & MSR-CA & GAP2SEQ & 99.76 & 34.2 & 0.441872
staph & SGA & GAP2SEQ & 347.27 & 167.56 & 0.461204
staph & SOAPdenovo & GAP2SEQ & 75.22 & 19.96 & 0.4742
staph & Velvet & GAP2SEQ & 117.63 & 41.4 & 0.474644
staph & TOTAL & GAP2SEQ & 1000.8 & 371.1 & 3.7
AVERAGE & TOTAL & GAP2SEQ & 125.1 & 46.4 & 0.5
staph & ABySS & GAPCLOSER & 91.23 & 26.8 & 0.702776
staph & ABySS2 & GAPCLOSER & 54.31 & 28.57 & 0.70278
staph & Allpaths-LG & GAPCLOSER & 209.56 & 137.57 & 0.716308
staph & Bambus2 & GAPCLOSER & 39.98 & 27.79 & 0.746912
staph & MSR-CA & GAPCLOSER & 38.05 & 29.57 & 0.7385
staph & SGA & GAPCLOSER & 52.56 & 25.12 & 0.7411
staph & SOAPdenovo & GAPCLOSER & 52.81 & 30.98 & 0.70278
staph & Velvet & GAPCLOSER & 39.49 & 25.33 & 0.729824
staph & TOTAL & GAPCLOSER & 578.0 & 331.7 & 5.8
AVERAGE & TOTAL & GAPCLOSER & 72.2 & 41.5 & 0.7
staph & ABySS & GAPFILLER_BWA & 739.31 & 339.06 & 0.118116
staph & ABySS2 & GAPFILLER_BWA & 693.45 & 309.67 & 0.118184
staph & Allpaths-LG & GAPFILLER_BWA & 708.37 & 319.95 & 0.1181
staph & Bambus2 & GAPFILLER_BWA & 854.76 & 448.7 & 0.118156
staph & MSR-CA & GAPFILLER_BWA & 792.37 & 371.78 & 0.132756
staph & SGA & GAPFILLER_BWA & 1522.84 & 861.73 & 0.194892
staph & SOAPdenovo & GAPFILLER_BWA & 600.06 & 253.28 & 0.118016
staph & Velvet & GAPFILLER_BWA & 851.35 & 417.09000000000003 & 0.118176
staph & TOTAL & GAPFILLER_BWA & 6762.5 & 3321.3 & 1.0
AVERAGE & TOTAL & GAPFILLER_BWA & 845.3 & 415.2 & 0.1"""

# Separately run alignment time and memory stats
bwa = {'ABySS': ['2:56.90', '2007081'],
    'ABySS2': ['2:20.10', '1783194'],
    'Allpaths-LG': ['2:00.31', '1797091'],
    'Bambus2': ['1:59.51', '1789254'],
    'MSR-CA': ['2:01.91', '1822535'],
    'SGA': ['2:09.15', '1753986'],
    'SOAPdenovo': ['2:12.91', '1759993'],
    'Velvet': ['2:07.47', '1814879']}

# Fix formatting for BWA stats
for key in bwa.keys():
    # Format times as minutes
    time = bwa[key][0].split(':')
    minutes = int(time[0])
    seconds = float(time[1])
    bwa[key][0] = minutes + (seconds / 60.)

    # Format memory usage as Gb
    bwa[key][1] = int(bwa[key][1]) / (1000*1000)

results = defaultdict(lambda: defaultdict(list))
for line in text.split('\n'):
    data = line.split(' & ')
    if data[0] == 'AVERAGE' or data[1] == 'TOTAL': continue
    results[data[2]][data[1]] = data[3:]

    # Add BWA stats to GAP2SEQ21
    if data[2] == 'GAP2SEQ21':
        time = float(results[data[2]][data[1]][1])
        results[data[2]][data[1]][1] = str(time + bwa[data[1]][0])
        #mem = float(results[data[2]][data[1]][2])
        #results[data[2]][data[1]][2] = str(mem + bwa[data[1]][1])

print(dict(results))
latexify(columns=1.7)

ind = np.arange(len(results['GAP2SEQ'].keys()))
width = 0.2

fig, ax = plt.subplots()
ax = format_axes(ax)

ax.set_xticks(ind + 2*width)
ax.set_xticklabels(results['GAP2SEQ'].keys())

# Time plot
ax.set_ylabel('Runtime (min, log)')
ax.set_ylim([0.15, 1000.0])
plot(ax, results, 1, bwa=bwa)
plt.tight_layout()
plt.savefig(sys.argv[1] + "_time.pgf")
fig.clf()

fig, ax = plt.subplots()
ax = format_axes(ax, log=False)
ax.set_xticks(ind + 2*width)
ax.set_xticklabels(results['GAP2SEQ'].keys())

# Memory plot
ax.set_ylabel('Memory (GB)')
ax.set_ylim([0.0, 1.0])
plot(ax, results, 2, log=False)
plt.tight_layout()
plt.savefig(sys.argv[1] + "_memory.pgf")
fig.clf()
