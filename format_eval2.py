text="""
ABySS & ORIGINAL & 5 & 10587 & 7935 & 31079 & 69 & 55885
ABySS2 & ORIGINAL & 5 & 10312 & 0 & 106796 & 35 & 9393
Allpaths-LG & ORIGINAL & 0 & 5991 & 0 & 110168 & 48 & 9900
Bambus2 & ORIGINAL & 0 & 24570 & 0 & 40233 & 99 & 29205
MSR-CA & ORIGINAL & 10 & 17276 & 0 & 64114 & 81 & 10353
SGA & ORIGINAL & 2 & 13811 & 0 & 9541 & 654 & 300607
SOAPdenovo & ORIGINAL & 2 & 35433 & 4055 & 69834 & 9 & 4857
Velvet & ORIGINAL & 25 & 24160 & 1270 & 46087 & 128 & 17688
TOTAL & ORIGINAL & 49.0 & 142140.0 & 13260.0 & 477852.0 & 1123.0 & 437888.0
AVERAGE & ORIGINAL & 6.1 & 17767.5 & 1657.5 & 59731.5 & 140.4 & 54736.0
ABySS & GAP2SEQ21 & 8 & 34696 & 4526 & 29538 & 1 & 3052
ABySS2 & GAP2SEQ21 & 5 & 19305 & 0 & 106796 & 0 & 0
Allpaths-LG & GAP2SEQ21 & 0 & 16956 & 0 & 72745 & 0 & 0
Bambus2 & GAP2SEQ21 & 0 & 48777 & 0 & 33807 & 0 & 0
MSR-CA & GAP2SEQ21 & 7 & 25347 & 0 & 59152 & 0 & 0
SGA & GAP2SEQ21 & 3 & 348507 & 0 & 4507 & 0 & 0
SOAPdenovo & GAP2SEQ21 & 2 & 36516 & 0 & 69834 & 0 & 0
Velvet & GAP2SEQ21 & 25 & 44248 & 643 & 43039 & 2 & 20
TOTAL & GAP2SEQ21 & 50.0 & 574352.0 & 5169.0 & 419418.0 & 3.0 & 3072.0
AVERAGE & GAP2SEQ21 & 6.2 & 71794.0 & 646.1 & 52427.2 & 0.4 & 384.0
ABySS & SEALER & 7 & 12751 & 7935 & 31079 & 40 & 42373
ABySS2 & SEALER & 6 & 7345 & 0 & 106796 & 30 & 8484
Allpaths-LG & SEALER & 0 & 7791 & 0 & 113102 & 29 & 5466
Bambus2 & SEALER & 0 & 21130 & 0 & 42735 & 76 & 23776
MSR-CA & SEALER & 10 & 17392 & 0 & 70419 & 61 & 8671
SGA & SEALER & 1 & 16345 & 0 & 10903 & 443 & 277170
SOAPdenovo & SEALER & 2 & 36592 & 0 & 69834 & 7 & 2897
Velvet & SEALER & 27 & 19012 & 1270 & 48120 & 85 & 15051
TOTAL & SEALER & 53.0 & 138358.0 & 9205.0 & 492988.0 & 771.0 & 383888.0
AVERAGE & SEALER & 6.6 & 17294.8 & 1150.6 & 61623.5 & 96.4 & 47986.0
ABySS & GAP2SEQ & 8 & 18211 & 4526 & 31180 & 9 & 3097
ABySS2 & GAP2SEQ & 7 & 7508 & 0 & 137725 & 7 & 515
Allpaths-LG & GAP2SEQ & 0 & 6644 & 0 & 149758 & 14 & 529
Bambus2 & GAP2SEQ & 0 & 24395 & 0 & 47162 & 30 & 4655
MSR-CA & GAP2SEQ & 8 & 16492 & 0 & 96378 & 35 & 3049
SGA & GAP2SEQ & 1 & 10745 & 0 & 30015 & 130 & 83850
SOAPdenovo & GAP2SEQ & 2 & 35949 & 0 & 69834 & 4 & 280
Velvet & GAP2SEQ & 27 & 15461 & 643 & 79866 & 40 & 3332
TOTAL & GAP2SEQ & 53.0 & 135405.0 & 5169.0 & 641918.0 & 269.0 & 99307.0
AVERAGE & GAP2SEQ & 6.6 & 16925.6 & 646.1 & 80239.8 & 33.6 & 12413.4
"""

def percent_of(n, o):
    if float(n) == float(o):
        return '+0%'

    if float(n) == 0.:
        return '-100%'

    percent = ((float(n) / float(o)) - 1) * 100
    sign = '+' if float(n) > float(o) else ''
    return sign + str(round(percent, 1)) + '%'

def percent_max(l):
    high = [(float(v[1:-1]), i) for i, v in enumerate(l) if v[0] == '+']
    low = [(float(v[1:-1]), i) for i, v in enumerate(l) if v[0] == '-']

    if len(high) != 0:
        return l[max(high)[1]]
    else:
        return l[min(low)[1]]

def percent_min(l):
    high = [(float(v[1:-1]), i) for i, v in enumerate(l) if v[0] == '+']
    low = [(float(v[1:-1]), i) for i, v in enumerate(l) if v[0] == '-']

    if len(low) != 0:
        return l[max(low)[1]]
    else:
        return l[min(high)[1]]

def color_stripe(tools, index, direction=-1):
    list = [str(tool[index]) for tool in tools]

    max, min = percent_max(list), percent_min(list)

    if min == max:
        return ' & '.join(list)

    high, low = '', ''
    if direction == 1:
        high, low = 'ForestGreen', 'Red'
    elif direction == -1:
        high, low = 'Red', 'ForestGreen'

    for i, v in enumerate(list):
        if v == max:
            list[i] = '{\\color{%s} %s}' % (high, max)
        elif v == min:
            list[i] = '{\\color{%s} %s}' % (low, min)

    return ' & '.join(list)

def latex_print(assembly, original, tools, header=True):
    table = "\\multirow{8}{*}{\\rotatebox[origin=c]{90}{\\bf %s}} && " % (assembly)
    if header:
        table += 'Original & ' + ' & '.join(['{\\bf %s}' % tool[0] for tool in tools]) + ' \\\\ \\bottomrule \n'
    else:
        table += '&' * len(tools) + ' \\\\ \\bottomrule \n'

    #stripe = lambda n, t: ' & '.join([str(tool[n]) for tool in t])

    # +1 = higher is better, -1 = lower is better
    labels = [('Misassemblies', -1), ('Erroneous length', -1), ('Unaligned length', -1),
        ('NGA50', 1), ('Number of gaps', -1), ('Total gap length', -1)]
    for i, label in enumerate(labels):
        table += '& %s & %s & ' % (label[0], original[i]) + color_stripe(tools, i+1, label[1]) + ' \\\\ \n'

    print(table.replace('%', '\%')[:-1])

from collections import defaultdict

original = {}
tools = defaultdict(lambda: defaultdict(list))
for line in text.split('\n'):
    data = line.split(' & ')
    if len(data) < 2: continue
    if data[1] == 'ORIGINAL':
        original[data[0]] = data[2:]
    else:
        tools[data[1]][data[0]] = data[2:]

assemblies = ['ABySS', 'ABySS2', 'Allpaths-LG', 'Bambus2', 'MSR-CA', 'SGA', 'SOAPdenovo', 'Velvet', 'TOTAL', 'AVERAGE']:
for i, assembly in enumerate(assemblies):
    latex_print(assembly, original[assembly], [[tool] + \
        [percent_of(tools[tool][assembly][j], original[assembly][j]) for j in range(len(tools[tool][assembly]))] \
        for tool in tools.keys()], (i == 0))
