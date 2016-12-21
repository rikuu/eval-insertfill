def percent_of(n, o):
    if float(n) == float(o):
        return '+0%'

    if float(n) == 0.:
        return '-100%'

    if float(o) == 0.:
        o = 1.

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

latex_table = 'quality_table_staph.tex'

original = {}
tools = defaultdict(lambda: defaultdict(list))
with open(latex_table, 'r') as f:
    for line in f:
        data = line.split(' & ')
        data[-1] = data[-1][:data[-1].find(' ')]
        
        if len(data) < 2: continue
        if data[1] == 'ORIGINAL':
            original[data[0]] = data[2:]
        else:
            tools[data[1]][data[0]] = data[2:]

# assemblies = ['ABySS', 'ABySS2', 'Allpaths-LG', 'Bambus2', 'MSR-CA', 'SGA', 'SOAPdenovo', 'Velvet', 'TOTAL', 'AVERAGE']:
# for i, assembly in enumerate(assemblies):
#     latex_print(assembly, original[assembly], [[tool] + \
#         [percent_of(tools[tool][assembly][j], original[assembly][j]) for j in range(len(tools[tool][assembly]))] \
#         for tool in tools.keys()], (i == 0))

# Split into 3 figures (lex-order)
split_assemblies = [['ABySS', 'ABySS2', 'Allpaths-LG'], ['Bambus2', 'MSR-CA', 'SGA'], ['SOAPdenovo', 'Velvet']]
labels = [('GAPCLOSER', 'GapCloser'), ('GAPFILLER_BWA', 'GapFiller'),
    ('GAP2SEQ', 'Gap2Seq'), ('GAP2SEQ21', 'Gap2Seq (with filter)')]
for assemblies in split_assemblies:
    print('\\begin{table}[h]\n\t\\centering\n\t\\begin{tabu}{clrrrrr}\n\t\t\\toprule')
    for j, assembly in enumerate(assemblies):
        res = [[tool[1]] + [percent_of(tools[tool[0]][assembly][j], original[assembly][j]) for j in range(len(tools[tool[0]][assembly]))] for tool in labels]
        latex_print(assembly, original[assembly], res, (j == 0))
    print('\t\t\\bottomrule\n\t\t\\end{tabu}\n\t\\caption{Quality of filled gaps on Human14 dataset.}\n\t\\label{fig:human14}\n\\end{table}')
