text = \
"""ABySS & ORIGINAL & 3 & 190458 & 262068 & 1320 & 1061 & 585628
ABySS2 & ORIGINAL & 99 & 555099 & 157759 & 11869 & 2820 & 949137
Allpaths-LG & ORIGINAL & 95 & 667229 & 36941 & 34534 & 4307 & 3227193
Bambus2 & ORIGINAL & 1584 & 11114542 & 161358 & 3045 & 11809 & 10370362
CABOG & ORIGINAL & 91 & 615239 & 2506 & 46665 & 3043 & 231078
MSR-CA & ORIGINAL & 1110 & 5412965 & 318421 & 5704 & 30622 & 6097928
SGA & ORIGINAL & 8 & 1580489 & 1160159 & 2644 & 21459 & 12840408
SOAPdenovo & ORIGINAL & 1250 & 8449941 & 1306173 & 6592 & 8544 & 10255930
Velvet & ORIGINAL & 9308 & 12531431 & 23484076 & 1793 & 51567 & 63559964
TOTAL & ORIGINAL & 13548.0 & 41117393.0 & 26889461.0 & 114166.0 & 135232.0 & 108117628.0
AVERAGE & ORIGINAL & 1505.3 & 4568599.2 & 2987717.9 & 12685.1 & 15025.8 & 12013069.8
ABySS & GAP2SEQ21 & 6 & 403969 & 129417 & 1310 & 148 & 117718
ABySS2 & GAP2SEQ21 & 122 & 1164541 & 17066 & 10439 & 0 & 0
Allpaths-LG & GAP2SEQ21 & 225 & 2488444 & 59 & 29389 & 0 & 0
Bambus2 & GAP2SEQ21 & 1572 & 3988976 & 55515 & 5363 & 0 & 0
CABOG & GAP2SEQ21 & 87 & 751103 & 2506 & 52655 & 9 & 180
MSR-CA & GAP2SEQ21 & 1099 & 9847021 & 106899 & 3731 & 76 & 17513
SGA & GAP2SEQ21 & 44 & 13968873 & 20388 & 2158 & 0 & 0
SOAPdenovo & GAP2SEQ21 & 1317 & 7812433 & 129790 & 7443 & 11 & 8815
Velvet & GAP2SEQ21 & 11625 & 14285276 & 368622 & 1419 & 2838 & 3543807
TOTAL & GAP2SEQ21 & 16097.0 & 54710636.0 & 830262.0 & 113907.0 & 3082.0 & 3688033.0
AVERAGE & GAP2SEQ21 & 1788.6 & 6078959.6 & 92251.3 & 12656.3 & 342.4 & 409781.4"""

# From the paper
old_tools = [
['GapCloser',
    {'ABySS': ['+133.3%', '+18.2%', '-16.6%', '+1.0%', '-5.9%', '-24.5%'],
    'SGA': ['+375%', '+21.1%', '-83.9%', '+244.2%', '-56.7%', '-53.5%'],
    'ABySS2': ['+18.2%', '+15.9%', '-21.4%', '+4.1%', '-14.5%', '-34.9%'],
    'Allpaths-LG': ['-6.3%', '+34.5%', '-14.3%', '+48.3%', '-35.1%', '-37.9%'],
    'CABOG': ['+16.5%', '+19.0%', '+0%', '+16.2%', '-18.9%', '-50.3%'],
    'Bambus2': ['+3.1%', '-9.6%', '-42.7%', '+34.8%', '-16.4%', '-45.6%'],
    'MSR-CA': ['+14.5%', '+2.8%', '-30.5%', '+73.3%', '-34.9%', '-49.3%'],
    'SOAPdenovo': ['+17.1%', '-1.3%', '-28.8%', '+17.4%', '-25.2%', '-21.3%'],
    'Velvet': ['+26.3%', '-10.4%', '-58.4%', '+104.4%', '-43.4%', '-22.8%']}],
['GapFiller',
    {'ABySS': ['+0%', '+6.1%', '-34.2%', '+0.7%', '-32.5%', '-27.6%'],
    'SGA': ['+112.5%', '-14.6%', '-86.5%', '+238.0%', '-49.9%', '-55.4%'],
    'ABySS2': ['+3.0%', '+3.3%', '-36.4%', '+3.1%', '-38.5%', '-25.3%'],
    'Allpaths-LG': ['+10.5%', '+6.6%', '-11.5%', '+22.5%', '-20.6%', '-17.3%'],
    'CABOG': ['+5.5%', '-0.2%', '+0%', '+64.7%', '-51.9%', '-42.0%'],
    'Bambus2': ['+4.4%', '+0.6%', '-37.8%', '+16.8%', '-2.4%', '-27.1%'],
    'MSR-CA': ['+17.6%', '+9.6%', '-30.6%', '+78.2%', '-47.8%', '-45.4%'],
    'SOAPdenovo': ['+6.1%', '+3.1%', '-27.2%', '+4.1%', '-5.2%', '-15.9%'],
    'Velvet': ['+31.1%', '+41.1%', '-64.5%', '+62.2%', '-26.0%', '-19.2%']}],
['Gap2Seq 1.0',
    {'ABySS': ['+100%', '-9.4%', '-8.4%', '+1.3%', '-33.4%', '-25.5%'],
    'SGA': ['+287.5%', '-24.6%', '-38.6%', '+149.1%', '-51.5%', '-30.2%'],
    'ABySS2': ['+5.1%', '+3.5%', '-15.4%', '+2.4%', '-23.9%', '-36.8%'],
    'Allpaths-LG': ['+14.7%', '-3.0%', '+26.8%', '+23.3%', '-29.8%', '-16.0%'],
    'CABOG': ['+7.7%', '-3.5%', '+0%', '+8.9%', '-13.8%', '-22.6%'],
    'Bambus2': ['+2.1%', '-0.6%', '+2.5%', '+1.8%', '-6.6%', '-4.6%'],
    'MSR-CA': ['+6.8%', '-6.1%', '-13.1%', '+32.9%', '-27.8%', '-18.1%'],
    'SOAPdenovo': ['+11.7%', '-0.5%', '-14.1%', '+4.0%', '-18.8%', '-6.3%'],
    'Velvet': ['+13.3%', '-0.5%', '-20.7%', '+27.0%', '-26.6%', '-4.2%']}]
]

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
        table += '{\\bf Original} & ' + ' & '.join(['{\\bf %s}' % tool[0] for tool in tools]) + ' \\\\ \\bottomrule \n'
    else:
        table += '&' * len(tools) + ' \\\\ \\midrule \n'

    #stripe = lambda n, t: ' & '.join([str(tool[n]) for tool in t])

    # +1 = higher is better, -1 = lower is better
    labels = [('Misassemblies', -1), ('Erroneous length', -1), ('Unaligned length', -1),
        ('NGA50', 1), ('Number of gaps', -1), ('Total gap length', -1)]
    for i, label in enumerate(labels):
        table += '& %s & %s & ' % (label[0], original[i]) + color_stripe(tools, i+1, label[1]) + ' \\\\ \n'

    print(table.replace('%', '\%')[:-1])

original = {}
new = {}
for line in text.split('\n'):
    data = line.split(' & ')
    if data[1] == 'ORIGINAL':
        original[data[0]] = data[2:]
    else:
        new[data[0]] = data[2:]

# Split into 3 figures (lex-order)
split_assemblies = [['ABySS', 'ABySS2', 'Allpaths-LG'],
    ['Bambus2', 'CABOG', 'MSR-CA'],
    ['SGA', 'SOAPdenovo', 'Velvet']] #, 'TOTAL', 'AVERAGE']:

for assemblies in split_assemblies:
    print('\\begin{table}[h]\n\t\\centering\n\t\\begin{tabu}{clrrrrr}\n\t\t\\toprule')
    for j, assembly in enumerate(assemblies):
        gap2seq2 = ['Gap2Seq 2.1']+[percent_of(new[assembly][i], original[assembly][i]) for i in range(len(new[assembly]))]
        tools = [[tool[0]] + tool[1][assembly] for tool in old_tools] + [gap2seq2]
        latex_print(assembly, original[assembly], tools, (j == 0))
    print('\t\t\\bottomrule\n\t\t\\end{tabu}\n\t\\caption{Quality of filled gaps on Human14 dataset.}\n\t\\label{fig:human14}\n\\end{table}')
