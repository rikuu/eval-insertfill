#!/usr/bin/env python
import subprocess, multiprocessing

# Required tools
gapmerger = '/cs/work/scratch/riqwalve/Gap2Seq/build/GapMerger'
gapcutter = '/cs/work/scratch/riqwalve/Gap2Seq/build/GapCutter'
gap2seq = '/cs/work/scratch/riqwalve/Gap2Seq/build/Gap2Seq'
extract = '/cs/work/scratch/riqwalve/extract/extract'

# An object for holding all the data for a library of short reads
class Library:
    def __init__(self, bam, read_length, mean_insert_size, std_dev):
        self.bam = bam
        self.len = read_length
        self.mu = mean_insert_size
        self.sd = std_dev

        # TODO: Assert bam.bai exists

    def data(self):
        return '%s %i %i %i' % (self.bam, self.len, self.mu, self.sd)

class Gap:
    def __init__(self, bed, gap, comment, id):
        (self.scaffold, self.start, self.end), \
            (self.left, self.right, self.length), \
            self.comment, self.id = bed, gap, comment, id

    def data(self):
        return '%s %i %i %i' % (self.scaffold, self.start, self.end, self.length)

    def filler_data(self):
        return ['-left', self.left,
            '-right', self.right,
            '-length', str(self.length)]

# Callback function to gather statistics on gaps
filled_gaps = []
def count_filled(result):
    global filled_gaps
    filled_gaps.append(result)

# Runs all the read filtering and gap filling for a single gap
def fill_gap(libraries, gap, k, fuz, solid, threshold, max_mem):
    # Cleanup, just to be sure
    subprocess.check_call(['rm', '-f', 'tmp.reads.' + gap.id + '.fasta'])

    # Extract reads
    with open('tmp.extract.' + gap.id + '.log', 'w') as f:
        unmapped = (threshold != -1)
        exact = 1 # Exact is currently faster and more accurate

        extract_params = '%i %i %i' % (exact, unmapped, threshold)
        for lib in libraries:
            subprocess.check_call('%s %s %s %s >> tmp.reads.%s.fasta' %
                (extract, lib.data(), gap.data(), extract_params, gap.id),
                shell=True, stderr=f)

    # Run Gap2Seq on the gap with the filtered reads
    fill = ''
    with open('tmp.gap2seq.' + gap.id + '.log', 'w') as f:
        fill = subprocess.check_output([gap2seq,
            '-k', str(k),
            '-fuz', str(fuz),
            '-solid', str(solid),
            '-nb-cores', '1',
            '-max-mem', str(max_mem),
            '-reads', 'tmp.reads.' + gap.id + '.fasta'] + gap.filler_data(),
            stderr=f)

    # Gap2Seq outputs 'Gap2Seq\n' to stdout when exiting
    filled = False
    fill = fill.split(b'\n')
    if len(fill) > 1:
        filled = True
        fill = ''.join([seq.decode() for seq in fill[:-2]])
    else:
        fill = gap.left + ('N' * gap.length) + gap.right

    # Cleanup reads and graph
    subprocess.check_call(['rm', '-f',
        'tmp.reads.' + gap.id + '.fasta',
        'tmp.reads.' + gap.id + '.h5'])

    return (filled, gap.comment, fill)

# NOTE: Assumes gaps and the bed file are both sorted by position
def parse_gap(bed, gap, id, flank_length):
    gap = gap.split('\n')
    comment = gap[0]

    gap = ''.join(gap[1:])

    # TODO: Not assume static flank lengths
    left = gap[:flank_length]
    right = gap[-flank_length:]
    length = len(gap) - 2*flank_length

    # Parse gap data from bed file
    gap_data = bed.readline().rstrip().split('\t')
    scaffold = gap_data[0]
    start = int(gap_data[1])
    end = int(gap_data[2])

    return Gap((scaffold, start, end), (left, right, length), comment, id)

# Starts multiple gapfilling processes in parallel
def start_fillers(bed, gaps, libraries, pool=None, async=False,
        k=31, fuz=10, solid=2, threshold=25, max_mem=20):
    if async: assert(pool != None)

    flank_length = k+fuz

    gap_id = 0

    seq = ''
    for gap in gaps:
        if gap[0] == '>' and seq != '':
            gap_object = parse_gap(bed, seq, str(gap_id), flank_length)

            if async:
                pool.apply_async(fill_gap, args=([libraries, gap_object,
                    k, fuz, solid, threshold, max_mem]),
                    callback=count_filled)
            else:
                count_filled(fill_gap(libraries, gap_object,
                    k, fuz, solid, threshold, max_mem))

            gap_id += 1
            seq = ''
        seq += gap

    gap_object = parse_gap(bed, seq, str(gap_id), flank_length)

    if async:
        pool.apply_async(fill_gap, args=([libraries, gap_object,
            k, fuz, solid, threshold, max_mem]),
            callback=count_filled)
    else:
        count_filled(fill_gap(libraries, gap_object,
            k, fuz, solid, threshold, max_mem))

    return gap_id+1

def join_filled(out='filled.fasta'):
    num_success = 0
    with open(out, 'w') as f:
        for (filled, comment, fill) in filled_gaps:
            if filled: num_success += 1
            f.write(comment + '\n' + fill + '\n')
    return num_success

def cut_gaps(scaffolds):
    contigs_file = 'tmp.contigs'
    gap_file = 'tmp.gaps'
    bed_file = 'tmp.bed'

    subprocess.check_call(['rm', '-f', contigs_file, gap_file, bed_file])

    subprocess.check_call([gapcutter,
            '-scaffolds', scaffolds,
            '-gaps', gap_file,
            '-contigs', contigs_file,
            '-bed', bed_file],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return open(bed_file, 'r'), open(gap_file, 'r')

def merge_gaps(filled, merged):
    contigs_file = 'tmp.contigs'

    subprocess.check_call([gapmerger,
            '-scaffolds', merged,
            '-gaps', filled,
            '-contigs', contigs_file],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # subprocess.check_call(['rm', '-f', contigs_file, filled])

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run Gap2Seq with read filtering.')
    parser.add_argument('-o', '--out', type=str, default='filled.fasta')
    parser.add_argument('-l', '--libraries', required=True)
    parser.add_argument('-i', '--index', type=int, default=-1)
    parser.add_argument('-t', '--threads', type=int, default=1)
    parser.add_argument('-u', '--unmapped-threshold', type=int, default=25)

    parser.add_argument('-k', type=int, default=31)
    parser.add_argument('--fuz', type=int, default=10)
    parser.add_argument('--solid', type=int, default=2)
    parser.add_argument('--max-mem', type=int, default=20)

    parser.add_argument('-s', '--scaffolds', type=str)

    parser.add_argument('-b', '--bed', type=argparse.FileType('r'))
    parser.add_argument('-g', '--gaps', type=argparse.FileType('r'))

    args = vars(parser.parse_args())

    # Short read libraries aligned to the scaffolds
    libraries = []
    with open(args['libraries'], 'r') as f:
        for lib in f:
            arg = lib.split('\t')
            libraries += [Library(arg[0], int(arg[1]), int(arg[2]), int(arg[3]))]

    # Use only 1 library
    if args['index'] != -1:
        libraries = [libraries[args['index']]]

    scaffolds_cut = False
    if args['bed'] == None or args['gaps'] == None:
        if args['scaffolds'] != None:
            print('Cutting gaps')
            args['bed'], args['gaps'] = cut_gaps(args['scaffolds'])
            scaffolds_cut = True
        else:
            parser.print_help()
            print('Either [-b/--bed and -g/--gaps] or [-s/--scaffolds] are required.')
            exit(1)

    pool = multiprocessing.Pool(args['threads'])

    # Gap2Seq divides the max mem evenly between threads, but as we run multiple
    # parallel instances with 1 thread, we need to pre-divide
    max_mem = int(args['max_mem'] / args['threads'])

    print('Starting gapfillers')
    gaps = start_fillers(args['bed'], args['gaps'], libraries,
        pool=pool, async=(args['threads'] > 1),
        k=args['k'], fuz=args['fuz'], solid=args['solid'],
        threshold=args['unmapped_threshold'],
        max_mem=max_mem)

    args['bed'].close()
    args['gaps'].close()

    pool.close()
    pool.join()

    num_success = 0
    if scaffolds_cut:
        print('Merging filled gaps and contigs')
        num_success = join_filled('tmp.filled')
        merge_gaps('tmp.filled', args['out'])
    else:
        print('Joining filled sequences')
        num_success = join_filled(args['out'])

    print('Filled %i out of %i gaps' % (num_success, gaps))
