#!/usr/bin/env python

import subprocess, multiprocessing
import sys, argparse

# Required tools
samtools = '~/samtools/samtools'
bwa = '~/bwa/bwa'
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

    def data(self):
        return ' ' + self.bam + ' ' + str(self.len) + ' ' + \
            str(self.mu) + ' ' + str(self.sd) + ' '

# Callback function to gather statistics on gaps
filled_gaps = []
def count_filled(result):
    global filled_gaps
    filled_gaps.append(result)
    print

# Runs all the read filtering and gap filling for a single gap
def fill_gap(libraries, seq, id,
        scaffold, start, end, length,
        solid, k, threshold):
    # Cleanup, just to be sure
    subprocess.check_call(['rm', '-f', 'tmp.reads.' + id + '.fasta'])

    # Extract reads
    with open('tmp.extract.' + id + '.log', 'w') as f:
        unmapped = (threshold != -1)
        exact = 1 # Exact is currently faster and more accurate

        gap_data = '%s %i %i %i' % (scaffold, start, end, length)
        extract_params = '%i %i %i' % (exact, unmapped, threshold)
        for lib in libraries:
            subprocess.check_call(extract + lib.data() + gap_data + \
                extract_params + ' >> tmp.reads.' + id + '.fasta', shell=True,
                stderr=f)

    # Run Gap2Seq on the gap with the filtered reads
    with open('tmp.gap2seq.' + id + '.log', 'w') as f:
        subprocess.check_call([gap2seq,
            '-solid', str(solid),
            '-k', str(k),
            '-nb-cores', '1',
            '-filled', 'tmp.fill.' + id + '.fasta',
            '-scaffolds', 'tmp.gap.' + id + '.fasta',
            '-reads', 'tmp.reads.' + id + '.fasta'],
            stderr=subprocess.STDOUT, stdout=f)

    filled = False
    try:
        with open('tmp.gap2seq.' + id + '.log', 'r') as f:
            filled = f.readlines()[-2][7:8] == '1'
    except:
        pass

    subprocess.check_call(['rm', '-f',
        'tmp.reads.' + id + '.fasta',
        'tmp.gap.' + id + '.fasta',
        'tmp.reads.' + id + '.h5'])

    return filled

def parse_gap(bed, gap, id, gap_lengths=None):
    subprocess.check_call(['rm', '-f', 'tmp.gap.' + id + '.fasta'])

    # Write the current gap to a file
    gap_out = open('tmp.gap.' + id + '.fasta', 'w')
    gap_out.write(gap)
    gap_out.close()

    # Parse gap data from bed file
    gap_data = bed.readline().rstrip().split('\t')
    scaffold = gap_data[0]
    start = int(gap_data[1])
    end = int(gap_data[2])

    # Parse gap length
    length = end - start
    if gap_lengths != None:
        gap_data = gap_lengths.readline().rstrip().split('\t')
        length = int(gap_data[2]) - int(gap_data[1])

    return scaffold, start, end, length

# Starts multiple gapfilling processes in parallel
def start_fillers(bed, gaps, libraries, pool=None, async=False,
        gap_lengths=None, k=31, solid=2, threshold=25):
    if async: assert(pool != None)

    gap_id = 0

    seq = ''
    for gap in gaps:
        if gap[0] == '>' and seq != '':
            # Force python to copy the string. This might not affect
            # the result. Best to be sure when multiprocessing
            copy = ''.join(seq)

            scaffold, start, end, length = \
                parse_gap(bed, copy, str(gap_id), gap_lengths)

            if async:
                pool.apply_async(fill_gap, args=([libraries, copy, str(gap_id),
                    scaffold, start, end, length, solid, k, threshold]),
                    callback=count_filled)
            else:
                count_filled(fill_gap(libraries, copy, str(gap_id), scaffold,
                    start, end, length, solid, k, threshold))

            gap_id += 1
            seq = ''
        seq += gap

    scaffold, start, end, length = \
        parse_gap(bed, copy, str(gap_id), gap_lengths)

    if async:
        pool.apply_async(fill_gap, args=([libraries, seq, str(gap_id),
            scaffold, start, end, length, solid, k, threshold]),
            callback=count_filled)
    else:
        count_filled(fill_gap(libraries, seq, str(gap_id), scaffold,
            start, end, length, solid, k, threshold))

    return gap_id+1

def join_filled(gaps, out='filled.fasta'):
    subprocess.check_call('rm -f ' + out + ' && touch ' + out, shell=True)
    for i in range(gaps):
        subprocess.check_call('cat tmp.fill.' + str(i) + \
            '.fasta >> ' + out, shell=True)

        subprocess.check_call(['rm', '-f', 'tmp.fill.' + str(i) + '.fasta'])

def cut_gaps(scaffolds):
    contigs_file = 'tmp.contigs'
    gap_file = 'tmp.gaps'
    bed_file = 'tmp.bed'

    subprocess.check_call(['rm', '-f', contigs_file, gap_file, bed_file])

    subprocess.check_call([gapcutter,
        '-scaffolds', scaffolds,
        '-gaps', gap_file,
        '-contigs', contigs_file,
        '-bed', bed_file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return bed_file, gap_file

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Gap2Seq with read filtering.')
    parser.add_argument('-o', '--out', type=str, default='filled.fasta')
    parser.add_argument('-l', '--libraries', required=True)
    parser.add_argument('-i', '--index', type=int, default=-1)

    parser.add_argument('-k', type=int, default=31)
    parser.add_argument('--solid', type=int, default=2)

    parser.add_argument('-s', '--scaffolds', type=argparse.FileType('r'))

    parser.add_argument('-b', '--bed', type=argparse.FileType('r'))
    parser.add_argument('-g', '--gaps', type=argparse.FileType('r'))

    parser.add_argument('-G', '--gap-lengths', type=argparse.FileType('r'))
    parser.add_argument('-u', '--unmapped-threshold', type=int, default=25)

    parser.add_argument('--async', action='store_true', default=False)
    parser.add_argument('-t', '--threads', type=int, default=1)

    args = vars(parser.parse_args())

    # Short read libraries aligned to the scaffolds
    libraries = []
    with open(args['libraries'], 'r') as f:
        for lib in f:
            arg = lib.split('\t')
            libraries += [Library(arg[0], int(arg[1]), int(arg[2]), int(arg[3]))]

    if args['index'] != -1:
        libraries = [libraries[args['index']]]

    if args['bed'] == None or args['gaps'] == None:
        if args['scaffolds'] != None:
            print('Cutting gaps')
            args['bed'], args['gaps'] = cut_gaps(args['scaffolds'])
        else:
            parser.print_help()
            print('Either [-b/--bed and -g/--gaps] or [-s/--scaffolds] are required.')
            sys.exit(1)

    pool = None
    if args['threads'] < 2:
        args['async'] = False
    else:
        args['async'] = True
        pool = multiprocessing.Pool(args['threads'])

    print('Starting gapfillers')
    gaps = start_fillers(args['bed'], args['gaps'], libraries,
        pool=pool, async=args['async'], gap_lengths=args['gap_lengths'],
        k=args['k'], solid=args['solid'], threshold=args['unmapped_threshold'])

    bed.close()
    gaps.close()
    gap_lengths.close()

    if pool != None:
        pool.close()
        pool.join()

    print('Joining filled sequences')
    join_filled(gaps, args['out'])

    # Merge gaps and contigs back to scaffolds here

    success = 0
    for fg in filled_gaps:
        if fg == True:
            success += 1

    print('Filled ' + str(success) + ' out of ' + str(len(filled_gaps)) + ' gaps')
