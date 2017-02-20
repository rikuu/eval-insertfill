import sys, pysam

samfile = pysam.AlignmentFile(sys.argv[1], "rb")

contig = sys.argv[2]
gap_start = int(sys.argv[3])
gap_end = int(sys.argv[4])

max_distance = int(sys.argv[5])

for read in samfile.fetch(contig, gap_start-max_distance, gap_end+max_distance):
    if read.is_paired and read.mate_is_unmapped:
        pair = 2 if read.is_read1 else 1
        print('>%s/%i' % (read.query_name, pair))
