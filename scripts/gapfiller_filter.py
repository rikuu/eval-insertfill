import sys, pysam

samfile = pysam.AlignmentFile(sys.argv[1], "rb")

gap_start = int(sys.argv[2])
gap_end = int(sys.argv[3])

max_distance = int(sys.argv[4])

for read in samfile.fetch():
    if read.is_paired and read.mate_is_unmapped:
        distance = min(abs(gap_start - read.query_alignment_start), abs(gap_end - read.query_alignment_start))
        if distance < max_distance:
            print('>'+read.query_name)
