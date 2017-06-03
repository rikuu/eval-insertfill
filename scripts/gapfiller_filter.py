import sys, pysam

samfile = pysam.AlignmentFile(sys.argv[1], "rb")

region = sys.argv[2]
colon, dash = region.find(':'), region.find('-')

contig = region[:colon]
start = int(region[colon+1:dash])
end = int(region[dash+1:])

max_distance = int(sys.argv[3])

for read in samfile.fetch(contig, start-max_distance, end+max_distance):
    if read.is_paired and read.mate_is_unmapped:
        pair = 2 if read.is_read1 else 1
        print('>%s/%i' % (read.query_name, pair))
