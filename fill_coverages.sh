#! /bin/sh

COVERAGES=(1 2 5 10 15 20 25 30 35 40 50)
MEAN=159
STDDEV=18
GENOMELENGTH=10000
GAPLENGTHS=(10 20 50 100 125 150 175 200 225 250 275 500 750 1000 1250 1500 1750 2000 2250 2500 2750 5000)

# Exit if any command fails
set -e

# Generate a genome
rm -f genome.fa
svsim simulate_genome -c contig_1 -e $GENOMELENGTH genome.fa

for i in "${COVERAGES[@]}"; do
  # Generate reads with the different parameters
  rm -f reads_pe1.fa reads_pe2.fa
  svsim simulate_reads -t dwgsim -c $i -m $MEAN -s $STDDEV genome.fa reads

  for j in "${GAPLENGTHS[@]}"; do
    # Generate a gapped genome
    rm -f tmp.* gapped_genome.fa
    python gap.py genome.fa contig_1 $j > gapped_genome.fa

    ../Gap2Seq-2.0/build/GapCutter -bed tmp.bed -gaps tmp.gaps -contigs tmp.contigs -scaffolds gapped_genome.fa
    ../Gap2Seq-2.0/build/Gap2Seq -scaffolds tmp.gaps -reads reads_pe1.fa,reads_pe2.fa -solid 2 -k 31 -filled tmp.filled 2> /dev/null

    FILLED=$(python filled.py tmp.filled)
    echo $i $j $FILLED >> results_fill
  done
done

rm -f genome.fa gapped_genome.fa tmp.* reads_pe*.fa
