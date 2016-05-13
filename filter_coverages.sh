#! /bin/sh

# Parameters inspired by Human14 reads from GAGE
COVERAGE=30
MEANS=(100 150 200 250 300 400 500 600 700 800 900 1000 1500 2000 2500 3000 3500 4000 4500 5000)
STDDEVS=(10 15 20 25 30 40 50 60 70 80 90 100 150 200 250 300 350 400 450 500)
GENOMELENGTH=10000
GAPLENGTHS=(10 20 50 100 125 150 175 200 225 250 275 500 750 1000 1250 1500 1750 2000 2250 2500 2750 5000)

# Exit if any command fails
set -e

# Generate a genome
rm -f genome.fa
svsim simulate_genome -c contig_1 -e $GENOMELENGTH genome.fa

for j in "${GAPLENGTHS[@]}"; do
  # Generate a gapped genome
  rm -f gapped_genome.fa* tmp.*
  python gap.py genome.fa contig_1 $j > gapped_genome.fa
  bwa index gapped_genome.fa

  # Create a bed file
  ../Gap2Seq-2.0/build/GapCutter -bed tmp.bed -gaps tmp.gaps -contigs tmp.contigs -scaffolds gapped_genome.fa
  CONTIG=$(cat tmp.bed | cut -f1)
  START=$(cat tmp.bed | cut -f2)
  END=$(cat tmp.bed | cut -f3)

  for ((i=0;i<${#MEANS[@]};++i)); do
    # Generate reads with the different parameters
    rm -f reads_pe*.fa aln.* overlap.fa unmapped.fa filter.fa
    svsim simulate_reads -t dwgsim -c $COVERAGE -m ${MEANS[i]} -s ${STDDEVS[i]} genome.fa reads

    # Map-sort-index reads to the gapped genome/assembly
    # svsim map_reads -pe1 reads_pe1.fa -pe2 reads_pe2.fa gapped_genome.fa aln
    bwa mem -t 4 -I ${MEANS[i]},${STDDEVS[i]} gapped_genome.fa reads_pe1.fa reads_pe2.fa | \
    samtools view -Shu - | \
    samtools sort - | \
    samtools rmdup -s - - > aln.bam

    samtools index aln.bam

    # Extract all overlapping reads
    samtools view -F4 -b -L tmp.bed aln.bam | \
    samtools fasta - > overlap.fa

    # Extract all unmapped reads
    samtools view -f4 -b aln.bam | \
    samtools fasta - > unmapped.fa

    # Extract filtered reads
    ../bam_bed_extract aln.bam 101 ${MEANS[i]} ${STDDEVS[i]} $CONTIG $START $END > filter.fa

    # Evaluate the schemes
    OVERLAP=$(python coverage.py $START $END overlap.fa)
    UNMAPPED=$(python coverage.py $START $END unmapped.fa)
    FILTER=$(python coverage.py $START $END filter.fa)
    echo ${MEANS[i]} $j $OVERLAP $UNMAPPED $FILTER >> results_coverage
  done
done

rm -f genome.fa reads_pe*.fa aln.* overlap.fa unmapped.fa filter.fa
rm -f tmp.* aln.bam* gapped_genome.fa*
