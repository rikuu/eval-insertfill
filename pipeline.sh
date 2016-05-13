#!/usr/local/bin/bash

COVERAGE=30
MEANS=(150 1500 3000)
STDDEVS=(15 150 300)
GENOMELENGTH=10000
# GAPLENGTHS=(10 20 30 40 50 100 110 120 130 140 150 175 200 210 220 230 240 250 260 270 280 290 300 \
#  400 450 500 550 600 650 700 750 800 850 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 \
#  2250 2500 2750 3000 3500 4000 4500 5000)

# Exit if any command fails
set -e

# Generate a genome
rm -f genome.fa reads_pe1.fa reads_pe2.fa known_aln*
svsim simulate_genome -c contig_1 -e $GENOMELENGTH genome.fa

for ((i=0;i<${#MEANS[@]};++i)); do
  # Generate reads with the different parameters
  svsim simulate_reads -t dwgsim -c $COVERAGE -m ${MEANS[i]} -s ${STDDEVS[i]} genome.fa reads"$i"

  # Map-sort-index the reads
  svsim map_reads -pe1 reads"$i"_pe1.fa -pe2 reads"$i"_pe2.fa genome.fa known_aln"$i"
  samtools index known_aln"$i".bam
done

#for GAPLENGTH in "${GAPLENGTHS[@]}"; do
for ((GAPLENGTH=10; GAPLENGTH<=5000; GAPLENGTH+=10)); do
  # Generate a gapped genome
  rm -f aln.* tmp.* known.* overlap.* unmapped.* filter.* gapped_genome.*
  python gap.py genome.fa contig_1 $GAPLENGTH > gapped_genome.fa
  bwa index gapped_genome.fa

  # Create a bed file (simulates a real pipeline, I guess)
  ../Gap2Seq-2.0/build/GapCutter -bed tmp.bed -gaps tmp.gaps -contigs tmp.contigs -scaffolds gapped_genome.fa
  CONTIG=$(cat tmp.bed | cut -f1)
  START=$(cat tmp.bed | cut -f2)
  END=$(cat tmp.bed | cut -f3)

  rm -f aln.* known.* overlap.* unmapped.* filter.*
  for ((i=0;i<${#MEANS[@]};++i)); do
    # Map-sort-index reads to the gapped genome/assembly
    rm -f aln.fa overlap.fa unmapped.fa
    for PE in pe1 pe2; do
      bwa mem -t 4 -I ${MEANS[i]},${STDDEVS[i]} gapped_genome.fa reads"$i"_"$PE".fa | \
        samtools view -Shu - | \
        samtools sort - | \
        samtools rmdup -s - - > "$PE".bam
      samtools index "$PE".bam

      # Extract all reads
      samtools fasta "$PE".bam >> aln.fa

      # Extract all overlapping reads
      samtools view -F4 -b -L tmp.bed "$PE".bam | \
        samtools fasta - >> overlap.fa

      # Extract all unmapped reads
      samtools view -f4 -b "$PE".bam | \
        samtools fasta - >> unmapped.fa
    done

    # Extract all known gap-covering reads
    samtools view -b -L tmp.bed known_aln"$i".bam | \
      samtools fasta - > known.fa

    # Extract filtered reads
    ../extract/extract pe1.bam pe2.bam 101 ${MEANS[i]} ${STDDEVS[i]} $CONTIG $START $END > filter.fa

    # Evaluate the schemes
    OVERLAP=$(python evaluate.py aln.fa known.fa overlap.fa)
    UNMAPPED=$(python evaluate.py aln.fa known.fa unmapped.fa)
    FILTER=$(python evaluate.py aln.fa known.fa filter.fa)
    echo $GAPLENGTH ${MEANS[i]} ${STDDEVS[i]} $OVERLAP $UNMAPPED $FILTER >> results
  done
done

rm -f genome.fa reads*_pe1.fa reads*_pe2.fa known_aln*.bam*
rm -f aln.* tmp.* known_aln.* known.* overlap.* unmapped.* filter.* gapped_genome.fa
