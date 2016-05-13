#!/usr/local/bin/bash

# Required tools
DWGSIM=dwgsim
BEDTOOLS=bedtools
SAMTOOLS=samtools
BWA=bwa
GAP2SEQ=../Gap2Seq-2.0/build/Gap2Seq
GAP2CUTTER=../Gap2Seq-2.0/build/GapCutter
EXTRACT=../extract/extract

# Parameters for the testing
READLENGTH=100
COVERAGE=30
MEANS=(150 1500 3000)
STDDEVS=(15 150 300)
GENOMELENGTH=10000
# GAPLENGTHS=(10 20 30 40 50 100 110 120 130 140 150 175 200 210 220 230 240 250 260 270 280 290 300 \
#  400 450 500 550 600 650 700 750 800 850 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 \
#  2250 2500 2750 3000 3500 4000 4500 5000)

# Exit if any command fails
set -e

for ((i=0;i<${#MEANS[@]};++i)); do
  # Generate reads with the different parameters
  if [ ! -f reads"$i"_pe1.fa ]; then
    $DWGSIM -i -1 $READLENGTH -2 $READLENGTH -d ${MEANS[i]} -s ${STDDEVS[i]} \
      -C $COVERAGE chr17.fa reads"$i"
    rm reads"$i".bfast*

    mv reads"$i".bwa.read1.fastq reads"$i"_pe1.fq
    mv reads"$i".bwa.read2.fastq reads"$i"_pe2.fq
  fi

  # Map-sort-index the reads
  if [ ! -f known_aln"$i".bam ]; then
    $BWA mem -t 4 -I ${MEANS[i]},${STDDEVS[i]} chr17.fa reads"$i"_pe1.fq \
        reads"$i"_pe2.fq | \
      $SAMTOOLS view -Shu - | \
      $SAMTOOLS sort - | \
      $SAMTOOLS rmdup -s - - > aln.bam
    $SAMTOOLS index known_aln"$i".bam
  fi
done

#for GAPLENGTH in "${GAPLENGTHS[@]}"; do
for ((GAPLENGTH=10; GAPLENGTH<=5000; GAPLENGTH+=10)); do
  rm -f aln.* tmp.* known.* overlap.* unmapped.* filter.* gapped_genome.*

  # Generate a gapped genome
  START=$(((RANDOM % 8000000) + 100000))
  END=$((START+GAPLENGTH))
  echo -e "$CONTIG\t$START\t$END" > tmp.bed
  $BEDTOOLS maskfasta -fi chr17.fa -bed tmp.bed -fo gapped_genome.fa
  $BWA index gapped_genome.fa

  # Create a bed file (simulates a real pipeline, I guess)
  $GAPCUTTER -bed tmp.bed -gaps tmp.gaps -contigs tmp.contigs \
    -scaffolds gapped_genome.fa
  CONTIG=$(cat tmp.bed | cut -f1)
  START=$(cat tmp.bed | cut -f2)
  END=$(cat tmp.bed | cut -f3)

  for ((i=0;i<${#MEANS[@]};++i)); do
    rm -f aln.bam aln.fa known.fa overlap.fa unmapped.fa filter.fa

    # Map-sort-index reads to the gapped genome/assembly
    $BWA mem -t 4 -I ${MEANS[i]},${STDDEVS[i]} gapped_genome.fa \
        reads"$i"_pe1.fq reads"$i"_pe2.fq | \
      $SAMTOOLS view -Shu - | \
      $SAMTOOLS sort - | \
      $SAMTOOLS rmdup -s - - > aln.bam
    $SAMTOOLS index aln.bam

    # Extract all reads
    $SAMTOOLS fasta aln.bam > aln.fa

    # Extract all overlapping reads
    $SAMTOOLS view -F4 -b -L tmp.bed aln.bam | \
      $SAMTOOLS fasta - > overlap.fa

    # Extract all unmapped reads
    $SAMTOOLS view -f4 -b aln.bam | \
      $SAMTOOLS fasta - > unmapped.fa

    # Extract all known gap-covering reads
    $SAMTOOLS view -b -L tmp.bed known_aln"$i".bam | \
      $SAMTOOLS fasta - > known.fa

    # Extract filtered reads
    $EXTRACT aln.bam 101 ${MEANS[i]} ${STDDEVS[i]} $CONTIG $START $END \
      > filter.fa

    # Evaluate the schemes
    OVERLAP=$(python evaluate.py aln.fa known.fa overlap.fa)
    UNMAPPED=$(python evaluate.py aln.fa known.fa unmapped.fa)
    FILTER=$(python evaluate.py aln.fa known.fa filter.fa)
    echo $GAPLENGTH ${MEANS[i]} ${STDDEVS[i]} $OVERLAP $UNMAPPED $FILTER \
      >> results
  done
done

rm -f aln.* tmp.* known_aln.* known.* overlap.* unmapped.* filter.* gapped_genome.fa
