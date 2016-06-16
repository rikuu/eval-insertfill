#!/usr/local/bin/bash

# Required tools
DWGSIM=dwgsim
BEDTOOLS=bedtools
SAMTOOLS=samtools
BWA=bwa
GAP2SEQ=../Gap2Seq-2.0/build/Gap2Seq
GAPCUTTER=../Gap2Seq-2.0/build/GapCutter
EXTRACT=../extract/extract

# Parameters for the testing
READLENGTH=100
COVERAGE=30
CONTIG=chr17
MEANS=(150 1500 3000)
STDDEVS=(15 150 300)

# Exit if any command fails
set -e

# Index the reference genome
if [ ! -f chr17.fa.amb ]; then
  $BWA index chr17.fa
fi

for ((i=0;i<${#MEANS[@]};++i)); do
  # Generate reads with the different parameters
  if [ ! -f reads"$i"_pe1.fq ]; then
    $DWGSIM -i -1 $READLENGTH -2 $READLENGTH -d ${MEANS[i]} -s ${STDDEVS[i]} \
      -C $COVERAGE chr17.fa reads"$i"
    rm reads"$i".bfast*

    mv reads"$i".bwa.read1.fastq reads"$i"_pe1.fq
    mv reads"$i".bwa.read2.fastq reads"$i"_pe2.fq
  fi

  # Map-sort-index the reads
  if [ ! -f known_aln"$i".bam ]; then
    $BWA mem -t 16 -I ${MEANS[i]},${STDDEVS[i]} chr17.fa reads"$i"_pe1.fq \
        reads"$i"_pe2.fq | \
      $SAMTOOLS view -Shu - | \
      $SAMTOOLS sort - | \
      $SAMTOOLS rmdup -s - - > known_aln"$i".bam
    $SAMTOOLS index known_aln"$i".bam
  fi
done

while read BED; do
  CONTIG=$(echo $BED | cut -f1 -d' ')
  START=$(echo $BED | cut -f2 -d' ')
  END=$(echo $BED | cut -f3 -d' ')

  echo -e "$CONTIG\t$START\t$END" > tmp.bed

  GAPLENGTH=$(($END - $START - 82))

  for ((i=0;i<${#MEANS[@]};++i)); do
    # Extract all reads
    $SAMTOOLS fasta aln."${MEANS[i]}".bam > aln.fa

    # Extract all overlapping reads
    $SAMTOOLS view -F4 -b -L tmp.bed aln."${MEANS[i]}".bam | \
      $SAMTOOLS fasta - > overlap.fa

    # Extract all unmapped reads
    $SAMTOOLS view -f4 -b aln."${MEANS[i]}".bam | \
      $SAMTOOLS fasta - > unmapped.fa

    # Extract all known gap-covering reads
    $SAMTOOLS view -b -L tmp.bed known_aln"$i".bam | \
      $SAMTOOLS fasta - > known.fa

    # Extract filtered reads
    LD_PRELOAD=~/htslib/libhts.so.1 $EXTRACT aln."${MEANS[i]}".bam \
      $READLENGTH ${MEANS[i]} ${STDDEVS[i]} \
      $CONTIG $START $END > filter.fa

    LD_PRELOAD=~/htslib/libhts.so.1 $EXTRACT1 aln."${MEANS[i]}".bam \
      $READLENGTH ${MEANS[i]} ${STDDEVS[i]} \
      $CONTIG $START $END > filter1.fa

    # Evaluate the schemes
    OVERLAP=$(python evaluate.py aln.fa known.fa overlap.fa)
    UNMAPPED=$(python evaluate.py aln.fa known.fa unmapped.fa)
    FILTER=$(python evaluate.py aln.fa known.fa filter.fa)
    FILTER1=$(python evaluate.py aln.fa known.fa filter1.fa)
    echo $GAPLENGTH ${MEANS[i]} ${STDDEVS[i]} $OVERLAP $UNMAPPED $FILTER $FILTER1 \
      >> results
  done
done < tmp.gaps.bed
