#!/bin/bash
set -e

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
SCRIPTS=$DIR/scripts

source $DIR/configuration.sh

source $SCRIPTS/gen_reads.sh
source $SCRIPTS/gen_reference.sh

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
    OVERLAP=$(python $SCRIPTS/evaluate.py aln.fa known.fa overlap.fa)
    UNMAPPED=$(python $SCRIPTS/evaluate.py aln.fa known.fa unmapped.fa)
    FILTER=$(python $SCRIPTS/evaluate.py aln.fa known.fa filter.fa)
    FILTER1=$(python $SCRIPTS/evaluate.py aln.fa known.fa filter1.fa)
    echo $GAPLENGTH ${MEANS[i]} ${STDDEVS[i]} $OVERLAP $UNMAPPED $FILTER $FILTER1 \
      >> results
  done
done < gaps.bed
