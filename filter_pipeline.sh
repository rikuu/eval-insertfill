#!/bin/bash
set -e

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
SCRIPTS=$DIR/scripts

source $DIR/configuration.sh
OUT=$PREFIX/filter/

source $SCRIPTS/gen_reads.sh
source $SCRIPTS/gen_reference.sh

cd $OUT

for ((i=0;i<${#MEANS[@]};++i)); do
  if [ ! -f aln."${MEANS[i]}".fa ]; then
    # Extract all reads
    $SAMTOOLS fasta $DATA/aln."${MEANS[i]}".bam > aln."${MEANS[i]}".fa

    # Extract all known gap-covering reads
    $SAMTOOLS view -b -L tmp.bed $DATA/known_aln"$i".bam | \
      $SAMTOOLS fasta - > known."${MEANS[i]}".fa

    # Extract all unmapped reads
    $SAMTOOLS view -f4 -b $DATA/aln."${MEANS[i]}".bam | \
      $SAMTOOLS fasta - > unmapped."${MEANS[i]}".fa
  fi
done

while read BED; do
  CONTIG=$(echo $BED | cut -f1 -d' ')
  START=$(echo $BED | cut -f2 -d' ')
  END=$(echo $BED | cut -f3 -d' ')

  echo -e "$CONTIG\t$START\t$END" > tmp.bed

  GAPLENGTH=$(($END - $START - 82))

  for ((i=0;i<${#MEANS[@]};++i)); do
    # Extract all overlapping reads
    $SAMTOOLS view -F4 -b -L tmp.bed $DATA/aln."${MEANS[i]}".bam | \
      $SAMTOOLS fasta - > overlap.fa

    # Extract filtered reads
    LD_PRELOAD=~/htslib/libhts.so.1 $EXTRACT $DATA/aln."${MEANS[i]}".bam \
      $READLENGTH ${MEANS[i]} ${STDDEVS[i]} \
      $CONTIG $START $END > filter.fa

    LD_PRELOAD=~/htslib/libhts.so.1 $EXTRACT1 $DATA/aln."${MEANS[i]}".bam \
      $READLENGTH ${MEANS[i]} ${STDDEVS[i]} \
      $CONTIG $START $END > filter1.fa

    # Evaluate the schemes
    OVERLAP=$(python $SCRIPTS/evaluate.py aln."${MEANS[i]}".fa known."${MEANS[i]}".fa overlap.fa)
    UNMAPPED=$(python $SCRIPTS/evaluate.py aln."${MEANS[i]}".fa known."${MEANS[i]}".fa unmapped.fa)
    FILTER=$(python $SCRIPTS/evaluate.py aln."${MEANS[i]}".fa known."${MEANS[i]}".fa filter.fa)
    FILTER1=$(python $SCRIPTS/evaluate.py aln."${MEANS[i]}".fa known."${MEANS[i]}".fa filter1.fa)
    echo $GAPLENGTH ${MEANS[i]} ${STDDEVS[i]} $OVERLAP $UNMAPPED $FILTER $FILTER1 \
      >> results
  done
done < $DATA/gaps.bed
