#!/bin/bash
set -e

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
SCRIPTS=$DIR/scripts

source $DIR/configuration.sh
OUT=$PREFIX/filter/

source $SCRIPTS/gen_reads.sh
source $SCRIPTS/gen_reference.sh

mkdir -p $OUT || true
cd $OUT

for ((i=0;i<${#MEANS[@]};++i)); do
  # Extract all reads
  if [ ! -f aln."${MEANS[i]}".fa ]; then
    $SAMTOOLS fasta $DATA/aln."${MEANS[i]}".bam > aln."${MEANS[i]}".fa
  fi

  # Extract all unmapped reads
  if [ ! -f unmapped."${MEANS[i]}".fa ]; then
    $SAMTOOLS view -f4 -u $DATA/aln."${MEANS[i]}".bam | \
      $SAMTOOLS fasta - > unmapped."${MEANS[i]}".fa
  fi
done

while read BED; do
  CONTIG=$(echo $BED | cut -f1 -d' ')
  START=$(echo $BED | cut -f2 -d' ')
  END=$(echo $BED | cut -f3 -d' ')

  FLANKLENGTH=41
  GAPLENGTH=$(($END - $START - 82))
  BREAKPOINT=$(($START + $FLANKLENGTH))

  for ((i=0;i<${#MEANS[@]};++i)); do
    # Extract all known gap-covering reads
    if [ ! -f known."$GAPLENGTH"."${MEANS[i]}".fa ]; then
      $SAMTOOLS view -u $DATA/known_aln"$i".bam "$CONTIG:$START-$END" | \
        $SAMTOOLS fasta - > known."$GAPLENGTH"."${MEANS[i]}".fa
    fi

    # Extract all overlapping reads
    if [ ! -f overlap."$GAPLENGTH"."${MEANS[i]}".fa ]; then
      $SAMTOOLS view -u $DATA/aln."${MEANS[i]}".bam "$CONTIG:$START-$END" | \
        $SAMTOOLS fasta - > overlap."$GAPLENGTH"."${MEANS[i]}".fa
    fi

    # Filter reads
    if [ ! -f filter."$GAPLENGTH"."${MEANS[i]}".fa ]; then
      $EXTRACT -bam $DATA/aln."${MEANS[i]}".bam \
        -read-length $READLENGTH -mean ${MEANS[i]} -std-dev ${STDDEVS[i]} \
        -scaffold $CONTIG -breakpoint $BREAKPOINT -flank-length $FLANKLENGTH \
        -gap-length $GAPLENGTH -unmapped $THRESHOLD \
        -reads filter."$GAPLENGTH"."${MEANS[i]}".fa
    fi

    # Evaluate the schemes
    RESULTS=$(python3 $SCRIPTS/evaluate.py \
      aln."${MEANS[i]}".fa \
      known."$GAPLENGTH"."${MEANS[i]}".fa \
      overlap."$GAPLENGTH"."${MEANS[i]}".fa \
      unmapped."${MEANS[i]}".fa \
      filter."$GAPLENGTH"."${MEANS[i]}".fa)

    echo $GAPLENGTH ${MEANS[i]} ${STDDEVS[i]} $RESULTS >> results
  done
done < $DATA/gaps.bed
