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

GAPLENGTHS=()
while read BED; do
  CONTIG=$(echo $BED | cut -f1 -d' ')
  START=$(echo $BED | cut -f2 -d' ')
  END=$(echo $BED | cut -f3 -d' ')

  GAPLENGTH=$(($END - $START - 82))
  GAPLENGTHS+=($GAPLENGTH)

  for ((i=0;i<${#MEANS[@]};++i)); do
    # Extract all known gap-covering reads
    if [ ! -f known."$GAPLENGTH"."${MEANS[i]}".fa ]; then
      $SAMTOOLS view -u $DATA/known_aln"$i".bam "$CONTIG:$START-$END" | \
        $SAMTOOLS fasta - > known."$GAPLENGTH"."${MEANS[i]}".fa
    fi
  done
done < $DATA/gaps.bed

I=0
while read BED; do
  CONTIG=$(echo $BED | cut -f1 -d' ')
  START=$(echo $BED | cut -f2 -d' ')
  END=$(echo $BED | cut -f3 -d' ')

  GAPLENGTH=${GAPLENGTHS[I]}
  I=$(($I + 1))

  for ((i=0;i<${#MEANS[@]};++i)); do
    # Extract all overlapping reads
    if [ ! -f overlap."$GAPLENGTH"."${MEANS[i]}".fa ]; then
      $SAMTOOLS view -u $DATA/aln."${MEANS[i]}".bam "$CONTIG:$START-$END" | \
        $SAMTOOLS fasta - > overlap."$GAPLENGTH"."${MEANS[i]}".fa
    fi

    # Extract filtered reads
    if [ ! -f filter."$GAPLENGTH"."${MEANS[i]}".fa ]; then
      $EXTRACT $DATA/aln."${MEANS[i]}".bam \
        $READLENGTH ${MEANS[i]} ${STDDEVS[i]} \
        $CONTIG $START $END $GAPLENGTH 1 1 25 \
        > filter."$GAPLENGTH"."${MEANS[i]}".fa
    fi

    # Evaluate the schemes
    OVERLAP=$(python $SCRIPTS/evaluate.py \
      aln."${MEANS[i]}".fa \
      known."$GAPLENGTH"."${MEANS[i]}".fa \
      overlap."$GAPLENGTH"."${MEANS[i]}".fa)

    UNMAPPED=$(python $SCRIPTS/evaluate.py \
      aln."${MEANS[i]}".fa \
      known."$GAPLENGTH"."${MEANS[i]}".fa \
      unmapped."${MEANS[i]}".fa)

    FILTER=$(python $SCRIPTS/evaluate.py \
      aln."${MEANS[i]}".fa \
      known."$GAPLENGTH"."${MEANS[i]}".fa \
      filter."$GAPLENGTH"."${MEANS[i]}".fa)

    echo $GAPLENGTH ${MEANS[i]} ${STDDEVS[i]} \
      $OVERLAP $UNMAPPED $FILTER >> results
  done
done < $DATA/breakpoints.bed
