#!/bin/bash
shopt -s extglob
set -eu
set -o pipefail

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
    $SAMTOOLS fasta $DATA/aln."${MEANS[i]}".bam | \
      grep '^>' > aln."${MEANS[i]}".fa
  fi

  # Extract all unmapped reads
  if [ ! -f unmapped."${MEANS[i]}".fa ]; then
    $SAMTOOLS view -f4 -u $DATA/aln."${MEANS[i]}".bam | $SAMTOOLS fasta - | \
      grep '^>' > unmapped."${MEANS[i]}".fa
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
        $SAMTOOLS fasta - | grep '^>' > known."$GAPLENGTH"."${MEANS[i]}".fa
    fi
  done
done < $DATA/gaps.bed

I=0
while read BED; do
  # Using breakpoints means END = START+82
  CONTIG=$(echo $BED | cut -f1 -d' ')
  START=$(echo $BED | cut -f2 -d' ')
  # END=$(echo $BED | cut -f3 -d' ')

  GAPLENGTH=${GAPLENGTHS[I]}
  I=$(($I + 1))

  echo -e "$I / ${#GAPLENGTHS[@]}"

  END=$(($START + $GAPLENGTH + 41))
  REGION="$CONTIG:$START-$END"

  for ((i=0;i<${#MEANS[@]};++i)); do
    # Overlapping reads
    if [ ! -f overlap."$GAPLENGTH"."${MEANS[i]}".fa ]; then
      $SAMTOOLS view -u $DATA/aln."${MEANS[i]}".bam $REGION | \
        $SAMTOOLS fasta - > overlap."$GAPLENGTH"."${MEANS[i]}".fa
    fi

    # GapFiller-style filtering
    if [ ! -f gapfiller."$GAPLENGTH"."${MEANS[i]}".fa ]; then
      python3 $SCRIPTS/gapfiller_filter.py $DATA/aln."${MEANS[i]}".bam \
        $REGION $((${MEANS[i]} / 4)) \
        > gapfiller."$GAPLENGTH"."${MEANS[i]}".fa
    fi

    # Filter reads
    if [ ! -f filter."$GAPLENGTH"."${MEANS[i]}".fa ]; then
      $EXTRACT -bam $DATA/aln."${MEANS[i]}".bam -insertion \
        -read-length $READLENGTH -mean ${MEANS[i]} -std-dev ${STDDEVS[i]} \
        -region $REGION -reads filter."$GAPLENGTH"."${MEANS[i]}".fa
    fi

    # Add unmapped reads to filtered
    if [ ! -f filter2."$GAPLENGTH"."${MEANS[i]}".fa ]; then
      cp filter."$GAPLENGTH"."${MEANS[i]}".fa filter2."$GAPLENGTH"."${MEANS[i]}".fa
      FILTERED=$(grep '^[^>;]' filter."$GAPLENGTH"."${MEANS[i]}".fa | wc -c)
      if (( $(echo "($FILTERED / ($GAPLENGTH + 82) ) < $THRESHOLD" | bc -l) )); then
        cat unmapped."${MEANS[i]}".fa >> filter2."$GAPLENGTH"."${MEANS[i]}".fa
      fi
    fi

    # Evaluate the schemes
    RESULTS=$(python3 $SCRIPTS/evaluate.py \
      aln."${MEANS[i]}".fa \
      known."$GAPLENGTH"."${MEANS[i]}".fa \
      unmapped."${MEANS[i]}".fa \
      overlap."$GAPLENGTH"."${MEANS[i]}".fa \
      gapfiller."$GAPLENGTH"."${MEANS[i]}".fa \
      filter."$GAPLENGTH"."${MEANS[i]}".fa \
      filter2."$GAPLENGTH"."${MEANS[i]}".fa)

    echo $GAPLENGTH ${MEANS[i]} ${STDDEVS[i]} $RESULTS >> results
  done
done < $DATA/breakpoints.bed
