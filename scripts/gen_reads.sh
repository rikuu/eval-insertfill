#!/bin/bash
set -e

source $DIR/configuration.sh

if [ ! -f $GENOME ]; then
  echo "Reference genome (set to \"$GENOME\") missing!"
  exit 1
fi

mkdir -p $DATA
cd $DATA

for ((i=0;i<${#MEANS[@]};++i)); do
  # Generate reads with the different parameters
  if [ ! -f reads"$i"_pe1.fq ]; then
    $ART -i $GENOME -l $READLENGTH -f $COVERAGE -m ${MEANS[i]} -s ${STDDEVS[i]} -sam -o reads"$i"

    $SAMTOOLS sort reads"$i".sam | \
      $SAMTOOLS view - -b --threads $THREADS > known_aln"$i".bam

    $SAMTOOLS index known_aln"$i".bam

    rm -f reads"$i"*.aln reads"$i".sam

    mv reads"$i"1.fq reads"$i"_pe1.fq
    mv reads"$i"2.fq reads"$i"_pe2.fq
  fi
done
