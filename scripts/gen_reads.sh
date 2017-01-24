#!/bin/bash
set -e

source $DIR/configuration.sh

cd $DATA

# Index the source genome
if [ ! -f $GENOME.amb ]; then
  $BWA index $GENOME
fi

for ((i=0;i<${#MEANS[@]};++i)); do
  # Generate reads with the different parameters
  if [ ! -f reads"$i"_pe1.fq ]; then
    $DWGSIM -i -1 $READLENGTH -2 $READLENGTH -d ${MEANS[i]} -s ${STDDEVS[i]} \
      -C $COVERAGE $GENOME reads"$i"
    rm reads"$i".bfast*

    mv reads"$i".bwa.read1.fastq reads"$i"_pe1.fq
    mv reads"$i".bwa.read2.fastq reads"$i"_pe2.fq
  fi

  # Map-sort-index the reads
  if [ ! -f known_aln"$i".bam ]; then
    $BWA mem -t $THREADS -I ${MEANS[i]},${STDDEVS[i]} $GENOME reads"$i"_pe1.fq \
        reads"$i"_pe2.fq | \
      $SAMTOOLS view -Shu - | \
      $SAMTOOLS sort - | \
      $SAMTOOLS rmdup -s - - > known_aln"$i".bam
    $SAMTOOLS index known_aln"$i".bam
  fi
done
