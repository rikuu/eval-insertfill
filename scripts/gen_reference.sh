#!/bin/bash
set -e

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)/..
source $DIR/configuration.sh

cd $DATA

if [ ! -f inserts.bed ]; then
  python $SCRIPTS/gap_bed.py inserts.bed gaps.bed breakpoints.bed
fi

if [ ! -f reference.fa ]; then
  # Make sure no indexes are left from previous runs
  rm -f assembly.fa* gaps.fa* inserts.fa* reference.fa*

  # Create a fake assembly (masked donor genome in SV context)
  $BEDTOOLS maskfasta -fi $GENOME -bed inserts.bed -fo assembly.fa

  # Cut the inserts and gaps from the assembly
  $BEDTOOLS getfasta -fi assembly.fa -bed gaps.bed -fo gaps.fa
  $BEDTOOLS getfasta -fi $GENOME -bed gaps.bed -fo inserts.fa

  # Remove all inserts from the masked donor to create a reference genome
  python $SCRIPTS/reference.py assembly.fa reference.fa
  $BWA index reference.fa
  $SAMTOOLS faidx reference.fa
fi

rm -f pindel.txt libraries.txt
for ((i=0;i<${#MEANS[@]};++i)); do
  if [ ! -f aln."${MEANS[i]}".bam ]; then
    # Map-sort-index reads to the reference genome
    echo -e "Aligning reads (${MEANS[i]})"
    $BWA mem -t $THREADS -I ${MEANS[i]},${STDDEVS[i]} reference.fa \
        reads"$i"_pe1.fq reads"$i"_pe2.fq | \
      $SAMTOOLS view -Shu - | \
      $SAMTOOLS sort - | \
      $SAMTOOLS rmdup -s - - > aln."${MEANS[i]}".bam
    $SAMTOOLS index aln."${MEANS[i]}".bam
  fi

  echo -e "$DATA/aln.${MEANS[i]}.bam\t${MEANS[i]}\tALN${MEANS[i]}" >> pindel.txt
  echo -e "$DATA/aln.${MEANS[i]}.bam\t$READLENGTH\t${MEANS[i]}\t${STDDEVS[i]}" >> libraries.txt
done
