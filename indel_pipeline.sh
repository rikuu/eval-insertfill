#!/usr/local/bin/bash

# Required tools
DWGSIM=dwgsim
BEDTOOLS=bedtools
SAMTOOLS=samtools
BWA=bwa
GAP2SEQ=../Gap2Seq-2.0/build/Gap2Seq
EXTRACT=../extract/extract

# Parameters for the testing
READLENGTH=100
COVERAGE=30
MEANS=(150 1500 3000)
STDDEVS=(15 150 300)
CONTIG=chr17

# Exit if any command fails
set -e

# Generate reads with the different parameters
for ((i=0;i<${#MEANS[@]};++i)); do
  if [ ! -f reads"$i"_pe1.fq ]; then
    $DWGSIM -i -1 $READLENGTH -2 $READLENGTH -d ${MEANS[i]} -s ${STDDEVS[i]} \
      -C $COVERAGE chr17.fa reads"$i"
    rm reads"$i".bfast*

    mv reads"$i".bwa.read1.fastq reads"$i"_pe1.fq
    mv reads"$i".bwa.read2.fastq reads"$i"_pe2.fq

    # Combine reads, Gap2Seq crashes with multiple reads (?!)
    cat reads"$i"_pe1.fq reads"$i"_pe2.fq > reads"$i".fq
  fi

  # TODO: Create DBG for the reads here
done

if [ ! -f gaps.fa ]; then
  # Create a fake assembly (masked donor genome in SV context)
  python gap_bed.py inserts.bed gaps.bed
  $BEDTOOLS maskfasta -fi chr17.fa -bed inserts.bed -fo assembly.fa

  # Cut the inserts and gaps from the assembly
  $BEDTOOLS getfasta -fi assembly.fa -bed gaps.bed -fo gaps.fa
  $BEDTOOLS getfasta -fi chr17.fa -bed gaps.bed -fo inserts.fa

  # Remove all inserts from the masked donor to create a reference genome
  python reference.py assembly.fa reference.fa
  $BWA index reference.fa

  rm -f pindel.txt libraries.txt aln.*.bam
  for ((i=0;i<${#MEANS[@]};++i)); do
    # Map-sort-index reads to the reference genome
    echo -e "Aligning reads (${MEANS[i]})"
    $BWA mem -t 16 -I ${MEANS[i]},${STDDEVS[i]} reference.fa \
        reads"$i"_pe1.fq reads"$i"_pe2.fq | \
      $SAMTOOLS view -Shu - | \
      $SAMTOOLS sort - | \
      $SAMTOOLS rmdup -s - - > aln."${MEANS[i]}".bam
    $SAMTOOLS index aln."${MEANS[i]}".bam

    # echo -e "aln.${MEANS[i]}.bam\t${MEANS[i]}\tALN${MEANS[i]}" >> pindel.txt
    echo -e "aln.${MEANS[i]}.bam\t$READLENGTH\t${MEANS[i]}\t${STDDEVS[i]}" >> libraries.txt
  done
done

for ((i=0;i<${#MEANS[@]};++i)); do
  echo -e "Running Gap2Seq (normal-${MEANS[i]})"
  $GAP2SEQ -filled tmp.filled."${MEANS[i]}".each -scaffolds gaps.fa \
    -reads reads"$i"_pe1.fq,reads"$i"_pe2.fq
done

# TODO: Don't hardcode these, figure out how not to
echo -e "Running Gap2Seq (normal-all)"
$GAP2SEQ -filled tmp.filled."${MEANS[i]}".all -scaffolds gaps.fa \
  -reads reads0_pe1.fq,reads0_pe2.fq,reads1_pe1.fq,reads1_pe2.fq,reads2_pe1.fq,reads2_pe2.fq

# Run filtered Gap2Seq
echo -e "Running filtered Gap2Seq"
python filler.py libraries.txt

python evaluate_fill.py inserts.fa \
  tmp.filled.150.normal tmp.filled.150.filter \
  tmp.filled.1500.normal tmp.filled.1500.filter \
  tmp.filled.3000.normal tmp.filled.3000.filter \
  tmp.filled.all.normal tmp.filled.all.filter > results_fills

# $PINDEL -f reference.fa -i pindel.txt -o pindel
