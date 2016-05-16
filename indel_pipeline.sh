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

# '('+' '.join([str(int(i)) for i in np.logspace(log(11, 10), log(5000, 10), 20)])+')'
GAPLENGTHS=(1 10 15 20 28 39 55 75 104 144 199 275 380 524 723 \
  999 1378 1902 2625 3623 4999)

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

for LENGTH in "${LENGTHS[@]}"; do
  rm -f assembly.fa tmp.*

  # Create a fake assembly/reference genome
  START=$(((RANDOM % 80000) + 1000))
  END=$((START+LENGTH))
  echo -e "$CONTIG\t$START\t$END" > tmp.bed
  $BEDTOOLS maskfasta -fi chr17.fa -bed tmp.bed -fo assembly.fa
  $BWA index assembly.fa

  # Cut the gap from the assembly
  echo -e "chr17\t$((START-41))\t$((END+41))" > tmp.bed
  $BEDTOOLS getfasta -fi assembly.fa -bed tmp.bed -fo gap.fa

  # Cut the insert from the assembly
  echo -e "chr17\t$((START-41))\t$((END+41))" > tmp.bed
  $BEDTOOLS getfasta -fi chr17.fa -bed tmp.bed -fo insert.fa

  for ((i=0;i<${#MEANS[@]};++i)); do
    rm -f filter.fa aln.bam* tmp.filled.*

    # Map-sort-index reads to the genome/assembly
    $BWA mem -t 16 -I ${MEANS[i]},${STDDEVS[i]} assembly.fa \
        reads"$i"_pe1.fq reads"$i"_pe2.fq | \
      $SAMTOOLS view -Shu - | \
      $SAMTOOLS sort - | \
      $SAMTOOLS rmdup -s - - > aln.bam
    $SAMTOOLS index aln.bam

    # Extract filtered reads
    $EXTRACT aln.bam $READLENGTH ${MEANS[i]} ${STDDEVS[i]} $CONTIG $START $END \
      > filter.fa

    # Genotype the insert
    $GAP2SEQ -filled tmp.filled."${MEANS[i]}" -scaffolds gap.fa \
      -reads reads"$i".fq # -graph graph.h5
    NORMAL=$(python evaluate_fill.py insert.fa tmp.filled."${MEANS[i]}")

    # Genotype with filtered reads if any, Gap2Seq crashes with empty read set
    FILTER=0
    if [[ -s filter.fa ]] ; then
      $GAP2SEQ -filled tmp.filled."${MEANS[i]}".filter -scaffolds gap.fa \
        -reads filter.fa
      FILTER=$(python evaluate_fill.py insert.fa tmp.filled."${MEANS[i]}".filter)
    fi

    echo "$LENGTH ${MEANS[i]} ${STDDEVS[i]} $FILTER $NORMAL" >> results_fills
  done
done

rm -f aln.bam* tmp.* filter.fa assembly.fa
