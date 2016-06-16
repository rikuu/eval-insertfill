#!/usr/local/bin/bash

# Required tools
DWGSIM=dwgsim
BEDTOOLS=bedtools
SAMTOOLS=samtools
BWA=bwa
GAP2SEQ=../Gap2Seq-2.0/build/Gap2Seq
GAPCUTTER=../Gap2Seq-2.0/build/GapCutter
EXTRACT=../extract/extract

# Parameters for the testing
READLENGTH=100
COVERAGE=30
CONTIG=chr17
MEANS=(150 1500 3000)
STDDEVS=(15 150 300)

# '('+' '.join([str(int(i)) for i in np.logspace(log(11, 10), log(5000, 10), 20)])+')'
GAPLENGTHS=(1 10 15 20 28 39 55 75 104 144 199 275 380 524 723 \
  999 1378 1902 2625 3623 4999)

# Exit if any command fails
set -e

# Index the reference genome
if [ ! -f chr17.fa.amb ]; then
  $BWA index chr17.fa
fi

for ((i=0;i<${#MEANS[@]};++i)); do
  # Generate reads with the different parameters
  if [ ! -f reads"$i"_pe1.fq ]; then
    $DWGSIM -i -1 $READLENGTH -2 $READLENGTH -d ${MEANS[i]} -s ${STDDEVS[i]} \
      -C $COVERAGE chr17.fa reads"$i"
    rm reads"$i".bfast*
    # rm reads"$i".mutations.vcf

    mv reads"$i".bwa.read1.fastq reads"$i"_pe1.fq
    mv reads"$i".bwa.read2.fastq reads"$i"_pe2.fq
  fi

  # Map-sort-index the reads
  if [ ! -f known_aln"$i".bam ]; then
    $BWA mem -t 16 -I ${MEANS[i]},${STDDEVS[i]} chr17.fa reads"$i"_pe1.fq \
        reads"$i"_pe2.fq | \
      $SAMTOOLS view -Shu - | \
      $SAMTOOLS sort - | \
      $SAMTOOLS rmdup -s - - > known_aln"$i".bam
    $SAMTOOLS index known_aln"$i".bam
  fi
done

for GAPLENGTH in "${GAPLENGTHS[@]}"; do
  rm -f tmp.* gapped_genome.fa*

  # Generate a gapped genome
  START=$(((RANDOM % 80000) + 1000))
  END=$((START+GAPLENGTH))
  echo -e "$CONTIG\t$START\t$END" > tmp.bed
  $BEDTOOLS maskfasta -fi chr17.fa -bed tmp.bed -fo gapped_genome.fa
  $BWA index gapped_genome.fa

  # Create a bed file (simulates a real pipeline, I guess)
  $GAPCUTTER -bed tmp.bed -gaps tmp.gaps -contigs tmp.contigs \
    -scaffolds gapped_genome.fa
  CONTIG=$(cat tmp.bed | cut -f1)
  START=$(cat tmp.bed | cut -f2)
  END=$(cat tmp.bed | cut -f3)

  for ((i=0;i<${#MEANS[@]};++i)); do
    rm -f pe*.bam aln.fa known.fa overlap.fa unmapped.fa filter.fa

    # Map-sort-index reads to the gapped genome/assembly
    for PE in pe1 pe2; do
      $BWA mem -t 16 -I ${MEANS[i]},${STDDEVS[i]} gapped_genome.fa \
          reads"$i"_"$PE".fq | \
        $SAMTOOLS view -Shu - | \
        $SAMTOOLS sort - | \
        $SAMTOOLS rmdup -s - - > "$PE".bam
      $SAMTOOLS index "$PE".bam

      # Extract all reads
      $SAMTOOLS fasta "$PE".bam >> aln.fa

      # Extract all overlapping reads
      $SAMTOOLS view -F4 -b -L tmp.bed "$PE".bam | \
        $SAMTOOLS fasta - >> overlap.fa

      # Extract all unmapped reads
      $SAMTOOLS view -f4 -b "$PE".bam | \
        $SAMTOOLS fasta - >> unmapped.fa
    done

    # Extract all known gap-covering reads
    $SAMTOOLS view -b -L tmp.bed known_aln"$i".bam | \
      $SAMTOOLS fasta - > known.fa

    # Extract filtered reads
    $EXTRACT pe1.bam pe2.bam $READLENGTH ${MEANS[i]} ${STDDEVS[i]} \
      $CONTIG $START $END > filter.fa

    # Evaluate the schemes
    OVERLAP=$(python evaluate.py aln.fa known.fa overlap.fa)
    UNMAPPED=$(python evaluate.py aln.fa known.fa unmapped.fa)
    FILTER=$(python evaluate.py aln.fa known.fa filter.fa)
    echo $GAPLENGTH ${MEANS[i]} ${STDDEVS[i]} $OVERLAP $UNMAPPED $FILTER >> results_npe
  done
done

rm -f aln.* tmp.* known.* overlap.* unmapped.* filter.* gapped_genome.fa
