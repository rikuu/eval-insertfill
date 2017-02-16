#!/bin/bash

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)

# Required tools
ART=~/art_bin_MountRainier/art_illumina
BEDTOOLS=/cs/work/scratch/riqwalve/bedtools2/bin/bedtools
SAMTOOLS=~/samtools/samtools
BWA=~/bwa/bwa
BOWTIE=~/bowtie2-2.3.0/bowtie2
BOWTIEBUILD=~/bowtie2-2.3.0/bowtie2-build

GAP2SEQ=$DIR/Gap2Seq/build/Gap2Seq
GAPCUTTER=$DIR/Gap2Seq/build/GapCutter

EXTRACT=$DIR/extract/build/extract

PINDEL=/cs/work/scratch/riqwalve/pindel/pindel
PINDEL2VCF=/cs/work/scratch/riqwalve/pindel/pindel2vcf
MINDTHEGAP=/cs/work/scratch/riqwalve/MindTheGap/build/bin/MindTheGap

# Number of threads used in all tools
THREADS=16

# Maximum memory usage for Gap2Seq (in Gb)
MAXMEM=20

# Read simulation parameters
READLENGTH=100
COVERAGE=30
MEANS=(150 1500 3000)
STDDEVS=(15 150 300)

# Gap simulation parameters
GAPNUM=30
MINGAPLEN=11
MAXGAPLEN=10000

# Threshold for using unmapped reads in filtering
# Should probably be close but less than $COVERAGE
THRESHOLD=25

# TODO: Run fix_chr17.py ?
GENOME=$DIR/chr17.fa
CONTIG=chr17

PREFIX=$DIR
DATA=$PREFIX/data
