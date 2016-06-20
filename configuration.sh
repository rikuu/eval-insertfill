#!/bin/bash

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)

# Required tools
DWGSIM=/cs/work/scratch/riqwalve/DWGSIM/dwgsim
BEDTOOLS=/cs/work/scratch/riqwalve/bedtools2/bin/bedtools
SAMTOOLS=~/samtools/samtools
BWA=~/bwa/bwa

GAP2SEQ=/cs/work/scratch/riqwalve/Gap2Seq/build/Gap2Seq
GAPCUTTER=/cs/work/scratch/riqwalve/Gap2Seq/build/GapCutter

EXTRACT=/cs/work/scratch/riqwalve/extract/extract

PINDEL=/cs/work/scratch/riqwalve/pindel/pindel
PINDEL2VCF=/cs/work/scratch/riqwalve/pindel/pindel2vcf
MINDTHEGAP=cs/work/scratch/riqwalve/MindTheGap/build/bin/MindTheGap

THREADS=16

# Parameters for reads
READLENGTH=100
COVERAGE=30
MEANS=(150 1500 3000)
STDDEVS=(15 150 300)

# TODO: Run fix_chr17.py ?
GENOME=$DIR/chr17.fa
CONTIG=chr17

PREFIX=$DIR
DATA=$PREFIX/data
