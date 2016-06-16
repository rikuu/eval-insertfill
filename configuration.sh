#!/bin/bash

# Required tools
DWGSIM=dwgsim
BEDTOOLS=bedtools
SAMTOOLS=samtools
BWA=bwa
GAP2SEQ=/cs/work/scratch/riqwalve/Gap2Seq-2.0/build/Gap2Seq
GAPCUTTER=/cs/work/scratch/riqwalve/Gap2Seq-2.0/build/GapCutter
EXTRACT=/cs/work/scratch/riqwalve/extract/extract
PINDEL=/cs/work/scratch/riqwalve/pindel/pindel

# Parameters for reads
READLENGTH=100
COVERAGE=30
MEANS=(150 1500 3000)
STDDEVS=(15 150 300)
CONTIG=chr17
GENOME=chr17.fa
# TODO: Run fix_chr17.py ?

PREFIX=$DIR
#DATA=$PREFIX/data/
