#!/bin/bash

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)

############################ TOOLS ##########################################

# https://www.niehs.nih.gov/research/resources/software/biostatistics/art/
ART=~/art_bin_MountRainier/art_illumina

# http://bedtools.readthedocs.io/en/latest/
BEDTOOLS=~/bedtools2/bin/bedtools

# https://github.com/samtools/samtools
SAMTOOLS=~/samtools/samtools

# http://bio-bwa.sourceforge.net
BWA=~/bwa/bwa

# http://bowtie-bio.sourceforge.net/index.shtml
BOWTIE=~/bowtie2-2.3.0/bowtie2
BOWTIEBUILD=~/bowtie2-2.3.0/bowtie2-build

# https://github.com/genome/pindel
PINDEL=~/pindel/pindel
PINDEL2VCF=~/pindel/pindel2vcf

# https://github.com/GATB/MindTheGap
MINDTHEGAP=~/MindTheGap/build/bin/MindTheGap

GAPFILLER=~/gapfiller-1.10/GapFiller.pl
GAPCLOSER=~/GapCloser/GapCloser
SEALER=~/abyss-1.9.0/Sealer/abyss-sealer

# GNU time
TIME=/usr/bin/time

# Included as git submodule
GAP2SEQ=$DIR/Gap2Seq/build/Gap2Seq-core
GAPCUTTER=$DIR/Gap2Seq/build/GapCutter
EXTRACT=$DIR/Gap2Seq/build/ReadFilter
GAP2SEQFILTER=$DIR/Gap2Seq/build/Gap2Seq

############################ PARAMETERS ######################################

# Number of threads used in all tools
THREADS=16

# Maximum memory usage for Gap2Seq (in Gb)
MAXMEM=20

# Read simulation parameters
READLENGTH=100
COVERAGE=30
MEANS=(150 1500 3000)
STDDEVS=(15 150 300)

# Aligner, "bowtie" or "bwa"
ALIGN="bowtie"

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
