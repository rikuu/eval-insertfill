#!/bin/bash

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

# Exit if any command fails
set -e
