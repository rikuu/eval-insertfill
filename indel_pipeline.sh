#!/bin/bash

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
source $DIR/configuration.sh

source $DIR/gen_reads.sh
source $DIR/gen_reference.sh

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
