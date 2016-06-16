#!/bin/bash
set -e

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
SCRIPTS=$DIR/scripts

source $DIR/configuration.sh

DATA=$PREFIX/data/
OUT=$PREFIX/indel/

source $SCRIPTS/gen_reads.sh
source $SCRIPTS/gen_reference.sh

cd $OUT

for ((i=0;i<${#MEANS[@]};++i)); do
  echo -e "Running Gap2Seq (normal-${MEANS[i]})"
  $GAP2SEQ -filled tmp.filled."${MEANS[i]}".each -scaffolds $DATA/gaps.fa \
    -reads $DATA/reads"$i"_pe1.fq,$DATA/reads"$i"_pe2.fq

  echo -e "Running Gap2Seq (filter-${MEANS[i]})"
  python $SCRIPTS/filler.py -l $DATA/libraries.txt -g $DATA/gaps.fa -bed $DATA/gaps.bed -i "$i"
done

# TODO: Don't hardcode these, figure out how not to
echo -e "Running Gap2Seq (normal-all)"
$GAP2SEQ -filled tmp.filled."${MEANS[i]}".all -scaffolds $DATA/gaps.fa \
  -reads $DATA/reads0_pe1.fq,$DATA/reads0_pe2.fq,$DATA/reads1_pe1.fq,$DATA/reads1_pe2.fq,$DATA/reads2_pe1.fq,$DATA/reads2_pe2.fq

echo -e "Running Gap2Seq (filter-all)"
python $SCRIPTS/filler.py -l libraries.txt -g $DATA/gaps.fa -bed $DATA/gaps.bed

python $SCRIPTS/evaluate_fill.py $DATA/inserts.fa \
  tmp.filled.150.normal tmp.filled.150.filter \
  tmp.filled.1500.normal tmp.filled.1500.filter \
  tmp.filled.3000.normal tmp.filled.3000.filter \
  tmp.filled.all.normal tmp.filled.all.filter > results_fills

# $PINDEL -f reference.fa -i pindel.txt -o pindel
