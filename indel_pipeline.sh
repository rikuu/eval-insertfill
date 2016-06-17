#!/bin/bash
set -e

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
SCRIPTS=$DIR/scripts

source $DIR/configuration.sh
OUT=$PREFIX/indel/

source $SCRIPTS/gen_reads.sh
source $SCRIPTS/gen_reference.sh

cd $OUT

for ((i=0;i<${#MEANS[@]};++i)); do
  if [ ! -f tmp.filled."${MEANS[i]}".normal ]; then
    echo -e "Running Gap2Seq (normal-${MEANS[i]})"
    $GAP2SEQ -filled tmp.filled."${MEANS[i]}".normal -scaffolds $DATA/gaps.fa \
      -reads $DATA/reads"$i"_pe1.fq,$DATA/reads"$i"_pe2.fq
  fi

  if [ ! -f tmp.filled."${MEANS[i]}".filter ]; then
    echo -e "Running Gap2Seq (filter-${MEANS[i]})"
    python $SCRIPTS/filler.py -l $DATA/libraries.txt -g $DATA/gaps.fa \
      -bed $DATA/gaps.bed -i "$i" -o tmp.filled."${MEANS[i]}".filter
  fi
done

if [ ! -f tmp.filled.all.filter ]; then
  # TODO: Don't hardcode these, figure out how not to
  echo -e "Running Gap2Seq (normal-all)"
  $GAP2SEQ -filled tmp.filled."${MEANS[i]}".all -scaffolds $DATA/gaps.fa \
    -reads $DATA/reads0_pe1.fq,$DATA/reads0_pe2.fq,$DATA/reads1_pe1.fq,$DATA/reads1_pe2.fq,$DATA/reads2_pe1.fq,$DATA/reads2_pe2.fq
fi

if [ ! -f tmp.filled.all.filter ]; then
  echo -e "Running Gap2Seq (filter-all)"
  python $SCRIPTS/filler.py -l libraries.txt -g $DATA/gaps.fa \
    -bed $DATA/gaps.bed -o tmp.filled.all.filter
fi

python $SCRIPTS/evaluate_fill.py $DATA/inserts.fa \
  tmp.filled.150.normal tmp.filled.150.filter \
  tmp.filled.1500.normal tmp.filled.1500.filter \
  tmp.filled.3000.normal tmp.filled.3000.filter \
  tmp.filled.all.normal tmp.filled.all.filter > results_fills

if [ ! -f tmp.filled.pindel ]; then
  $PINDEL -f reference.fa -i pindel.txt -o pindel
  $PINDEL2VCF -p pindel_LI -r reference.fa -R chr17_1 -d 20160612
  python $SCRIPTS/vcf2filled.py inserts.bed pindel_LI.vcf > tmp.filled.pindel
fi

python $SCRIPTS/evaluate_fill.py $DATA/inserts.fa \
  tmp.filled.all.normal tmp.filled.all.filter tmp.filled.pindel > results_tools
