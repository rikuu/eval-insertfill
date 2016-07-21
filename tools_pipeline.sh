#!/bin/bash
set -e

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
SCRIPTS=$DIR/scripts

source $DIR/configuration.sh
OUT=$PREFIX/tools/

source $SCRIPTS/gen_reads.sh
source $SCRIPTS/gen_reference.sh

mkdir -p $OUT || true
cd $OUT

READS="$DATA/reads0_pe1.fq,$DATA/reads0_pe2.fq"
for ((i=1;i<${#MEANS[@]};++i)); do
  READS="$READS,$DATA/reads"$i"_pe1.fq,$DATA/reads"$i"_pe2.fq"
done

if [ ! -f tmp.filled.all.normal ]; then
  echo -e "Running Gap2Seq (normal-all)"
  $GAP2SEQ -filled tmp.filled.all.normal -scaffolds $DATA/gaps.fa \
    -reads $READS -nb-cores $THREADS
fi

if [ ! -f tmp.filled.all.filter ]; then
  echo -e "Running Gap2Seq (filter-all)"
  python $SCRIPTS/filler.py -l $DATA/libraries.txt -g $DATA/gaps.fa \
    -b $DATA/breakpoints.bed -o tmp.filled.all.filter -t $THREADS
fi

if [ ! -f tmp.filled.pindel ]; then
  $PINDEL -f $DATA/reference.fa -i $DATA/pindel.txt -o pindel

  $PINDEL2VCF -p pindel_LI -r $DATA/reference.fa -R chr17_1 -d 20160612
  python $SCRIPTS/vcf2filled.py $DATA/inserts.bed pindel_LI.vcf > tmp.filled.pindel

  $PINDEL2VCF -p pindel_SI -r $DATA/reference.fa -R chr17_1 -d 20160612
  python $SCRIPTS/vcf2filled.py $DATA/inserts.bed pindel_SI.vcf >> tmp.filled.pindel
fi

if [ ! -f tmp.filled.mtg ]; then
  #$MINDTHEGAP find -in $READS -ref $DATA/reference.fa -out mtg
  $MINDTHEGAP fill -in $READS -bkpt $DATA/mtg.gaps.fa -out mtg
  python $SCRIPTS/mtg2filled.py $DATA/gaps.fa mtg.insertions.fa > tmp.filled.mtg
fi

python $SCRIPTS/evaluate_fill.py $DATA/inserts.fa \
  tmp.filled.all.normal tmp.filled.all.filter \
  tmp.filled.pindel tmp.filled.mtg > results_tools
