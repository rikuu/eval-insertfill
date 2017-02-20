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
  $TIME -v $GAP2SEQ -filled tmp.filled.all.normal -scaffolds $DATA/gaps.fa \
    -reads $READS -nb-cores $THREADS -max-mem $MAXMEM \
    2> gap2seq.stderr 1> gap2seq.stdout
fi

if [ ! -f tmp.filled.all.filter ]; then
  echo -e "Running Gap2Seq (filter-all)"
  $TIME -v python3 $SCRIPTS/filler.py -l $DATA/libraries.txt -g $DATA/gaps.fa \
    -b $DATA/breakpoints.bed -o tmp.filled.all.filter -t $THREADS \
    --max-mem $MAXMEM 2> gap2seq.filter.stderr 1> gap2seq.filter.stdout
fi

if [ ! -f tmp.filled.pindel ]; then
  $TIME -v $PINDEL -f $DATA/reference.fa -i $DATA/pindel.txt -o pindel \
    2> pindel.stderr 1> pindel.stdout

  $PINDEL2VCF -p pindel_LI -r $DATA/reference.fa -R chr17_1 -d 20160612
  python3 $SCRIPTS/vcf2filled.py $DATA/inserts.bed pindel_LI.vcf > tmp.filled.pindel

  $PINDEL2VCF -p pindel_SI -r $DATA/reference.fa -R chr17_1 -d 20160612
  python3 $SCRIPTS/vcf2filled.py $DATA/inserts.bed pindel_SI.vcf >> tmp.filled.pindel
fi

if [ ! -f tmp.filled.mtg ]; then
  #$MINDTHEGAP find -in $READS -ref $DATA/reference.fa -out mtg
  $TIME -v $MINDTHEGAP fill -in $READS -bkpt $DATA/mtg.gaps.fa -out mtg \
    2> mindthegap.stderr 1> mindthegap.stdout
  python3 $SCRIPTS/mtg2filled.py $DATA/gaps.fa mtg.insertions.fasta > tmp.filled.mtg
fi

python3 $SCRIPTS/evaluate_fill.py $DATA/inserts.fa \
  tmp.filled.all.normal tmp.filled.all.filter \
  tmp.filled.pindel tmp.filled.mtg > results
