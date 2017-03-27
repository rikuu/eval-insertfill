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

READS_SPACE="$DATA/reads0_pe1.fq $DATA/reads0_pe2.fq"
READS_COMMA="$DATA/reads0_pe1.fq,$DATA/reads0_pe2.fq"
for ((i=1;i<${#MEANS[@]};++i)); do
  READS_SPACE="$READS_SPACE $DATA/reads"$i"_pe1.fq $DATA/reads"$i"_pe2.fq"
  READS_COMMA="$READS_COMMA,$DATA/reads"$i"_pe1.fq,$DATA/reads"$i"_pe2.fq"
done

if [ ! -f tmp.filled.all.normal ]; then
  echo -e "Running Gap2Seq (normal-all)"
  $TIME -v $GAP2SEQ -filled tmp.filled.all.normal -scaffolds $DATA/gaps.fa \
    -reads $READS_COMMA -nb-cores $THREADS -max-mem $MAXMEM \
    2> gap2seq.stderr 1> gap2seq.stdout
fi

if [ ! -f tmp.filled.all.filter ]; then
  echo -e "Running Gap2Seq (filter-all)"
  $TIME -v python3 $SCRIPTS/filler.py -l $DATA/libraries.txt -g $DATA/gaps.fa \
    -b $DATA/breakpoints.bed -o tmp.filled.all.filter -t $THREADS \
    --max-mem $MAXMEM 2> gap2seq.filter.stderr 1> gap2seq.filter.stdout
fi

if [ ! -f tmp.filled.pindel ]; then
  echo -e "Running Pindel (all)"
  $TIME -v $PINDEL -f $DATA/reference.fa -i $DATA/pindel.txt -o pindel \
    -T $THREADS --report_inversions false --report_duplications false \
    --report_long_insertions true --genotyping 2> pindel.stderr 1> pindel.stdout

  $PINDEL2VCF -p pindel_LI -r $DATA/reference.fa -R chr17_1 -d 20160612
  python3 $SCRIPTS/vcf2filled.py $DATA/inserts.bed pindel_LI.vcf > tmp.filled.pindel

  $PINDEL2VCF -p pindel_SI -r $DATA/reference.fa -R chr17_1 -d 20160612
  python3 $SCRIPTS/vcf2filled.py $DATA/inserts.bed pindel_SI.vcf >> tmp.filled.pindel
fi

if [ ! -f tmp.filled.mtg ]; then
  echo -e "Running MindTheGap (all)"
  #$MINDTHEGAP find -in $READS -ref $DATA/reference.fa -out mtg
  $TIME -v $MINDTHEGAP fill -in $READS_COMMA -bkpt $DATA/mtg.gaps.fa -out mtg \
    2> mindthegap.stderr 1> mindthegap.stdout

  python3 $SCRIPTS/mtg2filled.py $DATA/gaps.fa mtg.insertions.fasta > tmp.filled.mtg
fi

if [ ! -f tmp.filled.gapfiller ]; then
  echo -e "Running GapFiller (all)"
  $TIME -v $GAPFILLER -T $THREADS -s $DATA/gaps.fa -l $DIR/GapFiller.config \
    -b gapfiller 2> gapfiller.stderr 1> gapfiller.stdout

  mv gapfiller/gapfiller.gapfilled.final.fa tmp.filled.gapfiller
  rm -r gapfiller
fi

if [ ! -f tmp.filled.gapcloser ]; then
  echo -e "Running GapCloser (all)"
  $TIME -v $GAPCLOSER -a $DATA/gaps.fa -b $DIR/GapCloser.config \
    -o tmp.filled.gapcloser -t $THREADS 2> gapcloser.stderr 1> gapcloser.stdout
fi

if [ ! -f tmp.filled.sealer ]; then
  echo -e "Running Sealer (all)"
  $TIME -v $SEALER -j $THREADS -o sealer -S $DATA/gaps.fa \
    -k34 -k33 -k32 -k31 -k30 -k29 -k28 -P 5 $READS_SPACE \
    2> sealer.stderr 1> sealer.stdout

  mv sealer_scaffold.fa tmp.filled.sealer
fi

# Gather wall clock times for all tools
tail -n19 $DATA/align.150.stderr | head -n1 | cut -c47- > times
tail -n19 $DATA/align.1500.stderr | head -n1 | cut -c47- >> times
tail -n19 $DATA/align.3000.stderr | head -n1 | cut -c47- >> times
tail -n19 gap2seq.stderr | head -n1 | cut -c47- >> times
tail -n19 gap2seq.filter.stderr | head -n1 | cut -c47- >> times
tail -n19 pindel.stderr | head -n1 | cut -c47- >> times
tail -n19 mindthegap.stderr | head -n1 | cut -c47- >> times
tail -n19 gapfiller.stderr | head -n1 | cut -c47- >> times
tail -n19 gapcloser.stderr | head -n1 | cut -c47- >> times
tail -n19 sealer.stderr | head -n1 | cut -c47- >> times

echo -e "Evaluating"
python3 $SCRIPTS/evaluate_fill.py $DATA/inserts.fa \
  tmp.filled.all.normal tmp.filled.all.filter \
  tmp.filled.pindel tmp.filled.mtg tmp.filled.gapfiller \
  tmp.filled.gapcloser tmp.filled.sealer > results
