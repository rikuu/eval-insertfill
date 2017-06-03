#!/bin/bash
set -e

################################################################################

REFERENCE=c_elegans.WS210.genomic.fa
READS_SPACE=SRR1259172_1.fastq.gz SRR1259172_2.fastq.gz
READS_COMMA=SRR1259172_1.fastq.gz,SRR1259172_2.fastq.gz
BAM=CB4856.bam
THREADS=16

################################################################################

TIME=/usr/bin/time
BWA=~/bwa/bwa
SAMTOOLS=~/samtools/samtools

################################################################################

BREAKDANCER=~/breakdancer-1.1.2/cpp/breakdancer-max
BAM2CFG=~/breakdancer-1.1.2/perl/bam2cfg

PINDEL=~/pindel/pindel
PINDEL2VCF=~/pindel/pindel2vcf

################################################################################

SCRIPTS=~/eval-insertfill/scripts
GAP2SEQ=~/Gap2Seq/build/Gap2Seq
MINDTHEGAP=~/MindTheGap/build/bin/MindTheGap
GAPFILLER=~/gapfiller-1.10/GapFiller.pl
GAPCLOSER=~/GapCloser/GapCloser
SEALER=~/abyss-1.9.0/Sealer/abyss-sealer

################################################################################

# Align donor reads to reference
$SAMTOOLS faidx $REFERENCE
$BWA index $REFERENCE
$TIME -v $BWA -t $THREADS $REFERENCE $READS_SPACE | $SAMTOOLS sort - | \
  $SAMTOOLS view -bh > $BAM 2> bwa.stderr

# Run breadancer-max to find SVs
perl $BAM2CFG $BAM > breakdancer.cfg
$TIME -v $BREAKDANCER -y 15 -q 15 breakdancer.cfg 1> breakdancer.sv 2> breakdancer.stderr

# Run pindel (with breakdancer SVs) to find more SVs
MEAN=$(cat breakdancer.cfg | cut -f9 | cut -c6-)
echo -e "$BAM\t$MEAN\tCB4856" > pindel.cfg
$TIME -v $PINDEL -f $REFERENCE -i pindel.cfg -o pindel -T $THREADS \
  --report_inversions false --report_duplications false --report_long_insertions true \
  --breakdancer breakdancer.sv --genotyping 2> pindel.stderr

$PINDEL2VCF -p pindel_SI -r $REFERENCE -R CB4856 -d 20170307
python3 $SCRIPTS/vcf2gaps.py $REFERENCE pindel_SI.vcf pindel_SI.gaps.fa pindel_SI.filled

################################################################################

# Gap2Seq
$TIME -v $GAP2SEQ -nb-cores $THREADS -reads $READS_COMMA \
  -filled filled.gap2seq -scaffolds pindel_SI.gaps.fa 2> gap2seq.stderr

# MindTheGap
python3 $SCRIPTS/gaps2mtg.py pindel_SI.gaps.fa > pindel_SI.mtg.fa
$TIME -v $MINDTHEGAP fill -nb-cores $THREADS -in $READS_COMMA \
  -bkpt pindel_SI.mtg.fa -out mtg 2> mindthegap.stderr
python3 $SCRIPTS/mtg2filled.py pindel_SI.gaps.fa mtg.insertions.fasta > filled.mindthegap

# GapFiller
$TIME -v $GAPFILLER -T $THREADS -s pindel_SI.gaps.fa -l GapFiller.config \
  -b gapfiller 2> gapfiller.stderr
mv gapfiller/gapfiller.gapfilled.final.fa filled.gapfiller
rm -r gapfiller

# GapCloser
$TIME -v $GAPCLOSER -a pindel_SI.gaps.fa -b GapCloser.config \
  -o filled.gapcloser -t $THREADS 2> gapcloser.stderr

# Sealer
$TIME -v $SEALER -j $THREDS -o sealer -S pindel_SI.gaps.fa \
  -k34 -k33 -k32 -k31 -k30 -k29 -k28 -P 5 $READS_SPACE 2> sealer.stderr
mv sealer_scaffold.fa filled.sealer

################################################################################

python3 $SCRIPTS/evaluate_fill.py pindel_SI.filled \
  filled.gap2seq filled.gap2seq.filter filled.mindthegap filled.gapfiller filled.gapcloser filled.sealer > filled

python3 $SCRIPTS/venn.py pindel_SI.vcf \
  filled.gap2seq filled.gap2seq.filter filled.mindthegap filled.gapfiller filled.gapcloser filled.sealer > venn
