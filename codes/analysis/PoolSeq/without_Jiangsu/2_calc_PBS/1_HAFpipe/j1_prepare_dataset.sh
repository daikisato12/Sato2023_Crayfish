#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

WORKDIR=${data}/mapped_reads/mpileup/hafpipe
cd $WORKDIR

REF=${data}/ref/GCA_020424385.2_ASM2042438v2_genomic_index_changed.fa
perl /codes/analysis/PoolSeq/2_calc_PBS/1_HAFpipe/j0_mkupper.pl $REF ./GCA_020424385.1_ASM2042438v2_genomic_index_changed_upper.fa

# subseq contigs larger than 500kbp
seqkit seq -m 500000 $REF > Ref_500k.fa
perl -lane 'if($_=~/>(\S+)/){print $1}' Ref_500k.fa > contig_500kbp.txt

echo ending at
date
