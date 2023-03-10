#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

WORK_DIR=${data}/${mapped_reads_pool}/mpileup/HAFpipe
cd $WORK_DIR
REF=${data}/${ref}/GCA_020424385.2_ASM2042438v2_genomic_index_changed.fa
perl j0_mkupper.pl $REF ./GCA_020424385.1_ASM2042438v2_genomic_index_changed_upper.fa

echo ending at
date
