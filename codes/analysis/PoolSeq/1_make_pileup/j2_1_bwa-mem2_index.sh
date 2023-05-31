#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

bwa-mem2 index ${data}/ref/GCA_020424385.2_ASM2042438v2_genomic_index_changed.fa

echo ending at
date
