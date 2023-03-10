#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

cd /home/daikisato177/Projects/Crayfish/Data_2022/Ref_Xu2021
bwa-mem2 index ./GCA_020424385.2_ASM2042438v2_genomic_index_changed.fa

echo ending at
date
