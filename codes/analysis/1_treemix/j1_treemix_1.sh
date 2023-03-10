#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Mapped_reads/mpileup
cd $WORK_DIR

FILE=pooled_7pops

vcftools --gzvcf $FILE.vcf.gz --max-missing 1 --recode --stdout | gzip > $FILE.noN.vcf.gz
/home/daikisato177/Projects/Crayfish/Build/ldPruning.sh $FILE.noN.vcf.gz
gzip $FILE.noN.LDpruned.vcf

echo ending at
date
