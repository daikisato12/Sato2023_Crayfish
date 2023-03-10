#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Mapped_reads/mpileup/phylip
cd $WORK_DIR

FILE=pooled_7pops

python /home/daikisato177/Projects/Crayfish/Build/vcf2phylip/vcf2phylip.py --input ../$FILE.vcf.gz

echo ending at
date
