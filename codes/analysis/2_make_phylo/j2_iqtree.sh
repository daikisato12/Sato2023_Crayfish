#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 10

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Mapped_reads/mpileup/phylip
cd $WORK_DIR

FILE=pooled_7pops

iqtree -s $FILE.min4.phy -bb 1000 -nt 10

echo ending at
date
