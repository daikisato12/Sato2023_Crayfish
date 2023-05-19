#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


WORKDIR=${data}/mapped_reads/mpileup/phylip
cd $WORKDIR

POOLNAME="pooled_8pops"
iqtree -s ${POOLNAME}.min4.phy -bb 1000 -nt 10

echo ending at
date
