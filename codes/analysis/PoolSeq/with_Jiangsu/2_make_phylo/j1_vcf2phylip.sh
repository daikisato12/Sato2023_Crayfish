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
python /codes/software/vcf2phylip/vcf2phylip.py --input ../${POOLNAME}.vcf

echo ending at
date
