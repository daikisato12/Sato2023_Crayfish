#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


WORKDIR=${data}/mapped_reads/mpileup/phylip
cd $WORKDIR

POOLNAME="pooled_7pops"
python /codes/software/vcf2phylip/vcf2phylip.py --input ../${POOLNAME}.vcf.gz

echo ending at
date
