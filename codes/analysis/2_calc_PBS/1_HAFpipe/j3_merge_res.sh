#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

WORK_DIR=${data}/${mapped_reads_pool}/mpileup/HAFpipe
cd $WORK_DIR

perl j3_merge_res.pl

echo ending at
date
