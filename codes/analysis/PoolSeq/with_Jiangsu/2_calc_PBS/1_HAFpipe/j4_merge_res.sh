#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

WORK_DIR=${data}/poolseq/mapped_reads/mpileup/hafpipe
cd $WORK_DIR

perl j4_merge_res.pl

echo ending at
date
