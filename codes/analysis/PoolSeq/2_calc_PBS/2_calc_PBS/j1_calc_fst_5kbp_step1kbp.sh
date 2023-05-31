#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

WORKDIR=${data}/mapped_reads/mpileup/hafpipe
cd $WORKDIR

perl /codes/analysis/PoolSeq/without_Jiangsu/2_calc_PBS/2_calc_PBS/1_calc_pi.pl merge_result.tsv merge_result_pi.tsv

INFILE=${WORKDIR}/merge_result_pi.tsv
OUTFILE1=${WORKDIR}/fst/HAF_fst_5kbp_step1kbp.tsv
OUTFILE2=${WORKDIR}/pbs/HAF_PBS_5kbp_step1kbp.tsv
WINDOW_SIZE=5000
STEP_SIZE=1000

perl /codes/analysis/PoolSeq/with_Jiangsu/2_calc_PBS/2_calc_PBS/2_calc_fst_sliding_window.pl $INFILE $OUTFILE1 $WINDOW_SIZE $STEP_SIZE
perl /codes/analysis/PoolSeq/with_Jiangsu/2_calc_PBS/2_calc_PBS/3_calc_PBS.pl $OUTFILE1 $OUTFILE2

echo ending at
date
