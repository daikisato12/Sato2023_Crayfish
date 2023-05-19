#!/bin/sh
#$ -S /bin/bash
#$ -t 1-6
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=${data}/reseq/fastq/post_fastp
cd $FASTQDIR

FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=${SGE_TASK_ID} '{if (NR == line) print $0}'`
SAMPLE=$(echo ${FILENAME} | rev | cut -c 9- | rev | uniq)

MINDEPTH=8
NAME=$(eval echo 'dep${MINDEPTH}sites')

WORKDIR=${data}/reseq/fastq/mapped_reads/psmc
cd $WORKDIR

cat ${SAMPLE}.${NAME}.consensus.psmc ./bootstrap/${SAMPLE}.${NAME}.filtered.{1..100}.bs.psmc > ./bs_results/${SAMPLE}.${NAME}.filtered.combined.psmc
perl /codes/software/utils/psmc_plot.pl -X50000000 -p -g1 -R -x1000 -u4.59e-09 ./bs_results/${SAMPLE}.${NAME}.results ./bs_results/${SAMPLE}.${NAME}.filtered.combined.psmc

echo ending at
date
