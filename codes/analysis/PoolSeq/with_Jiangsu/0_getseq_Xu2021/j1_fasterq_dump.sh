#!/bin/sh
#$ -S /bin/bash
#$ -t 1-8
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

WORKDIR=${data}/poolseq/fastq/tmp_reseq
cd $WORKDIR

SAMPLE=`cat /codes/analysis/PoolSeq/0_getseq_Xu2021/sra_ids_Xu2021.txt | awk -v line=${SGE_TASK_ID} '{if (NR == line) print $0}'`
fasterq-dump ${SAMPLE} -O ./ -t ../tmp -e 10 -p
pigz -p 10 ${SAMPLE}_1.fastq
pigz -p 10 ${SAMPLE}_2.fastq

echo ending at
date
