#!/bin/sh
#$ -S /bin/bash
#$ -t 1-15
#$ -cwd

echo original path: $PATH
ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=${data}/rnaseq/fastq
cd $FASTQDIR

FILENAME=`ls -1 ./rawdata/*_1.fastq.gz | xargs -n1 basename | awk -v line=${SGE_TASK_ID} '{if (NR == line) print $0}'`
SAMPLE=$(echo ${FILENAME} | rev | cut -c 12- | rev | uniq)

FASTQ1=${FASTQDIR}/rawdata/${SAMPLE}_1.fastq.gz
FASTQ2=${FASTQDIR}/rawdata/${SAMPLE}_2.fastq.gz

## fastp
fastp -i ${FASTQ1} -I ${FASTQ2}\
	-o ./post_fastp/${SAMPLE}_1.fq.gz -O ./post_fastp/${SAMPLE}_2.fq.gz\
	-h ./post_fastp/report_${SAMPLE}.html -j ./post_fastp/report_${SAMPLE}.json -q 30 -u 30

echo ending at
date

