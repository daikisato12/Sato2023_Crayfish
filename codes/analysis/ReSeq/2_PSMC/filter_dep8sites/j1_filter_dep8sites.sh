#!/bin/sh
#$ -S /bin/bash
#$ -t 1-3
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
FILTER=$(eval echo 'MIN\(INFO/DP\)\>\=${MINDEPTH}')
NAME=$(eval echo 'dep${MINDEPTH}sites')

WORK_DIR=${data}/reseq/fastq/mapped_reads/psmc
cd $WORK_DIR

bcftools view -i ${FILTER} ${SAMPLE}.sort.filtered.vcf.gz | bgzip > ${SAMPLE}.${NAME}.sort.filtered.vcf.gz
bcftools index -f ${SAMPLE}.${NAME}.sort.filtered.vcf.gz

echo ending at
date
