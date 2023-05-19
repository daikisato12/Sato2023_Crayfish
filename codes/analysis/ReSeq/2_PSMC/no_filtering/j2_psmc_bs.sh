#!/bin/sh
#$ -S /bin/bash
#$ -t 1-600
#$ -tc 100
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=${data}/reseq/fastq/post_fastp
cd $FASTQDIR

TMP=$(( $(( ${SGE_TASK_ID} - 1 ))/ 100 + 1 ))
FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=${TMP} '{if (NR == line) print $0}'`
SAMPLE=$(echo ${FILENAME} | rev | cut -c 9- | rev | uniq)

REF=${data}/ref/GCA_020424385.2_ASM2042438v2_genomic_index_changed.fa

WORKDIR=${data}/reseq/fastq/mapped_reads/psmc
cd $WORKDIR

mkdir -p bootstrap
id=$(( ${SGE_TASK_ID} - $TMP * 100 + 100 ))
${build}/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o ./bootstrap/${SAMPLE}.filtered.${id}.bs.psmc ${SAMPLE}.split.filtered.psmcfa

echo ending at
date
