#!/bin/sh
#$ -S /bin/bash
#$ -cwd

echo original path: $PATH
ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=${data}/rnaseq/fastq/post_fastp
cd $FASTQDIR
FILES1=($(ls -1 $FASTQDIR/*_1.fq.gz/))
FILES1_TMP="${FILES1[*]}"
FILELIST1=$(eval echo '${FILES1_TMP//${IFS:0:1}/,}')

FILES2=($(ls -1 $FASTQDIR/*_2.fq.gz/))
FILES2_TMP="${FILES2[*]}"
FILELIST2=$(eval echo '${FILES2_TMP//${IFS:0:1}/,}')

Trinity --full_cleanup --bypass_java_version_check --seqType fq --max_memory 32G --CPU 16 \
    --left ${FILELIST1} \
    --right ${FILELIST2}

echo ending at
date
