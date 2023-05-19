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

FILENAME=`ls -1 *_1.fastq.gz | xargs -n1 basename | awk -v line=${SGE_TASK_ID} '{if (NR == line) print $0}'`
SAMPLE=$(echo ${FILENAME} | rev | cut -c 12- | rev | uniq)
SEED=$(eval echo '-s${SGE_TASK_ID}')

cd ../rawdata

seqtk sample ${SEED} ${SAMPLE}_1.fastq.gz 20000000 | pigz -p 10 > ../subsample/${SAMPLE}_subsample20M_1.fastq.gz
seqtk sample ${SEED} ${SAMPLE}_2.fastq.gz 20000000 | pigz -p 10 > ../subsample/${SAMPLE}_subsample20M_2.fastq.gz

done

echo ending at
date