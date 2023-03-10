#!/bin/sh
#$ -S /bin/bash
#$ -t 1-800
#$ -tc 50
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=${data}/${fastq_afterQC_1ind}
cd $FASTQDIR

TMP=$(( $(( $SGE_TASK_ID - 1 ))/ 100 + 1 ))
FILENAME=`ls -1 *_R1.fastq.gz | awk -v line=$TMP '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 13- | rev | uniq)

WORK_DIR=${data}/${mapped_reads_1ind}/PSMC
cd $WORK_DIR

id=$(( $SGE_TASK_ID - $TMP * 100 + 100 ))
${build}/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o ${SAMPLE}.filtered.${id}.bs.psmc ${SAMPLE}.split.filtered.psmcfa

echo ending at
date