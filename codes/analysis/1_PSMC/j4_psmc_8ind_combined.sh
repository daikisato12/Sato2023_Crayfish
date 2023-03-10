#!/bin/sh
#$ -S /bin/bash
#$ -t 1-8
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=${data}/${fastq_afterQC_1ind}
cd $FASTQDIR

FILENAME=`ls -1 *_R1.fastq.gz | awk -v line=$SGE_TASK_ID '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 13- | rev | uniq)

WORK_DIR=${data}/${mapped_reads_1ind}/PSMC
cd $WORK_DIR

cat ${SAMPLE}.consensus.psmc ${SAMPLE}.filtered.{1..100}.bs.psmc > ${SAMPLE}.filtered.combined.psmc
/home/daikisato177/Software/psmc/utils/psmc_plot.pl -X50000000 -p -g1 -R -x1000 -u4.59e-09 ${SAMPLE}.results ${SAMPLE}.filtered.combined.psmc

echo ending at
date
