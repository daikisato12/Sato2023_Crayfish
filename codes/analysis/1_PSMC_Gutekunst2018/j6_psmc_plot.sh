#!/bin/sh
#$ -S /bin/bash
#$ -t 1-3
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=${data}/Gutekunst2018/post_fastp
cd $FASTQDIR

FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=${SGE_TASK_ID} '{if (NR == line) print $0}'`
SAMPLE=$(echo ${FILENAME} | rev | cut -c 9- | rev | uniq)

REF=${data}/${ref}/GCA_020424385.2_ASM2042438v2_genomic_index_changed.fa

WORK_DIR=${data}/Gutekunst2018/Mapped_reads/PSMC
cd $WORK_DIR

cat ${SAMPLE}.consensus.psmc ./bootstrap/${SAMPLE}.filtered.{1..100}.bs.psmc > ${SAMPLE}.filtered.combined.psmc
${build}/utils/psmc_plot.pl -X50000000 -p -g1 -R -x1000 -u4.59e-09 ${SAMPLE}_results ${SAMPLE}.filtered.combined.psmc

echo ending at
date
