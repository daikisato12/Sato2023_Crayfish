#!/bin/sh
#$ -S /bin/bash
#$ -t 1-3
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=/home/daikisato177/Projects/Crayfish/Data_2022/Gutekunst2018/post_fastp
cd $FASTQDIR

FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=${SGE_TASK_ID} '{if (NR == line) print $0}'`
SAMPLE=$(echo ${FILENAME} | rev | cut -c 9- | rev | uniq)

REF=/home/daikisato177/Projects/Crayfish/Data_2022/Ref_Xu2021/GCA_020424385.1_ASM2042438v1_genomic_index_changed.fa

WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Gutekunst2018/Mapped_reads/PSMC
cd $WORK_DIR

cat ${SAMPLE}.consensus.psmc ./bootstrap/${SAMPLE}.filtered.{1..100}.bs.psmc > ${SAMPLE}.filtered.combined.psmc
/home/daikisato177/Software/psmc/utils/psmc_plot.pl -X50000000 -p -g1 -R -x1000 -u4.59e-09 ${SAMPLE}_results ${SAMPLE}.filtered.combined.psmc

echo ending at
date
