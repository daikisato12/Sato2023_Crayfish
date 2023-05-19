#!/bin/sh
#$ -S /bin/bash
#$ -t 1-6
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=${data}/reseq/fastq/post_fastp
cd $FASTQDIR

FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=${SGE_TASK_ID} '{if (NR == line) print $0}'`
SAMPLE=$(echo ${FILENAME} | rev | cut -c 9- | rev | uniq)

REF=${data}/ref/GCA_020424385.2_ASM2042438v2_genomic_index_changed.fa

WORKDIR=${data}/reseq/fastq/mapped_reads/psmc
cd $WORKDIR

bcftools index ${SAMPLE}.sort.filtered.vcf.gz
bcftools consensus -f $REF ${SAMPLE}.sort.filtered.vcf.gz -I | gzip > ${SAMPLE}.consensus.fa.gz
/codes/software/psmc/utils/fq2psmcfa ${SAMPLE}.consensus.fa.gz > ${SAMPLE}.consensus.psmcfa
/codes/software/psmc/utils/splitfa ${SAMPLE}.consensus.psmcfa > ${SAMPLE}.split.filtered.psmcfa
/codes/software/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${SAMPLE}.consensus.psmc ${SAMPLE}.split.filtered.psmcfa
perl /codes/software/psmc/utils/psmc_plot.pl -X50000000 -p -g 1 -R -x1000 -u 4.59e-09 ${SAMPLE}.results ${SAMPLE}.consensus.psmc

echo ending at
date
