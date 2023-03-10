#!/bin/sh
#$ -S /bin/bash
#$ -t 1-3
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=${data}/Gutenkunst2018/post_fastp
cd $FASTQDIR

FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=${SGE_TASK_ID} '{if (NR == line) print $0}'`
SAMPLE=$(echo ${FILENAME} | rev | cut -c 9- | rev | uniq)

REF=${data}/${ref}/GCA_020424385.2_ASM2042438v2_genomic_index_changed.fa

WORK_DIR=${data}/Gutenkunst2018/Mapped_reads/PSMC
cd $WORK_DIR

bgzip ../${SAMPLE}.sort.filtered.vcf
mv ../${SAMPLE}.sort.filtered.vcf.gz ./
bcftools index ${SAMPLE}.sort.filtered.vcf.gz
bcftools consensus -f $REF ${SAMPLE}.sort.filtered.vcf.gz -I | gzip > ${SAMPLE}.consensus.fa.gz
${build}/psmc/utils/fq2psmcfa ${SAMPLE}.consensus.fa.gz > ${SAMPLE}.consensus.psmcfa
${build}/psmc/utils/splitfa ${SAMPLE}.consensus.psmcfa > ${SAMPLE}.split.filtered.psmcfa
${build}/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${SAMPLE}.consensus.psmc ${SAMPLE}.split.filtered.psmcfa
#${build}/psmc/utils/psmc_plot.pl -X50000000 -p -g 1 -R -x1000 -u 4.59e-09 ${SAMPLE}_results ${SAMPLE}.consensus.psmc

echo ending at
date
