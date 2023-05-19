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

WORKDIR=${data}/reseq/fastq/mapped_reads
cd $WORKDIR

samtools mpileup -C 50 -q 20 -Q 20 -d 200 -u -f $REF ${SAMPLE}.sort.filtered.bam \
    | bcftools call -m -v - | bcftools filter -g 10 -G 10 -O v \
    | bcftools filter -s FAIL -e '%QUAL<20 || INFO/DP<=3 || INFO/DP>=30 || INFO/MQ<=30 || INFO/MQB<1e-20 || INFO/RPB<0.0001 || INFO/BQB<0.0001 || (INFO/DP4[1]+INFO/DP4[2]<=2) || (INFO/DP4[3]+INFO/DP4[4]<=2)' -O v | bgzip > ./psmc/${SAMPLE}.sort.filtered.vcf.gz
bcftools index -f ./psmc/${SAMPLE}.sort.filtered.vcf.gz

echo ending at
date
