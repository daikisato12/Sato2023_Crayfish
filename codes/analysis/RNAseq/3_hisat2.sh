#!/bin/sh
#$ -S /bin/bash
#$ -t 1-15
#$ -cwd

echo original path: $PATH
ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=${data}/rnaseq/fastq/post_fastp
cd $FASTQDIR

FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=${SGE_TASK_ID} '{if (NR == line) print $0}'`
SAMPLE=$(echo ${FILENAME} | rev | cut -c 9- | rev | uniq)

GTF=${data}/ref/Procambarus_clarkii_v1.0.gtf
REF=${data}/ref/Procambarus_clarkii_genome_v1.0_only_with_genes.fa

WORKDIR=${data}/rnaseq/mapped_reads
cd $WORKDIR

hisat2 -p 2 -x ${REF} -1 ${FASTQDIR}/${SAMPLE}_1.fq.gz -2 ${FASTQDIR}/${SAMPLE}_2.fq.gz \
    -S ${SAMPLE}_hisat2b.sam --dta-cufflinks --no-discordant --no-mixed
samtools view -bS ${SAMPLE}_hisat2b.sam > ${SAMPLE}_hisat2b.bam
samtools sort ${SAMPLE}_hisat2b.bam -o ${SAMPLE}_hisat2b.s.bam
stringtie ${SAMPLE}_hisat2b.s.bam -e -G $GTF -o ${SAMPLE}_stringtie_non_novel.gtf
gffcompare -r $GTF ${SAMPLE}_stringtie_non_novel.gtf

echo ending at
date

