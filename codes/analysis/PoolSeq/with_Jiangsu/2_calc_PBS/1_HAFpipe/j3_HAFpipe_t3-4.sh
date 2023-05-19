#!/bin/sh
#$ -S /bin/bash
#$ -t 1-752 #94*8
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

#export LD_LIBRARY_PATH=/usr/lib64/:$LD_LIBRARY_PATH

FASTQDIR=${data}/poolseq/fastq/post_fastp
cd $FASTQDIR

TMP=$(( $(( $SGE_TASK_ID - 1 ))/ 94 + 1 ))
FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=$TMP '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 9- | rev | uniq)
BAM=$SAMPLE.sort.filtered.bam
echo $BAM
REF=${data}/mapped_reads/mpileup/hafpipe/GCA_020424385.1_ASM2042438v2_genomic_index_changed_upper.fa

TMP2=$(( $SGE_TASK_ID - $(( $TMP - 1 )) * 94 ))
CHR=`cat ${data}/mapped_reads/mpileup/hafpipe/contig_500kbp.txt | awk -v line=$TMP2 '{if (NR == line) print $0}'`

WORKDIR=${data}/poolseq/mapped_reads/mpileup/hafpipe/$SAMPLE
mkdir -p $WORKDIR
cd $WORKDIR
rm HAFpipe_${SAMPLE}_${CHR}.log
sh /codes/software/HAFpipe-line/HAFpipe_wrapper.sh -t 3,4 -b ${BAM} -g 15 -r ${REF} -l HAFpipe_${SAMPLE}_${CHR}.log -s ../HAFpipe_${CHR}.snpTable.simpute -c ${CHR} -e sanger

echo ending at
date
