#!/bin/sh
#$ -S /bin/bash
#$ -t 1-658
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

#export LD_LIBRARY_PATH=/usr/lib64/:$LD_LIBRARY_PATH

FASTQDIR=${data}/${fastq_afterQC_pool}
cd $FASTQDIR

TMP=$(( $(( $SGE_TASK_ID - 1 ))/ 94 + 1 ))
FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=$TMP '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 9- | rev | uniq)
BAM=$SAMPLE.sort.filtered.bam
echo $BAM
REF=${data}/${mapped_reads_pool}/mpileup/HAFpipe/GCA_020424385.1_ASM2042438v2_genomic_index_changed_upper.fa

TMP2=$(( $SGE_TASK_ID - $(( $TMP - 1 )) * 94 ))
CHR=`cat /codes/analysis/2_calc_PBS/1_HAFpipe/id_only.txt | awk -v line=$TMP2 '{if (NR == line) print $0}'`

WORK_DIR=${data}/${mapped_reads_pool}/mpileup/HAFpipe/$SAMPLE
#mkdir -p $WORK_DIR
cd $WORK_DIR
rm HAFpipe_${SAMPLE}_${CHR}.log
${build}/HAFpipe-line/HAFpipe_wrapper.sh -t 3,4 -b ${BAM} -g 15 -r ${REF} -l HAFpipe_${SAMPLE}_${CHR}.log -s ../HAFpipe_${CHR}.snpTable.simpute -c ${CHR} -e sanger

echo ending at
date
