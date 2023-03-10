#!/bin/sh
#$ -S /bin/bash
#$ -t 1-7
#$ -tc 3
#$ -cwd
#qsub -l medium -l s_vmem=20G -l mem_req=20G j1_pool_mpileup_eachpool_pi.sh

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

FASTQDIR=/home/daikisato177/Projects/Crayfish/Data/Pool/post_fastp
cd $FASTQDIR

FILENAME=`ls -1 *_1.fq.gz | xargs -n1 basename | awk -v line=$SGE_TASK_ID '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 9- | rev | uniq)

REF=/home/daikisato177/Projects/Crayfish/Data_2022/Ref_Xu2021/GCA_020424385.1_ASM2042438v1_genomic_index_changed.fa

WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Mapped_reads/mpileup
cd $WORK_DIR

samtools mpileup -C 50 -q 20 -Q 20 -f ${REF} ../${SAMPLE}.sort.filtered.bam > ${SAMPLE}.pileup
perl /home/daikisato177/Projects/Crayfish/Build/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input ${SAMPLE}.pileup --indel-window 5 --min-count 6  --output ${SAMPLE}.indels.gtf
perl /home/daikisato177/Projects/Crayfish/Build/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf ${SAMPLE}.indels.gtf --input ${SAMPLE}.pileup --output ${SAMPLE}.indelfiltered.pileup
rm ${SAMPLE}.pileup

echo ending at
date
