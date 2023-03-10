#!/bin/sh
#$ -S /bin/bash
#$ -t 1-7
#$ -cwd
#qsub j0_calc_depth.sh

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

bedtools genomecov -ibam ../${SAMPLE}.sort.filtered.bam -g $REF > ../${SAMPLE}.sort.filtered.depth
grep -n 'LG' ../${SAMPLE}.sort.filtered.depth | perl -lane 'if(defined $hash{$F[1]}){$hash{$F[1]} = $hash{$F[1]}+$F[4] if($F[1] != 0)}else{$hash{$F[1]} = $F[4] if($F[1] != 0)}; END{foreach my $key (sort {$hash{$b} <=> $hash{$a}} keys %hash){print $key; last}}' > ../${SAMPLE}.sort.filtered.depth.mode

TARGET_COV=`cat ../${SAMPLE}.sort.filtered.depth.mode | awk '{print $1}'`
MAX_COV=$((2 * TARGET_COV))

echo $SAMPLE
echo "target coverage: $TARGET_COV"
echo "max coverage: $MAX_COV"

if [ $SAMPLE = "triunfo2" ]; then
  POOLSIZE=7
elif [ $SAMPLE = "NewOrleans" ]; then
  POOLSIZE=6
else
  POOLSIZE=8
fi

#perl /home/daikisato177/Projects/Crayfish/Build/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf ./mpileup/pooled_7pops.indels.gtf --input ./mpileup/pooled_7pops.pileup --output ./mpileup/pooled_7pops.indelfiltered.pileup
perl /home/daikisato177/Projects/Crayfish/Build/popoolation_1.2.2/basic-pipeline/subsample-pileup.pl --method withoutreplace --max-coverage $MAX_COV --target-coverage $TARGET_COV --input ${SAMPLE}.indelfiltered.pileup --output ${SAMPLE}.indelfiltered.subsampled.pileup

echo ending at
date
