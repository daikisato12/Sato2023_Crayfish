#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


REF=/home/daikisato177/Projects/Crayfish/Data_2022/Ref_Xu2021/GCA_020424385.1_ASM2042438v1_genomic_index_changed.fa

WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Mapped_reads
cd $WORK_DIR

LIST=(`cat /home/daikisato177/Projects/Crayfish/Data/Pool/samplelist.txt`)
LIST2=""
for line in "${LIST[@]}"; do
	LIST2="$LIST2 $line.sort.filtered.bam";
done
echo $LIST2

sambamba mpileup -t 10 -o ./mpileup/pooled_7pops.pileup ${LIST2} --samtools -C50 -q20 -Q20 -f ${REF} 
perl /home/daikisato177/Projects/Crayfish/Build/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input ./mpileup/pooled_7pops.pileup --indel-window 5 --min-count 6  --output ./mpileup/pooled_7pops.indels.gtf
perl /home/daikisato177/Projects/Crayfish/Build/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf ./mpileup/pooled_7pops.indels.gtf --input ./mpileup/pooled_7pops.pileup --output ./mpileup/pooled_7pops.indelfiltered.pileup
#rm ./mpileup/pooled_7pops.pileup
/home/daikisato177/Build/jre1.8.0_241/bin/java -ea -Xmx7g -jar /home/daikisato177/Projects/Crayfish/Build/popoolation2_1201/mpileup2sync.jar --input ./mpileup/pooled_7pops.indelfiltered.pileup --output ./mpileup/pooled_7pops.indelfiltered.sync --fastq-type sanger --min-qual 20 --threads 10
rm ./mpileup/pooled_7pops.indelfiltered.pileup

echo ending at
date
