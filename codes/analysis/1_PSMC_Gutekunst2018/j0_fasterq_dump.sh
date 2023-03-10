#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Gutekunst2018/rawdata
cd $WORK_DIR

cat ../sra_ids.txt | while read line; do
	fasterq-dump $line -O ./ -t ../tmp -e 10 -p
	pigz -p 10 ${line}_1.fastq
	pigz -p 10 ${line}_2.fastq
done

#fasterq-dump SRR5115151 -O ./rawdata -t ./tmp -e 8 -p #WGS of Procambarus fallax: female1
#fasterq-dump SRR5115153 -O ./rawdata -t ./tmp -e 8 -p #WGS of Procambarus alleni: female1
#fasterq-dump SRR5115141 -O ./rawdata -t ./tmp -e 8 -p #WGS of Procambarus virginalis: Madagascar

echo ending at
date
