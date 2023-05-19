#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

WORKDIR=${data}/reseq/fastq/rawdata
cd $WORKDIR

# get sequence from Gutekunst et al. 2018
cat /codes/analysis/ReSeq/0_getseq/sra_ids_Gutekunst2018.txt | while read SAMPLE; do
	fasterq-dump ${SAMPLE} -O ./ -t ../tmp -e 10 -p
	pigz -p 10 ${SAMPLE}_1.fastq
	pigz -p 10 ${SAMPLE}_2.fastq
done

# get sequence from Xu et al. 2021
fasterq-dump SRR14457235 -O ./ -t ../tmp -e 8 -p 
pigz -p 10 SRR14457235_1.fastq
pigz -p 10 SRR14457235_2.fastq

echo ending at
date
