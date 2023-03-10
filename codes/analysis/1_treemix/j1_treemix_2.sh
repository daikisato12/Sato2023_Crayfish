#!/bin/sh
#$ -S /bin/bash
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date


WORK_DIR=/home/daikisato177/Projects/Crayfish/Data_2022/Mapped_reads/mpileup/treemix
cd $WORK_DIR

FILE=pooled_7pops.noN.LDpruned

#ln -s ../$FILE.vcf.gz ./
#perl mkclust.pl
#/home/daikisato177/Projects/Crayfish/Build/treemix_scripts/vcf2treemix.sh $FILE.vcf.gz $FILE.clust

for i in {0..7}
do
	/home/daikisato177/Projects/Crayfish/Build/treemix/src/treemix -i $FILE.treemix.frq.gz -m $i -o $FILE.$i -bootstrap -k 500 -noss -root NewOrleans
done

#/home/daikisato177/Projects/Crayfish/Build/treemix/src/treemix -i $FILE.treemix.frq.gz -o $FILE -k 500 -noss

echo ending at
date
