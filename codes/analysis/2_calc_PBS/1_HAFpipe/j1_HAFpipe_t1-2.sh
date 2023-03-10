#!/bin/sh
#$ -S /bin/bash
#$ -t 1-94
#$ -cwd

ulimit -s unlimited
echo running on `hostname`
echo starting at
date

CHR=`cat /codes/analysis/2_calc_PBS/1_HAFpipe/id_only.txt | awk -v line=$SGE_TASK_ID '{if (NR == line) print $0}'`

WORK_DIR=${data}/${mapped_reads_pool}/mpileup/HAFpipe
cd $WORK_DIR
${build}/HAFpipe-line/HAFpipe_wrapper.sh -t 1 -v ../pooled_7pops_HAFpipe.vcf -c ${CHR} -l HAFpipe_${CHR}.snpTable.log -s HAFpipe_${CHR}.snpTable -o ./
sed -i -e '1d' HAFpipe_${CHR}.snpTable
mv HAFpipe_${CHR}.snpTable ${CHR}.tmp2.txt
echo "${CHR},Ref,Aomori,atchafalaya2,kamakura2,NewOrleans,okinawa2,tonden2,triunfo2" > ${CHR}.tmp1.txt
cat ${CHR}.tmp1.txt ${CHR}.tmp2.txt > HAFpipe_${CHR}.snpTable
rm ${CHR}.tmp*

${build}/HAFpipe-line/HAFpipe_wrapper.sh -t 2 -l HAFpipe_${CHR}.snpTable.simpute.log -s HAFpipe_${CHR}.snpTable -i simpute

echo ending at
date
