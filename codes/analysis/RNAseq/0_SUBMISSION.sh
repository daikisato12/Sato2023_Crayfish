# fastp
qsub 1_fastp_qc.sh

# Trinity
qsub -l medium -pe def_slot 16 -l s_vmem=32G -l mem_req=32G 2_trinity.sh

# Hisat2 - StringTie
qsub -pe def_slot 2 -l medium -l s_vmem=16G -l mem_req=16G 3_hisat2.sh

# Estimate raw count
cd ${data}/rnaseq/count
prepDE.py -i crayfish_stringtie_non_novel_list.txt -l 100
perl -pe 's/\,/\t/g' gene_count_matrix.csv | perl -lane 'print "$1$3$2" if ($_=~/(\S.+\S)(\s+\S+\s+\S+\s+\S+)(\s+\S+\s+\S+\s+\S+)/)' > gene_count_matrix_non_novel.tsv 
perl -pe 's/\,/\t/g' transcript_count_matrix.csv | perl -lane 'print "$1$3$2" if ($_=~/(\S.+\S)(\s+\S+\s+\S+\s+\S+)(\s+\S+\s+\S+\s+\S+)/)' > transcript_count_matrix_non_novel.tsv 
perl -pe 's/transcript_id/gene_id/' transcript_count_matrix_non_novel.tsv > temp
mv temp transcript_count_matrix_non_novel.tsv 

#Sendai day 0 vs. Sendai day 7
perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]" if( $_=~/gene_id/ or $F[1]+$F[2]+$F[3]+$F[4]+$F[5]+$F[6] > 0)' gene_count_matrix_non_novel.tsv > gene_stringtie_non_novel_count_sen0_sen7

#Sendai day 0 vs. Sapporo day 0
perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[7]\t$F[8]\t$F[9]" if( $_=~/gene_id/ or $F[1]+$F[2]+$F[3]+$F[7]+$F[8]+$F[9] > 0)' gene_count_matrix_non_novel.tsv > gene_stringtie_non_novel_count_sen0_ton0

#Sendai day 7 vs. Sapporo day 7
perl -lane 'print "$F[0]\t$F[4]\t$F[5]\t$F[6]\t$F[10]\t$F[11]\t$F[12]" if( $_=~/gene_id/ or $F[4]+$F[5]+$F[6]+$F[10]+$F[11]+$F[12] > 0)' gene_count_matrix_non_novel.tsv > gene_stringtie_non_novel_count_sen7_ton7

#Sendai day 7 vs. Sapporo day 30
perl -lane 'print "$F[0]\t$F[4]\t$F[5]\t$F[6]\t$F[13]\t$F[14]\t$F[15]" if( $_=~/gene_id/ or $F[4]+$F[5]+$F[6]+$F[13]+$F[14]+$F[15] > 0)' gene_count_matrix_non_novel.tsv > gene_stringtie_non_novel_count_sen7_ton30

#Sapporo day 0 vs. Sapporo day 7
perl -lane 'print "$F[0]\t$F[7]\t$F[8]\t$F[9]\t$F[10]\t$F[11]\t$F[12]" if( $_=~/gene_id/ or $F[7]+$F[8]+$F[9]+$F[10]+$F[11]+$F[12] > 0)' gene_count_matrix_non_novel.tsv > gene_stringtie_non_novel_count_ton0_ton7

#Sapporo day 0 vs. Sapporo day 30
perl -lane 'print "$F[0]\t$F[7]\t$F[8]\t$F[9]\t$F[13]\t$F[14]\t$F[15]" if( $_=~/gene_id/ or $F[7]+$F[8]+$F[9]+$F[13]+$F[14]+$F[15] > 0)' gene_count_matrix_non_novel.tsv > gene_stringtie_non_novel_count_ton0_ton30

#Sapporo day 7 vs. Sapporo day 30
perl -lane 'print "$F[0]\t$F[10]\t$F[11]\t$F[12]\t$F[13]\t$F[14]\t$F[15]" if( $_=~/gene_id/ or $F[10]+$F[11]+$F[12]+$F[13]+$F[14]+$F[15] > 0)' gene_count_matrix_non_novel.tsv > gene_stringtie_non_novel_count_ton7_ton30

#Sendai day 0 vs. Sapporo day 7
perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[10]\t$F[11]\t$F[12]" if( $_=~/gene_id/ or $F[1]+$F[2]+$F[3]+$F[10]+$F[11]+$F[12] > 0)' gene_count_matrix_non_novel.tsv > gene_stringtie_non_novel_count_sen0_ton7

#Sendai day 0 vs. Sapporo 3day 0
perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[13]\t$F[14]\t$F[15]" if( $_=~/gene_id/ or $F[1]+$F[2]+$F[3]+$F[13]+$F[14]+$F[15] > 0)' gene_count_matrix_non_novel.tsv > gene_stringtie_non_novel_count_sen0_ton30

#Sapporo day 0 v.s. Sendai day 7 
perl -lane 'print "$F[0]\t$F[7]\t$F[8]\t$F[9]\t$F[4]\t$F[5]\t$F[6]" if( $_=~/gene_id/ or $F[7]+$F[8]+$F[9]+$F[4]+$F[5]+$F[6] > 0)' gene_count_matrix_non_novel.tsv > gene_stringtie_non_novel_count_ton0_sen7

########## TCC in R ################
crayfish_TCC.R

