qsub j0_mk_vcf.sh
qsub -l s_vmem=10G -l mem_req=10G j1_treemix_1.sh
qsub -l medium -l s_vmem=20G -l mem_req=20G j1_treemix_2.sh
