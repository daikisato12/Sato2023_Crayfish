#220111 再度v1
qsub -l epyc -l s_vmem=20G -l mem_req=20G -pe def_slot 8 
j1_bwa-mem2_mapping_1ind.sh
qsub -l epyc -l s_vmem=20G -l mem_req=20G -pe def_slot 8 
j1_bwa-mem2_mapping_pool.sh
qsub -l epyc -l s_vmem=20G -l mem_req=20G -hold_jid 14503878 
j1_pool_mpileup_eachpool_pi.sh

qsub -l epyc -l s_vmem=40G -l mem_req=40G -pe def_slot 8 j2_elprep_1ind.sh 
#やはりelprepのところがうまくいっていないみたいなので必要
#qsub -l epyc -l s_vmem=100G -l mem_req=100G j2_elprep_1ind.sh

#220121
qsub -l epyc -l s_vmem=20G -l mem_req=20G -pe def_slot 8 j3_pileup_1ind.sh
qsub -l epyc -l s_vmem=40G -l mem_req=40G -pe def_slot 8 
j1_bwa-mem2_mapping_pool.sh
qsub -l epyc -l s_vmem=40G -l mem_req=40G -hold_jid 14537097 
j1_pool_mpileup_eachpool_pi.sh

#220202 再度v1
# pool
qsub -l epyc -l s_vmem=40G -l mem_req=40G -pe def_slot 8 
j1_bwa-mem2_mapping_pool.sh
qsub -l epyc -l s_vmem=20G -l mem_req=20G -pe def_slot 8 -hold_jid 
14577342 j2_mkpooled_pileup.sh
qsub -l epyc -l s_vmem=40G -l mem_req=40G -hold_jid 14577342 
j1_pool_mpileup_eachpool_pi.sh

# 1 ind
qsub -l epyc -l s_vmem=40G -l mem_req=40G -pe def_slot 8 
j1_bwa-mem2_mapping_1ind.sh

#220317 pooled_7pops.pileupができていなかったので再度
qsub -l medium -l s_vmem=20G -l mem_req=20G -pe def_slot 10 
j2_mkpooled_pileup.sh
