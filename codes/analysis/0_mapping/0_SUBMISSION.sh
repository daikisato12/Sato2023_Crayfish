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
