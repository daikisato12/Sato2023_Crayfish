qsub -l medium -l s_vmem=50G -l mem_req=50G j1_mpileup_8ind.sh
qsub -l medium -l s_vmem=20G -l mem_req=20G -hold_jid 14909435 j2_psmc_8ind.sh
qsub -l medium -l s_vmem=5G -l mem_req=5G -hold_jid 14909452 j3_psmc_8ind_bs.sh
qsub -l medium -l s_vmem=5G -l mem_req=5G -hold_jid 14909453 j4_psmc_8ind_combined.sh
