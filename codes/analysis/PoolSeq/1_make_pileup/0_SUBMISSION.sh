qsub -l s_vmem=20G -l mem_req=20G j1_fastp_qc.sh
qsub -l s_vmem=100G -l mem_req=100G j2_1_bwa-mem2_index.sh
qsub -l s_vmem=40G -l mem_req=40G -pe def_slot 8 j2_2_bwa-mem2_mapping.sh
qsub -l s_vmem=20G -l mem_req=20G -pe def_slot 10 j3_mpileup.sh