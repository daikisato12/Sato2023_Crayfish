qsub -l s_vmem=10G -l mem_req=10G j1_vcf2phylip.sh
qsub -pe def_slot 10 j2_iqtree.sh