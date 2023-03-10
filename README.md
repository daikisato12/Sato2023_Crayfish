# Sato2023_Crayfish
IN this repository you find codes for analysis of molecular evolutionary, population genomic, and transcriptomic data generated in Sato et al. (2023) Genomic adaptive potential to cold environments in the invasive red swamp crayfish.

# Codes and data
All the scripts used to analyze data and make tables & figures are stored in `/codes/` directory. Codes used for CAFE analysis are available [here](https://github.com/ymat2/crayfish_cafe_analysis). Raw data are stored in `/data/raw_data/` directory as instructed below. Supplementary data described in the paper are available via [FigShare](https://figshare.com/s/22f447f7569e78fd7012), which includes Supplementary Data 1-5 and allele frequency data of 7 pools.

## Deposited raw data
`All_dAF_0.5over` : All the SNP coordinates with dAF above 0.5.

`All-Sweeps.bed` : The regions with signature of selective sweep in domestic rabbits, which can be found as `Database_S1` in the paper of Carneiro et al. 2014 _Science_.

`Deltas_vs_OryCun2.liftover.29mCons.mapped.chrUnFixed.sorted.DeltasOver0.8.bed_NOT_ANY_OTHER_EXON.bed` : The SNPs with dAF above 0.8 located in highly conserved non-coding regions.

`Dopamine_Ciliary_genes.tsv` : The list of genes relevant to dopamine or ciliary functions detected in DE analyses.

`Ens76~` or `Oryctolagus_cuniculus.OryCun2.0.76.gtf`: Reference information downloaded from Ensembl ver.76.

`Ribosomal_Proteins.tsv` : The list of genes coding ribosomal proteins.

## Normalization of count data
`run_TMM_scale_matrix.pl` : The Perl script to normalize the count data by TMM implemented in edgeR, downloaded from [trinity](https://github.com/trinityrnaseq/trinityrnaseq/blob/master/util/support_scripts/run_TMM_scale_matrix.pl). This can be run as follows.

```zsh
perl run_TMM_scale_matrix.pl --matrix Brain-R7-9-Ami-Hyp-count-s-3 > Brain-R7-9-Ami-Hyp-count-s-3_TMM_FPKM.matrix
```

This step has been already conducted and the normalized count data are available from Database S4.

## Differential expression analysis
`DEanalysis.R`: Conduct the DE analysis for each brain region with `exactTest`. Genes with CPM < 1 for three or more samples in a comparison were excluded. Differentially expressed genes (DEGs) were evaluated under three thresholds: 1) FDR < 0.01, 2) |Log2FC| > 1 and 3) SD/Mean of CPM in a group < 1 and no overlaps in CPM between (domestic and wild) groups.

`DatabaseS2.R`: Conduct the DE analysis for each brain region under different labeling set of wild and domestic.

## Gene enrichment analysis
DAVID (v.6.8): https://david.ncifcrf.gov

GREAT (v.4.0.4): http://great.stanford.edu/public/html/

## PPI network analysis
STRING (v.11.0): https://string-db.org
- Download the resultant network file as tsv.

Cytoscape (v.3.8.0): Settings of `MCODE` are as below.
- Uncheck "include loops"
- Degree Cutoff: 5
- "Fluff" and "Haircut"
- Node Density Cutoff: 0.7
- Node Score Cutoff: 0.7
- K-Core: 3
- Max. Depth: 100

## Tables and figures
You can reproduce all the results by running the given scripts corresponding to each analysis.
