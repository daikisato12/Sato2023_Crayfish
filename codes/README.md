# codes
IN this repository you find codes for analyzing data and visualizing the results. Information for used softwares can be found at below links.
- [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
- [seqtk](https://github.com/lh3/seqtk)
- [fastp](https://github.com/OpenGene/fastp)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
- [elPrep](https://github.com/ExaScience/elprep)
- [Samtools](http://www.htslib.org)
- [Bcftools](http://samtools.github.io/bcftools/)
- [bedtools](https://github.com/arq5x/bedtools2)
- [Popoolation](https://sourceforge.net/p/popoolation/wiki/Main/)
- [Popoolation2](https://sourceforge.net/p/popoolation2/wiki/Main/)
- [PSMC](https://github.com/lh3/psmc)
- [HAF-pipe](https://github.com/petrov-lab/HAFpipe-line)
- [CAFE](https://github.com/ymat2/crayfish_cafe_analysis)
- [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip)
- [IQ-TREE](http://www.iqtree.org)
- [HISAT2](http://daehwankimlab.github.io/hisat2/)
- [StringTie](https://github.com/gpertea/stringtie/tree/master)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- [TCC](https://bioconductor.org/packages/release/bioc/html/TCC.html)

The structure of subdirectories in `${data}` described in codes are expected as below:
```zsh
${data}
├── blastp
├── poolseq
│   ├── fastq
│   │   ├── post_fastp
│   │   │   └── post_fastqc
│   │   └── rawdata
│   └── mapped_reads
│       └── mpileup
│           ├── fst
│           ├── hafpipe
│           │   ├── fst
│           │   └── pbs
│           ├── phylip
│           └── pi
├── ref
├── reseq
│   ├── fastq
│   │   ├── post_fastp
│   │   │   └── post_fastqc
│   │   └── rawdata
│   └── mapped_reads
│       └── psmc
│           ├── bootstrap
│           └── bs_results
└── rnaseq
    ├── count
    ├── fastq
    │   ├── post_fastp
    │   └── rawdata
    └── mapped_reads
```

## analysis
Codes used for population genomic and transcriptomic analysis are stored. The `prepDE.py` used for making count matrix is attached to the `StringTie` software and available [here](https://github.com/gpertea/stringtie/blob/master/prepDE.py).

## figures
R codes used for analyzing and visualizing the results.
