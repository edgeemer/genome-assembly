<h1 align="center">Genome Assembly and Annotation (short)</h1>

|Programm|Source|Parameters|Purpose|
|:------:|:----:|:--------:|:-----:|
|**FastQC**|Perun environment|Default|Quality of raw reads|
|**MultiQC**|Perun environment|Default|Informative reports for FastQC data|
|**FastP**|Perun environment|--qualified_quality_phred 20 --length_required 50 |Trim reads to achieve better quality|
|**SPAdes**|Perun environment|-k 21,33,55,77 –careful|Assembly of preprocessed reads to draft genomes|
|**custom script**|This repository|Default (min len=200bp & min cov=5)|Additional trimming of the results|
|**QUAST**|Perun environment|Default|Evaluating the quality of the assemblies|
|**MultiQC**|Perun environment|Default|Informative reports for QUAST data|
|**CheckM**|NA|Default|Evaluating the completeness and quality of strains|
|**Prokka**|NA|Default|Annotating of the high-quality assemblies|

1. Use FastQC to assess the quality of raw reads.
2. Use FastP to quality trim the raw reads, removing adaptor sequences, using default settings, and removing parts with low quality and short length by parameters --qualified_quality_phred --length_required.
3. Use SPAdes to perform de novo assembly on the preprocessed reads with parameters -k 21,33,55,77 –careful.
4. Trim the resulting scaffolds.fasta files using a custom script (provided separately) with a minimum length parameter of 200 bp and a minimum coverage parameter of 5.
5. Evaluate the quality of the assemblies using QUAST to calculate assembly statistics.
6. Use CheckM to evaluate the completeness and quality of 753 assembled strains and filter out 39 E. coli isolates with strain heterogeneity lower than 50%.
7. Annotate the high-quality assemblies using Prokka with default parameters to identify open reading frames.

<h1 align="center">Genome Assembly and Annotation (long)</h1>

**FastQC** (Perun) (v0.11.7) (Andrews, 2010) was used to assess the quality of raw reads, which were then quality trimmed to remove adaptor sequences using **Trimmomatic** (is downloaded from the web) (v0.36) (Bolger et al., 2014) at default settings. Draft genomes were de novo assembled from preprocessed reads using the **SPAdes** (v3.11.0) algorithm (Bankevich et al., 2012) (parameters, -k 21,33,55,77 –careful). The resulting scaffolds.fasta files were furtherly trimmed by **a custom script** (provided separately) at the parameters of the minimum length at 200 bp and minimum coverage at 5. Then the quality of the assemblies was controlled using the assembly statistics calculated by the tool **QUAST** (v5.0.2) (Gurevich et al., 2013). Then **CheckM** (Parks et al., 2015) was used to evaluated the completeness and quality of 753 assembled strains. 39 E. coli isolates with the strain heterogeneity lower than 50% were filtered. Furtherly, high-quality assemblies were annotated using **Prokka** (v1.12) software (Seemann, 2014) (default parameters) to identify open reading frames.

Andrews, S. 2010. FastQC: a quality control tool for high throughput sequence data. Babraham Bioinformatics, Babraham Institute, Cambridge, United Kingdom.

Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S. & Prjibelski, A. D. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology, 19, 455-477.

Bolger, A. M., Lohse, M. & Usadel, B. 2014. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30, 2114-2120.

Gurevich, A., Saveliev, V., Vyahhi, N. & Tesler, G. 2013. QUAST: quality assessment tool for genome assemblies. Bioinformatics, 29, 1072-1075.

Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P. & Tyson, G. W. 2015. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research, 25, 1043-1055.

Seemann, T. 2014. Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30, 2068-2069.
