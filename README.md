# Genome Assembly and Annotation (short)
       Programm      Source                     Parameters                               Purpose
       **FastQC**         Perun enviroment           Default                                  Quality of raw reads
       **Trimmomatic**    Web                        Default                                  Removing adaptor sequences
       **SPAdes**                                    -k 21,33,55,77 –careful                  Assembly of preprocessed reads to draft genomes
       **custom script**  This repository            default (min len=200bp & min cov=5)      Additional trimming of the results
       **QUAST**                                     default                                  Evaluating the quality of the assemblies
       **CheckM**                                    default                                  Evaluating the completeness and quality of strains
       **Prokka**                                    default                                  Annotating of the high-quality assemblies

# Genome Assembly and Annotation (long)

**FastQC** (Perun) (v0.11.7) (Andrews, 2010) was used to assess the quality of raw reads, which were then quality trimmed to remove adaptor sequences using **Trimmomatic** (is downloaded from the web) (v0.36) (Bolger et al., 2014) at default settings. Draft genomes were de novo assembled from preprocessed reads using the **SPAdes** (v3.11.0) algorithm (Bankevich et al., 2012) (parameters, -k 21,33,55,77 –careful). The resulting scaffolds.fasta files were furtherly trimmed by **a custom script** (provided separately) at the parameters of the minimum length at 200 bp and minimum coverage at 5. Then the quality of the assemblies was controlled using the assembly statistics calculated by the tool **QUAST** (v5.0.2) (Gurevich et al., 2013). Then **CheckM** (Parks et al., 2015) was used to evaluated the completeness and quality of 753 assembled strains. 39 E. coli isolates with the strain heterogeneity lower than 50% were filtered. Furtherly, high-quality assemblies were annotated using **Prokka** (v1.12) software (Seemann, 2014) (default parameters) to identify open reading frames.

Andrews, S. 2010. FastQC: a quality control tool for high throughput sequence data. Babraham Bioinformatics, Babraham Institute, Cambridge, United Kingdom.

Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S. & Prjibelski, A. D. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology, 19, 455-477.

Bolger, A. M., Lohse, M. & Usadel, B. 2014. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30, 2114-2120.

Gurevich, A., Saveliev, V., Vyahhi, N. & Tesler, G. 2013. QUAST: quality assessment tool for genome assemblies. Bioinformatics, 29, 1072-1075.

Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P. & Tyson, G. W. 2015. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research, 25, 1043-1055.

Seemann, T. 2014. Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30, 2068-2069.
