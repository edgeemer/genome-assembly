<h2 align="center" id="Contents">Contents</h2>


<h2 align="center" id="Whole Genome Assembly Pipeline">Whole Genome Assembly Pipeline</h2>


---

<h2 align="center" id="Core Genome Assembly Pipeline">Core Genome Assembly Pipeline</h2>

---

<h2 align="center" id="Plasmid Genome Assembly Pipeline">Plasmid Genome Assembly Pipeline</h2>

---


/*Old Variant*/


<h3 align="center">Current pipleine</h3>
  
|Step|Software|Brief description|
|:--:|:--:|:--:|
| 1 | FastQC | Quality check |
| 2 | FastP | Polishing raw reads (len and coverage) |
| 3 | FastQC | Quality Check |
| 4 | SPAdes | Reads -> scaffolds |
| 5 | custom script | Polishing scaffolds (len and coverage) |
| 6 | QUAST + CheckM | Quality check |
| 7 | bowtie-build | create basenames |
| 8 | Bowtie2 | Mapping raw reads |
| 9 | samtools | SAM file to BAM |
| 10 | samtools | Sort and index the BAM file |
| 11 | bedtools | Convert the assembly to BED format |
| 12 | bedtools | Calculate coverage |
| 13 | IGV | Visualize coverage |
| 14 | Racoon | Polishing assembly |
| 15 | IGV | Visualize coverage |

1. FastQC is used to perform quality control checks on the raw reads before assembly.
1. FastP is used to polish the raw reads for length and coverage.
1. FastQC is again used to perform quality control checks on the polished reads.
1. SPAdes is used for assembling the reads into scaffolds.
1. A custom script is used to polish the scaffolds for length and coverage.
1. QUAST and CheckM are used for assessing the quality of the resulting assembly.
1. Bowtie-build is used to create basenames for the reads.
1. Bowtie2 is used to map the raw reads to the assembled scaffolds.
1. Samtools is used to convert the SAM files to BAM format.
1. Samtools is again used to sort and index the BAM files.
1. Bedtools is used to convert the assembly to BED format.
1. Bedtools is again used to calculate coverage.
1. IGV is used for visualizing coverage.
1. Racoon is used to further polish the assembly.
1. IGV is again used for visualizing coverage.

|Programm|Source|Parameters|Purpose|
|:------:|:----:|:--------:|:-----:|
|**FastQC**|Perun environment|Default|Quality of raw reads|
|**MultiQC**|Perun environment|Default|Informative reports for FastQC data|
|**FastP**|Perun environment|--qualified_quality_phred 20 --length_required 50 |Trim reads to achieve better quality|
|**SPAdes**|Perun environment|-k 21,33,55,77 –careful|Assembly of preprocessed reads to draft genomes|
|**custom script**|This repository|Default (min len=200bp & min cov=5)|Additional trimming of the results|
|**QUAST**|Perun environment|Default|Evaluating the quality of the assemblies|
|**MultiQC**|Perun environment|Default|Informative reports for QUAST data|
|**CheckM**|Perun environment|Default|Evaluating the completeness and quality of strains|
|**Prokka**|Perun environment|Default|Annotating of the high-quality assemblies|

|Programm|Example of a command|
|:------:|:--------:|
|**FastQC**|for i in ../SequenceData/RawSeqs/*; do; fastqc -o ../path/to/output $i; done|
|**MultiQC**|multiqc /path/to/fastqc_or_quast_output_directory -o /path/to/multiqc_output_directory|
| **Trimmomatic** | trimmomatic PE /input_R1.fastq.gz /input_R2.fastq.gz /output_R1_trimmed.fastq.gz /output_R1_unpaired.fastq.gz /output_R2_trimmed.fastq.gz output_R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 |
|**FastP**|fastp -i input_R1.fastq.gz -I input_R2.fastq.gz -o output_R1.fastq.gz -O output_R2.fastq.gz --qualified_quality_phred 15 --length_required 50|
|**SPAdes**|spades.py -k 21,33,55,77 --careful -1 input_R1.fastq.gz -2 input_R2.fastq.gz -o path/to/output/ |
|**custom script**||
|**QUAST**| quast /path/to/input/file.fasta --output-dir /path/to/output/ --reference /path/to/reference/genome |
|**CheckM**| checkm lineage_wf /path/to/input/ /path/to/output/ -x fasta |
|**Prokka**||
