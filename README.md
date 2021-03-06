# RNAseqData
These scripts detail how I perform cleaning and preliminary analysis on RNA sequencing data. 

Details on how to run these scripts, how to use Longleaf, more information on simple Unix commands, and the workflow of the analysis can be found in the document: RNA-seq Analysis Protocol_Helfrich_190708v1.2. This is a detailed explanation for my procedure, and should be followed closely to obtain the same results.

The FASTQC_Script should not be run in totality. This script contains a list of instructions and short commands to be run separately. DO NOT RUN AS A SCRIPT. 

The Alignent_Script2 details using Bowtie2 for alignment of short reads against the mouse reference genome. This should be run as an independent script.

The SortSAM_Script details using Samtools to convert SAM to BAM files, to sort the BAM files, and to index the BAM files. This should be run as an independent script.

The FeatureCounts scripts all do the same thing, but were modified for the maternal vs. fetal liver samples. These can be run as independent scripts.

Other scripts labeled "[Date]ModuleVersions" are lists of the versions of modules used in Longleaf on specific days of analysis.

The script entitled "FetMLiver_MDvsEtOH_DeSeq2.R" is the code used to run DESeq2 on the fetal mouse liver data. This should be run line by line, as not all lines will be needed and the code will not run cohesively. 

The R Markdown files in this folder detail the preliminary analysis of gene expression data, including 2 versions of overexpression analysis, gene set enrichment analysis, KEGG pathway annotation, and all of the associated graphs. The information on what the code does and the meaning of the results is included within the code. 
