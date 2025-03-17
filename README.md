# OConnorLab
Combinatorial selection in environmental hosts drives the evolution of a human pathogen

Jason M. Park1, Soma Ghosh1, Tamara J. O’Connor1,*

1Department of Biological Chemistry, The Johns Hopkins University School of Medicine, Baltimore, Maryland, USA

*Corresponding author

This README file provides instructions to execute the custom script for statistical analysis of the TnSeq data. The script uses fitness ratios, obtained after processing TnSeq sequencing data using analysis tools provided by the TUCF through Galaxy1 that are published as the Tn-seq data analysis tool MAGenTA2. Through Galaxy or MAGenTA, the raw reads are processed and a fitness ratio (output reads/input reads) is reported for each insertion site. The fitness tables obtained include gene annotations and fitness ratios for each chromosomal transposon insertion site: these are used by the TnseqSA script for subsequent analyses, based on an established set of statistical analyses developed by Chris Sassetti, Dana Boyd and Eric Rubin3-5. Details of the Galaxy output file can be found in McCoy KM, et. al., 2017. Please note that the TnseqSA script is generic and can be used to analyze fitness ratios obtained from any analysis software as long as the gene annotations and fitness ratios are generated for each insertion site and the formatting is consistent. The fitness tables are available in the folder designated ‘DataSets’ and can be used as example files.

Software Requirements:  R statistical package

Preparation of input file
•	Create a folder and download script, ‘TnSeqSA.r’ 
•	Create a sub-folder named ‘fitnesstables’. The code reads this name and thus has to be exact.
•	Fitness tables for each replicate sample should be placed in the sub-folder ‘fitnesstables’. (The fitness table files can be renamed as long as the code is modified to read the new file name. Files names that cannot contain special characters and spaces.)

Execution of the script
•	Before running the script for the first time, please ensure that all R packages are installed. The list of R packages required and the syntax is provided in the code from line # 7-15  
•	Open R console and execute the code. Make sure that the code is present in the current folder and the fitness files are in the sub-folder. The output files will be saved in the current folder.

Description of the Output files

1.	pool_distribution_b4_norm.jpg
Individual frequency distribution plots of the log transformed ratios for each replicate. 

2.	pool_distribution_after_norm.jpg
Frequency distribution plots of the log transformed ratios for each replicate after zero mean normalization.

3.	per_gene_stats_outlier_removed.csv
The table presents the fitness ratios, MAD and Zi scores for each insertion after outliers are removed. The following columns are present
gene	gene annotation 
insertions	total number of insertions for the gene
ratio	fitness ratio before normalization
new_ratio	fitness ratio after normalization
median	median (new_ratio) of all insertions/gene
dev	deviation from median (abs(median – new_ratio))
MAD	mean of absolute deviation
Ziscore	0.6745*(dev/MAD) for each insertion

4.	log transformed average of the fitness ratios with outliers removed.jpg
Overlay of the log transformed average of the fitness ratios lacking outliers and Gaussian distribution based on the population mean and all data points above the mean and their corresponding negative log transformed values to define where the data set deviates from a normal distribution.

5.	Zscores_per_gene.csv
This is the final output file that tabulates the number of insertions, fitness ratio, corresponding Z scores and the Student t test p values for each gene. The following columns are present
Gene	The gene targeted by the insertion site
Insertions	Total number of insertions for the gene
Average	Average of the fitness ratio 
Sdev	Standard deviation of the fitness ratios
Log_trans_ratio	log(average)
Zscore	Calculated Z score
pval	p value after performing Student t-test 
adpval	adjusted p value after Benjamini and Hochberg correction for multiple testing


References

1.	Afgan E, Baker D, Batut B, van den Beek M, Bouvier D. The Galaxy platform for accessible, reproducible and collaborative biomedical analyses: 2018 update. Nucleic Acids Res. 2018; 46:W537-W544.
2.	McCoy KM, Antonio ML, van Opijnen T. MAGenTA: a Galaxy implemented tool for complete Tn-seq analysis and data visualization. Bioinformatics. 2017; 33:2781-2783.
3.	Sassetti CM, Boyd DH, Rubin EJ. Comprehensive identification of conditionally essential genes in mycobacteria. Proc Natl Acad Sci USA. 2001; 98:12712-12717.
4.	Sassetti CM, Boyd DH, Rubin EJ. Genes required for mycobacterial growth defined by high density mutagenesis. Mol Microbiol. 2003; 48:77-84.
5.	Sassetti CM, Rubin EJ. Genetic requirements for mycobacterial survival during infection. Proc Natl Acad Sci USA. 2003; 100:12989-12994.


