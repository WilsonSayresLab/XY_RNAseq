# XY_RNAseq
Alignment and filtering effects on RNAseq analysis on the X and Y chromosome

### RNAseq work flow

### Differential expression work flow 

## SAMPLES
The Genotyping-Tissue Expression (GTEx) Project was initiated to give researchers a resource to analyze RNAseq data among human individuals across multiple tissues (GTEx Consortium 2015, 2013). GTEx Project includes 544 recently deceased donors over 53 tissues types for 8,555 total samples, with multiple tissues collected per individual. RNA was performed using a non-strand-specific with a poly-A selection using Illumina TrueSeq and resulted in an average of 50 million 76 base pairs (bp) paired-end reads per sample. 

## Trim for quality 
Since it has been found that raw untrimmed data leads to errors in read-mapping (Del Fabbro et al. 2013), we tested the effects of trimming versus no-trimming on read abundance, and 
approach that scans reads in the 5’-3’ direction and calculates the average quality of a group of 4 bases, read groups on the 3’-end whose quality scores were lower than the phred score 30 were removed.

## Aligning to the reference genomes 
using -q to specify reads are fastq, --phred33 to indicate that input qualities are ASCII chars equal to the Phred+33 encoding which is used by the GTEx Illumina processing pipeline. HISAT2 parameters -p 8 launched 8 number of parallel search threads which increased alignment throughput by approximately a multiple of the number of threads, and finally -x followed by the basename of the index for the reference genome being either Def, Y-masked or YPARs-masked.



All post alignment processing described above was completed for the brain cortex, lung and whole blood tissues that were aligned to both the default genome and to the reference genome informed on the sex chromosome complement of the subject.

## Generating gene read counts 


## Computing differential expression 
Designed to assign mapped reads or fragments from pair-end genomic features from genes, exons, and promoters, featureCounts with the limma/voom (Law et al. 2014) differential expression pipeline is highly rated as one of the best-performing pipelines for the analyses of RNAseq data (SEQC/MAQC-III Consortium 2014) and was therefore chosen for our analysis.  

A gene-level information file associated with the rows of the counts matrix was created using the Homo_sapiens.GRCh38.89.gtf gene annotation file, which was used in the subread featureCounts to generate the gene count data, the gene-level.csv file contains unique gene ids for each row and the corresponding chromosome location of the gene. The gene order is the same in both the annotation Homo_sapiens.GRCh38.89.gtf and the DGEList gene-level.csv data object. 

The limma/voom vebayesfit is an empirical Bayes moderation that takes information from all genes to obtain a more precise estimate of gene-wise variability and is recommended for RNAseq analysis (Law et al. 2014). 



