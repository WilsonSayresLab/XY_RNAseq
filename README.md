# XY_RNAseq
Alignment and filtering effects on RNAseq analysis on the X and Y chromosome

### RNAseq work flow

### Differential expression work flow 

# SAMPLES
## Download Data
The Genotyping-Tissue Expression (GTEx) Project was initiated to give researchers a resource to analyze RNAseq data among human individuals across multiple tissues (GTEx Consortium 2015, 2013). GTEx Project includes 544 recently deceased donors over 53 tissues types for 8,555 total samples, with multiple tissues collected per individual. RNA was performed using a non-strand-specific with a poly-A selection using Illumina TrueSeq and resulted in an average of 50 million 76 base pairs (bp) paired-end reads per sample. 

In sra (sequence read archive, known as short-read archive) format. Include GEO accession number 
`wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRR1850937.sra`

## Convert sra to fastq
Files will need to be converted from sra to fastq for downstream analysis. FASTQ format is a text-based format for storing both a biological sequence (usually nucleotide sequence) and its corresponding quality scores. Both the sequence letter and quality score are each encoded with a single ASCII character for brevity.

`fastq-dump sampleID.sra`
fastq-dump -> program part of the SRA toolkit that converts SRA files to fastq format
sampleID.sra -> path and name of sample.sra that you would like to convert to fastq format 

## Create fastqc reports
Fastqc reads raw sequence data from high throughput sequencers and runs a set of quality checks to produce a report. Best reports are those whose "per base sequence quality" are included in the green area of the graph & kmer content is good or average.

`fastqc sampleID.fastq`
fastqc -> Babraham bioinformatics program that that checks for quality of reads 
sampleID.fastq -> path and name of sampleID in fastq format, may also be in fastq.gz format

Move fastqc reports to desktop to visualize them as you can't open html in a terminal. Open new terminal as this will not work if logged into a HPC (high performance computing) cluster

`scp user@saguaro.a2c2.asu.edu:/Project/fastqc/sampleID_fastqc.html /Users/Desktop/`
scp	-> secure copy linux command                  
/Project/fastqc/sampleID_fastqc.html ->	path to where the files are located
/Users/Desktop -> path to where you would like to copy the files to 

# Trim for quality 
Since it has been found that raw untrimmed data leads to errors in read-mapping (Del Fabbro et al. 2013), we tested the effects of trimming versus no-trimming on read abundance, and 
approach that scans reads in the 5’-3’ direction and calculates the average quality of a group of 4 bases, read groups on the 3’-end whose quality scores were lower than the phred score 30 were removed.

The current trimming steps are:
ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
LEADING: Cut bases off the start of a read, if below a threshold quality
TRAILING: Cut bases off the end of a read, if below a threshold quality
CROP: Cut the read to a specified length
HEADCROP: Cut the specified number of bases from the start of the read
MINLEN: Drop the read if it is below a specified length
TOPHRED33: Convert quality scores to Phred-33
TOPHRED64: Convert quality scores to Phred-64
It works with FASTQ (using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used), either uncompressed or gzipp'ed FASTQ. Use of gzip format is determined based on the .gz extension.

The parameters selected were 
std: slidingwindow:4:30 leading10 trailing25 minlen40 phred33
int: slidingwindow:4:27 leading10 trailing25 minlen40 phred33
mod: slidingwindow:4:25 leading10 trailing25 minlen40 phred33

`java -jar /project/tools/trimmomatic-0.36.jar PE -phred33 /project/fastq/sampleID_input_1.fastq /project/fastq/sampleID_input_2.fastq /project/fastq/std_trim/sampleID_output_1_paired.fastq /project/fastq/std_trim/sampleID_output_1_unpaired.fastq /project/fastq/std_trim/sampleID_output_2_paired.fastq /project/fastq/std_trim/sampleID_output_2_unpaired.fastq ILLUMINACLIP:/project/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:25 SLIDINGWINDOW:4:30 MINLEN:40`
`java -jar /project/tools/trimmomatic-0.36.jar PE -phred33 /project/fastq/sampleID_input_1.fastq /project/fastq/sampleID_input_2.fastq /project/fastq/int_trim/sampleID_output_1_paired.fastq /project/fastq/int_trim/sampleID_output_1_unpaired.fastq /project/fastq/int_trim/sampleID_output_2_paired.fastq /project/fastq/int_trim/sampleID_output_2_unpaired.fastq ILLUMINACLIP:/project/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:25 SLIDINGWINDOW:4:27 MINLEN:40`
`java -jar /project/tools/trimmomatic-0.36.jar PE -phred33 /project/fastq/sampleID_input_1.fastq /project/fastq/sampleID_input_2.fastq /project/fastq/mod_trim/sampleID_output_1_paired.fastq /project/fastq/mod_trim/sampleID_output_1_unpaired.fastq /project/fastq/mod_trim/sampleID_output_2_paired.fastq /project/fastq/mod_trim/sampleID_output_2_unpaired.fastq ILLUMINACLIP:/project/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:25 SLIDINGWINDOW:4:25 MINLEN:40`

java -> indicates that this is a java program and will require java in order to run
-jar -> jar file to follow
trimmomatic-0.36.jar -> tool that will trim the raw fastq files
PE -> PE is for pair end reads. If single end then SE. 
-phred33 -> Using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used, either uncompressed or gzipp'ed FASTQ 
sampleID_input.fastq -> sampeID in fastq format
sampleID_output.fastq -> sampleID output file. Use a descriptive name such as sampleID_minlen50_sliding430_leading30_trailing40.fq
ILLUMINACLIP:TruSeq3-PE:2:30:10 -> Remove Illumina adapters provided in the TruSeq3-PE.fa file (provided). Initially Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches. These seeds will be extended and clipped if in the case of paired end reads a score of 30 is reached (about 50 bases), or in the case of single ended reads a score of 10, (about 17 bases).
LEADING:10 -> Cut bases off the start of a read, if below a threshold quality of 10
TRAILING:25 -> Cut bases off the end of a read, if below a threshold quality of 25
SLIDINGWINDOW:4:30 -> Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 30
MINLEN:40 -> Drop the read if it is below a specified length of 40
adapters-> add pathway to adapters directory




## Aligning to the reference genomes 
using -q to specify reads are fastq, --phred33 to indicate that input qualities are ASCII chars equal to the Phred+33 encoding which is used by the GTEx Illumina processing pipeline. HISAT2 parameters -p 8 launched 8 number of parallel search threads which increased alignment throughput by approximately a multiple of the number of threads, and finally -x followed by the basename of the index for the reference genome being either Def, Y-masked or YPARs-masked.



All post alignment processing described above was completed for the brain cortex, lung and whole blood tissues that were aligned to both the default genome and to the reference genome informed on the sex chromosome complement of the subject.

## Generating gene read counts 


## Computing differential expression 
Designed to assign mapped reads or fragments from pair-end genomic features from genes, exons, and promoters, featureCounts with the limma/voom (Law et al. 2014) differential expression pipeline is highly rated as one of the best-performing pipelines for the analyses of RNAseq data (SEQC/MAQC-III Consortium 2014) and was therefore chosen for our analysis.  

A gene-level information file associated with the rows of the counts matrix was created using the Homo_sapiens.GRCh38.89.gtf gene annotation file, which was used in the subread featureCounts to generate the gene count data, the gene-level.csv file contains unique gene ids for each row and the corresponding chromosome location of the gene. The gene order is the same in both the annotation Homo_sapiens.GRCh38.89.gtf and the DGEList gene-level.csv data object. 

The limma/voom vebayesfit is an empirical Bayes moderation that takes information from all genes to obtain a more precise estimate of gene-wise variability and is recommended for RNAseq analysis (Law et al. 2014). 



