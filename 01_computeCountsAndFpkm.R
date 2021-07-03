# This program uses the BAM files that are obtained after aligning the FASTQ files with the mouse reference genome
# The FASTQ files are available publicly on NCBI GEO database (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175642)
# The mouse reference genome mm10 from Ensembl was used (hftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz)
# The annotation file used here: ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
# Alignment tool used: STAR v2.7.7a
# Note: Sample 10 did not contain sufficient number of reads so it is excluded from downstream analysis
# Create a folder that only contains the BAM files for all the samples and the annotation file for mm10
# This repo contains the file "00_sampleType.txt" that has information about the neuron type for each sample.
# The annotation file can be obtained using the following link: https://www.dropbox.com/s/zl6upetipp53t2j/Mus_musculus.gtf?dl=0
# IMPORTANT: Add the annotation file "Mus_musculus" to the folder containing BAM Files

# Enter the location for the BAM Files
{
  # Clear workspace environment
  rm(list = ls())
  
  # Import required libraries
  library(Rsubread) # Rsubread version 2.4.3 
  
  fileLocation = readline(prompt = "Enter the path to the folder that contains the BAM Aligned Files: ")
  fileName = vector()
  # Import BAM files
  for (i in seq(1, length(list.files(path = fileLocation, pattern = "bam")))){
    if (i < 10){
      fileName[i] = c(paste0(fileLocation, "Sample", i))
    }
    else{
      fileName[i] = c(paste0(fileLocation, "Sample", i+1))
    }
    
  }
  
  # Quantify reads by assigning mapped sequencing reads to exon using featureCounts tool from Rsubread package
  fc = featureCounts(files=paste0(fileName, ".bam"), annot.ext=paste0(fileLocation, "Mus_musculus.gtf"), 
                     isGTFAnnotationFile = T, isPairedEnd = T, GTF.featureType="exon", GTF.attrType="gene_name", "length")
  
  # Export read counts to a text file
  fragmentCounts = data.frame(fc$counts, stringsAsFactors = F)
  colnames(fragmentCounts) = sub(".bam", "", colnames(fragmentCounts))
  # Rearrange the columns such that it contains sample in the following order: 7 Type-1 samples, 1 Unknown sample and 9 Type-2 samples
  fragmentCounts = fragmentCounts[, c(4:6, 14:17, 8, 1:3, 7, 9:13)]
  write.table(fragmentCounts, file = paste0(fileLocation, "fragmentCounts.txt"), quote = F, sep = "\t", row.names = T, col.names = T)
  
  # Compute fragments per million mapped reads
  fpm = lapply(fragmentCounts, function(x) x * 10^6 / sum(x))
  
  # Compute fragments per kilobase million mapped reads
  geneLength = data.frame(fc$annotation[, c("GeneID", "Length")])
  options(scipen = 999)
  fpkm = lapply(fpm, function(y) y * 10^3 / geneLength$Length)
  fpkm = data.frame(fpkm, row.names = geneLength$GeneID)
  # Export the FPKM values
  write.table(fpkm, file = paste0(fileLocation, "fpkm.txt"), quote = F, sep = "\t", row.names = T, col.names = T)
}



