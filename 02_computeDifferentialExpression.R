# This program finds the transcriptomic differences between Type-1 and Type-2 neurons
# This program need 3 text files:
# 1. fragment counts computed by the code, 01_computeCountsAndFpkm.R, named "fragmentCounts.txt"
# 2. Information about type of the neuron, named "sampleType.txt"
# 3. FPKM value for each gene in each sample computed by the code, 01_computeCountsAndFpkm.R, named "fpkm.txt"
# You'll need the location of the folder that contains the above three files.
# IMPORTANT: Keep all the text files in the same folder

# Clear workspace environment
rm(list = ls())

# Import the libraries
library(DESeq2) # DESeq2 Version used for the paper: 1.30.1
library(pheatmap) # pheatmap Version used for the paper: 1.0.12
library(RColorBrewer) # RColorBrewer Version used for the paper: 1.1-2


# Create a function to check if the unknown type neuron is present or not
# The Unknown neuron was not used for differential expression analysis
removeUnknown = function(checkData){
  # Remove the unknown type of neuron if present
  if ("Sample8" %in% colnames(checkData)){
    checkData = checkData[, -c(8), drop = F]
  }
  if ("Sample_8" %in% rownames(checkData)){
    checkData = checkData[-c(8), , drop = F]
  }
  return(checkData)
}

{
  # Enter the location for the text file containing fragment counts
  fileLocation = readline(prompt = "Enter the path to the folder that contains the fragment counts: ")
  
  # Import information about the electrophysiological type of samples and raw counts
  sampleInfoWithUnknown = read.table(paste0(fileLocation, "sampleType.txt"), header = F, sep = "\t", stringsAsFactors = F, row.names = 1)
  colnames(sampleInfoWithUnknown) = c("typeOfNeuron")
  sampleInfo = removeUnknown(sampleInfoWithUnknown)
  
  countsWithUnknown = read.table(paste0(fileLocation, "fragmentCounts.txt"), header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
  counts = removeUnknown(countsWithUnknown)
  
  # Import the FPKM values for visualization of the DE genes
  fpkmWithUnknown = read.table(paste0(fileLocation, "fpkm.txt"), header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
  fpkm = removeUnknown(fpkmWithUnknown)
  
  # Select genes with non-zero values across all the samples for DE-Seq analysis
  countsSelect = counts[rowSums(counts)!=0, , drop = F]
  
  # Create a design matrix 
  design = model.matrix(~0+factor(sampleInfo$typeOfNeuron))
  colnames(design) =levels(factor(sampleInfo$typeOfNeuron))
  
  # Set threshold for significance
  l2fc = 1.5
  padj = 0.01
  
  # Prepare the data for DE-Seq analysis
  ddsData = DESeqDataSetFromMatrix(countsSelect, colData = sampleInfo, design = ~typeOfNeuron)
  
  # Perform DESeq2 analysis
  dds = DESeq(ddsData, sfType="poscounts", minmu = 1e-8) 
  
  # Summarize the results from DESeq analysis
  res = results(dds)
  summary(res)
  mcols(res, use.names = T)
  
  # How many genes pass the significance threshold?
  paste0("Number of DE genes: ", sum(abs(res$log2FoldChange) > l2fc & res$padj < padj, na.rm = TRUE))
  
  # Order the results in the ascending order of adjusted p-value
  res= res[order(res$padj),]
  # Export all genes with their DE-Seq2 statistics
  write.table(res, file = paste0(fileLocation, "deSeqStatistics.txt"), sep = "\t", quote = F, col.names = T, row.names = T)
  
  # Create a plot of L2FC versus p-value for all the genes (Fig. 2A of the paper)
  # Set the color of the data as magenta, if up-regulated in Type-1, and blue, if up-regulated in Type-2
  res$color = "black"
  res$color[res$log2FoldChange > 0] = "#0076BA"
  res$color[res$log2FoldChange < 0] = "#D41876"
  # Create the plot
  plot(res$padj, res$log2FoldChange, pch = 20, col = res$color, ylim = c(-35, 35), bty = "n", las = 1, log = "x")
  abline(h = 0, lty = "solid", lwd = 1)
  # Add horizontal line at the L2FC threshold and vertical line at the adjusted p-value threshold
  abline(v = padj, lty = "dashed", lwd = 0.5)
  abline(h = 2^l2fc, lty = "dashed", lwd = 0.5)
  abline(h = -2^l2fc, lty = "dashed", lwd = 0.5)
  
  # Select significant genes for downstream analysis
  sigGenes = res[which(abs(res$log2FoldChange) > l2fc & res$padj < padj), ]
  sigGenes = sigGenes[order(sigGenes$log2FoldChange, decreasing = F), , drop = F]
  # Extract FPKM values of the DE genes
  sigFpkm = fpkm[match(rownames(sigGenes), rownames(fpkm)), , drop = F]
  # Compute log2(FPKM) for DE genes
  sigFpkm = log2(sigFpkm + 1)
  # Export log2(FPKM) for DE genes
  write.table(sigFpkm, file = paste0(fileLocation, "deLog2Fpkm.txt"), sep = "\t", quote = F, col.names = T, row.names = T)
  
  # Plot the log2(FPKM) in a heatmap
  geneSymbol = lapply(rownames(sigFpkm), function(x) bquote(italic(.(x)))) # Conventionally, gene symbols are written in italics
  # Arrange the samples such that all the Type-1 samples are followed Type-2 samples
  pheatmap(as.matrix(sigFpkm), scale="none", col = brewer.pal(9, "Oranges"),
           fontsize_col = 9, fontsize_row = 6,
           show_rownames = T, show_colnames = T,
           labels_row = as.expression(geneSymbol),
           cluster_rows = F, cluster_cols = F,
           main = expression("log"[2]*"(FPKM) for DE genes"), fontsize = 9,
           border_color = NA)
}




