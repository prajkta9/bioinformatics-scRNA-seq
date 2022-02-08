# This program plots the log2(FPKM) for the selected list of genes
# This program need 2 text files:
# 1. FPKM value for each gene in each sample computed by the code, 01_computeCountsAndFpkm.R, named "fpkm.txt"
# 2. Text file that contains the gene names of interest in each row, named "plotGenes.txt"
# You'll need the location of the folder that contains all the above files.
# IMPORTANT: Keep all the text files in the same folder

# Clear the workspace environment
rm(list = ls())

# Import required libraries
library(matrixStats)

{
  # Enter the location for the text file containing fragment counts
  fileLocation = readline(prompt = "Enter the path to the folder that contains the fragment counts: ")
  
  # Import the log2(FPKM) for visualization of the DE genes
  fpkm = read.table(paste0(fileLocation, "fpkm.txt"), header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
  # Set genes names in lower case for matching them with the target list of genes
  rownames(fpkm) = tolower(rownames(fpkm))
  
  # Import the list of genes that need to be plotted like K channel genes, Na channel genes, etc.
  targetGenes = read.table(paste0(fileLocation, "plotGenes.txt"), header = F, stringsAsFactors = F)
  rownames(targetGenes) = tolower(targetGenes[, 1])
  
  
  # Extract the log2(FPKM) for the target list of genes
  targetGenesFpkm = fpkm[match(rownames(targetGenes), rownames(fpkm)), , drop = F]
  
  # Keep genes that are expressed in >25% of the samples, here 4 samples
  keepGenes = targetGenesFpkm[apply(targetGenesFpkm, MARGIN = 1, FUN = function(x) {length(which(x != 0)) >= 4}), , drop =F]
  
  keepGenesPlot = log2(keepGenes + 1) # Compute log2(FPKM) for visualization
  keepGenesMax = max(as.matrix(keepGenesPlot)) # Calculate maximum value to set the y-axis limit for the barplot
  
  # Compute the group data statistics, mean and standard deviation, for Type-1 and Type-2 neurons 
  keepGenesPlot$type1Mean = log2(rowMeans(as.matrix(keepGenes[, c(1:7)])) + 1)
  keepGenesPlot$type2Mean = log2(rowMeans(as.matrix(keepGenes[, c(9:17)])) + 1)
  # Compute the length of s.d. bars = mean + s.d.
  keepGenesPlot$type1SdPos = log2(rowMeans(as.matrix(keepGenes[, c(1:7)])) + rowSds(as.matrix(keepGenes[, c(1:7)])) + 1)
  keepGenesPlot$type2SdPos = log2(rowMeans(as.matrix(keepGenes[, c(9:17)]) + rowSds(as.matrix(keepGenes[, c(9:17)]))) + 1)
  
  # Create an array with color coding: Magenta for Type-1, Green for Unknown Type and Blue for Type-2
  typeColor = c("#D41876", "#0076BA", 
                "#D41876", "#D41876", "#D41876", "#D41876", "#D41876", "#D41876", "#D41876",
                "#1DB100", "#0076BA", "#0076BA", "#0076BA", "#0076BA", "#0076BA", "#0076BA", "#0076BA", "#0076BA", "#0076BA")
  
  # Prepare the data for bar plot
  keepGenesPlot = data.frame(keepGenesPlot[, c(18, 19, 1:17, 20, 21)])
  
  # Get the total number of genes for the plot
  nGenes = ncol(keepGenesPlot)
  # Set the plot screen margins to fit the graph
  par(mar=c(0,2,0,0))
  par(mfcol = c(nGenes, 1))
  geneNames = c(rownames(keepGenesPlot))
  
  xLabels = rownames(keepGenesPlot)
  widthBars = c(0.015, 0.015, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
  spaceBtwnBars = c(0.1, 0.1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  
  for (j in seq(1, nGenes)) {
    plt = keepGenesPlot[1:19, j] # We need group data ahead of the individual samples
    sdBars = keepGenesPlot[c(20, 21), j]
    b = barplot(plt, width = widthBars, 
                xaxt = "n", col = typeColor, las=2, border = NA, 
                ylim = c(0, keepGenesMax), cex.axis = 0.75,
                space = spaceBtwnBars)
    
    axis(side = 2, col = "#555555", labels = FALSE)
    if (!is.null(sdBars[1])){
      arrows(b[1], plt[1], b[1], sdBars[1], code = 3, angle = 90, length = 0.1, col = "#3D3D3D")
    }
    if (!is.null(sdBars[2])){
      arrows(b[2], plt[2], b[2], sdBars[2], code = 3, angle = 90, length = 0.1, col = "#3D3D3D")
    }
  }
  
  # Is any of the gene in the target list significantly different between Type-1 and Type-2?
  deGene = read.csv(paste0(fileLocation, "deLog2Fpkm.txt"), header = T, sep = "\t", row.names = 1)
  rownames(deGene) = tolower(rownames(deGene))
  keepSig = keepGenes[match(rownames(deGene), rownames(keepGenes)), , drop = F]
  keepSig = na.omit(keepSig)
  if (nrow(keepSig) != 0){
    paste0("Number of DE genes in the input list of genes are: ", rownames(keepSig))
  }
}











