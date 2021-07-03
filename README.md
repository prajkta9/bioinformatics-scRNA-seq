# Gene Expression Profile

## Transcriptomes of Electrophysiologically Recorded Dbx1-derived Respiratory Neurons of the preBötzinger Complex in Neonatal Mice

- The gene expression profile shown in the manuscript (doi: ) were computed using the codes available in this repository.
- We performed patch-seq on 17 Dbx1 preBötC inspiratory neurons.
- Using electrophysiology, these neurons were classified into three categories: Type-1, Type-2 and Unknown (classification of neurons is explained in depth the paper).
- We have 7 Type-1 neurons, 9 Type-2 neurons and 1 Unknown neuron.
- The sample name of the 7 Type-1 neurons are: Sample_4, Sample_5, Sample_6, Sample_15, Sample_16, Sample_17, Sample_18.
- The sample name of the 9 Type-2 neurons are: Sample_1, Sample_2, Sample_3 Sample_7, Sample_9, Sample_11, Sample_12, Sample_13, Sample_14.
- The sample name of the 1 Unknown neuron is Sample_8.
- The FASTQ files for the above samples are available on NCBI GEO database: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175642
- The NCBI GEO database also contains the fragment counts, fpkm values and gene length for each sample, computed using the code 01_computeCountsAndFpkm.R.
- In the manuscript, we used STAR v2.7.7a to generate BAM files for the each sample after trimming the adapter reads and overrepresented sequences.
- The annotation file used in the manuscript is availble on the following link: https://www.dropbox.com/s/zl6upetipp53t2j/Mus_musculus.gtf?dl=0

### 01_computeCountsAndFpkm.R
- This code calculates fragment counts, gene length and fpkm values for each sample.
- If you have generated BAM files, then place them in a folder. You will need the location for this folder to run this code.
- You will also need the annotations (in gif format) for the mouse reference genome.
- Place the annotation file in the same folder that contains the BAM files. The name of the annotation must be *Mus_musculus.gtf*.

### 02_computeDifferentialExpression.R
- This code finds the transcriptomic differences between Type-1 and Type-2 neurons and plots the log<sub>2</sub>(FPKM) of the differentially expressed (DE) genes in a heatmap.
- This code will need the fragment counts, FPKM value and information about the type of the sample. All files should be in .txt format. Place all the text files in the same folder.

### 03_plotGenesOfInterest.R
- This code bar plots for the log<sub>2</sub>(FPKM) of user-picked list of genes.
- This code will need the FPKM value and a text file that contains the list of genes. Make sure to press "return/enter" after the last gene in the list. Place all the text files in the same folder.
