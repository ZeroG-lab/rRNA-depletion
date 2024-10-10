#RNA-Seq Pipeline for total RNA and polysome enriched samples of NAA treated and control plants 

# INSTALL / LOAD LIBRARIES ###################################################################################

# BiocManager::install("DESeq2")
# install.packages("devtools")
# install.packages("BiocManager")
# install.packages("ggplot2")
# install.packages("ggrastr")
# install.packages("plotly")
# install.packages("ashr")
# BiocManager::install("vsn")
# install.packages("RColorBrewer")
# BiocManager::install("pheatmap")
# BiocManager::install("sva") #Check citation of ComBat-seq
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.At.tair.db")
# install.packages("VennDiagram")
# BiocManager::install("GO.db")
# BiocManager::install("HDO.db")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
# install.packages("ggridges")
# BiocManager::install("pcaExplorer")
# BiocManager::install("Rsubread")
# BiocManager::install("BiocParallel")

#install_github(repo = "lcalviell/Ribo-seQC")

library("DESeq2")
library("BiocManager")
library("vsn")
library("RColorBrewer")
library("pheatmap")
library("sva")
library("ashr")

library("BiocParallel")
library("clusterProfiler")
library("org.At.tair.db")
library("pathview")
library("enrichplot")
library("ggridges")
library("pcaExplorer")

library("devtools")
library("RiboseQC")
library("Rsubread")

library("ggplot2")
library("ggrastr")
library("plotly")
library("VennDiagram") 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RNA-Seq #####################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# PREPARE RAW DATA ############################################################################################

#Concatenate files from different lanes into single file
#Rename files to a uniform format if not already done (eg: Experiment_Sample_Treatment_Replicate_Mate -> NAA_90min_Ctrl_1_1.fastq.gz)
#setwd("Y:/Omics/RiboSeq/PolysomeEnrichment/RawData/cleaned)
#file.rename(dir(), gsub("-MACE_D.*Seq", "",dir()))
#Only work with gzipped files for more efficient storage
#Keep all files in a single folder


# QUALITY CONTROL _____________________________________________________________________________________________

#Run FastQC on the read files
setwd("Y:/Omics/RiboSeq/PolysomeEnrichment/RawData/cleaned/") #switch to directory containing the reads
system(paste("wsl fastqc *.fastq.gz --outdir ../FastQC/")) #Run FastQC over all files and save results to specified folder

# SET DIRECTORIES _____________________________________________________________________________________________

#Set linux path to genome directory
genomeDir <- "/mnt/y/Omics/RiboSeq/PolysomeEnrichment/Annotations/STAR_Index"
#Set linux path to genome fasta file
genomeFastaFiles <- "/mnt/y/Omics/RiboSeq/PolysomeEnrichment/Annotations/TAIR10_chr_all.fas"
#Set linux path to annotation GTF file
GTFfile <- "/mnt/y/Omics/RiboSeq/PolysomeEnrichment/Annotations/Araport11_GTF_genes_transposons.20241001.gtf"

# STAR MAPPING ################################################################################################

#Generate Reference Genome using R to construct a command for the windows command line which in turn calls the linux subsystem (wsl)
system(paste("wsl ~/STAR-2.7.11b/source/STAR",
             "--runThreadN 30", #number of threads to use for computation
             "--runMode genomeGenerate", #run mode to generate index
             "--genomeDir", genomeDir, #directory in which to store the index (see above)
             "--genomeFastaFiles", genomeFastaFiles, #directory of the genome fasta (see above)
             "--sjdbGTFfile", GTFfile, #directory of the annotation (see above)
             "--sjdbOverhang 64", # max read length - 1
             "--genomeSAindexNbases 12")) #length of SA pre-indexing string, scaled to Arabidopsis genome according to manual: min(14, log2(GenomeLength)/2 - 1)

#Run STAR mapping in single-end mode 
setwd("Y:/Omics/RiboSeq/PolysomeEnrichment/RawData/cleaned/") #switch to directory containing the reads

for (k in list.files(pattern = ".fastq.gz$")) { 
  system(paste0("wsl echo \"Processing ", k, "\"; ~/STAR-2.7.11b/source/STAR ",
                "--readFilesCommand zcat ",
                "--genomeDir ", genomeDir, " ",
                "--readFilesIn ", k, " ",
                "--outFileNamePrefix ~/STAR_Output/", gsub("fastq\\.gz$", "", k), " ",
                "--runThreadN 30 ",
                "--quantMode GeneCounts"))
  print("Moving STAR Output from Linux filesystem to Windows...")
  system("wsl mv ~/STAR_Output/* /mnt/y/Omics/RiboSeq/PolysomeEnrichment/STAR_Output/")
  print("Done")
}

setwd("Y:/Omics/RiboSeq/PolysomeEnrichment") #switch back to original working directory

#MultiQC #####################################################################################################

#run MultiQC to get a summarised report
#run in Ubuntu shell because it's not working with the system("wsl ...") command,  
# cd /mnt/y/Omics/RiboSeq/PolysomeEnrichment/ #set working directory
# multiqc . -o MultiQC/ #performing MultiQC analysis with data from working directory, MultiQC report will be saved in the defined output directory 

# STAR DATA IMPORT ###########################################################################################

#Set working directory of data to be analyzed
setwd("Y:/Omics/RiboSeq/PolysomeEnrichment/")

#Set strandedness of the library (influences which counts to read from STAR count files)
# 2 = unstranded
# 3 = 1st read strand (stranded = yes)
# 4 = 2nd read strand (stranded = reverse)
strnd <- 2

#Read the counts of each samples ReadsPerGene file and cbind them
for (i in list.files("./STAR_Output/", pattern = "ReadsPerGene", full.names = TRUE)) {
  temp <- read.table(i)
  head(temp)
  temp <- temp[-c(1:4),]
  
  if(grep(i, list.files("./STAR_Output/", pattern = "ReadsPerGene", full.names = TRUE)) == 1){
    STAR.counts <- data.frame(temp[,strnd])
    row.names(STAR.counts) <- temp[,1]
  }else{
    STAR.counts <- cbind(STAR.counts, temp[,strnd])
    colnames(STAR.counts) <- gsub("\\.ReadsPerGene.*$","",list.files("./STAR_Output/", pattern = "ReadsPerGene")[1:length(colnames(STAR.counts))])
  }
  rm(temp,i)
}

# BUILD METADATA ______________________________________________________________________________________

#Generate sample table containing experiment information 
sample.table <- data.frame(sample = colnames(STAR.counts))
rownames(sample.table) <- sample.table$sample 
sample.table$treatment <- rep(c("Ctrl", "NAA"), each = 4)
sample.table$condition <- rep(c("Polysomes", "Total_RNA"), each = 8)
sample.table$replicate <- gsub("^.*_NAA_", "", sample.table$sample)
sample.table$replicate <- gsub("^.*_Ctrl_", "", sample.table$replicate)


#optional: write sample table to file for export
#write.table(sample.table, "./sample_table.txt", quote = FALSE, row.names = FALSE)

dds.NAA <- DESeqDataSetFromMatrix(STAR.counts,
                                  colData = sample.table,
                                  design = ~ replicate+condition+treatment)

# DATA EXPLORATION ##########################################################################################

#Perform Variance Stabilized Transformation on the count data to equalize variances across means

#Use vst for larger datasets, fast
dds.NAA.vst <- vst(dds.NAA)
#Use rlog for smaller datasets or large sequencing depth differences between samples, slower
dds.NAA.rlog <- rlog(dds.NAA)

#Visualize mean standard deviation to decide which transformation to use (flatter = better)
#Each plot is printed to a file for reference after drawing
#Combined Data
meanSdPlot(assay(dds.NAA), ranks = FALSE)
#dev.print(pdf, "./Plots/meanSd_raw.pdf")
meanSdPlot(assay(dds.NAA.vst), ranks = FALSE)
#dev.print(pdf, "./Plots/meanSd_vst.pdf")
meanSdPlot(assay(dds.NAA.rlog), ranks = FALSE)
#dev.print(pdf, "./Plots/meanSd_rlog.pdf")

#Use variance stabilized data to visualize sample distances in a heatmap (= how similar are the samples?)

#Calculate euclidian sample distance matrix from preferred stabilization
#NAA
euc.dist.NAA <- dist(t(assay(dds.NAA.vst)))
sampleDistMatrix.NAA <- as.matrix(euc.dist.NAA)

#Draw heatmap
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
#NAA
pheatmap(sampleDistMatrix.NAA,
         clustering_distance_rows = euc.dist.NAA,
         clustering_distance_cols = euc.dist.NAA,
         col = colors)
#Print out plot for reference
#dev.print(pdf, "./Plots/NAA/distheatmap.pdf") 

#Principal Component Analysis _____________________________________________________________________________

#modify DESeq's plotPCA function to gain access to more PCs than just the first two
plotPCA.ext <- function (object, intgroup = "treatment", ntop = 500, 
                         returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], PC5 = pca$x[, 5], PC6 = pca$x[, 6], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:6]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
                                                                                                            100), "% variance")) + coord_fixed()
} 

#Generate and export PCA data based on preferred stabilization
pcaData <- plotPCA.ext(dds.NAA.vst, intgroup = "treatment", returnData = TRUE) 


pcaData$condition <- rep(c("Polysome", "Total_RNA"), each = 8)

#Plot PCAs 1-6 in ggplot and save them
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcs <- data.frame("a" = c("PC1", "PC3", "PC5"), "b" = c("PC2", "PC4", "PC6"))
#Plot PCAs 1-6 in ggplot and save them
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcs <- data.frame("a" = c("PC1", "PC3", "PC5"), "b" = c("PC2", "PC4", "PC6"))
for (l in c(1:3)) {
  ggplot(pcaData, aes(.data[[pcs[l,1]]], .data[[pcs[l,2]]]))+
    geom_point(size = 3, aes(color = treatment, shape = condition))+
    #geom_text(aes(label = gsub("_C.*$|_T.*$","",name)), size = 2.5, nudge_y = -0.8)+
    xlab(paste0(pcs[l,1],": ", percentVar[l+(l-1)], "% variance")) +
    ylab(paste0(pcs[l,2], ": ", percentVar[l+l], "% variance")) +
    scale_color_viridis_d(option = "inferno", begin = 0.3, end = 0.8, name = NULL)+
    ggtitle("Principal Component Analysis")+
    theme_light(base_size = 9)+
    theme(axis.text = element_text(color = "black"), legend.position = "right")
  
  #Save PCA plot
  ggsave(filename = paste0("Y:/Omics/RiboSeq/PolysomeEnrichment/Plots/PCA_", l+(l-1), "+", l+l, ".pdf"), width = 12, height = 12, units = "cm")
}
