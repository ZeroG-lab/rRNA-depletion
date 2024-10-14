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
#install.packages("MetBrewer")
#install.packages("paletteer")
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
library("paletteer")

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


# DIFFERENTIAL GENE EXPRESSION - NAA Samples ####################################################################

#Run DESeq2 on the dataset
dds.NAA <- DESeq(dds.NAA)

#Normalize counts
dds.NAA <- estimateSizeFactors(dds.NAA)

#Write counts to file
write.csv(STAR.counts, "counts_NAA.csv", quote = FALSE)
#Write normalized counts to file
write.csv(counts(dds.NAA, normalized = TRUE), "counts_NAA_normalized.csv", quote = FALSE)

#Extract results with applied filters for p-value and log2 foldchange threshold
#Set p-value (set the variable here, so the results function and any manual filtering later rely on the same value)
p.val <- 0.05
fc.limit <- 1

#Extract results (contrast needs 3 values: Which independent variable to use, Numerator=Treatment, Denominator=Control)
res_NAA <- results(dds.NAA, alpha = p.val, contrast = c("treatment", "NAA", "Ctrl"))

#Print summary
summary(res_NAA)

#Convert to dataframe
res_NAA_df <- data.frame(res_NAA)

#Write results to file
write.csv(res_NAA_df, "Foldchanges_NAA.csv", quote = FALSE)

#Write only significantly regulated genes to file (log2FC > 1 $ p.adj < 0.05)
genes.sig <- subset(res_NAA_df, abs(log2FoldChange) >= fc.limit & padj <= p.val)
write.csv(genes.sig, "Foldchanges_sig_NAA.csv", quote = FALSE)
#write.table(genes.sig, "Y:/Omics/RiboSeq/PolysomeEnrichment/Significantly_regulated_genes.txt", quote = TRUE, row.names = FALSE, sep = "\t" )

#Add Metadata to significantly regulated genes and write out
genes.sig.meta <- genes.sig
genes.sig.meta$ID <- rownames(genes.sig.meta)
Gene.Metadata <- read.table("Y:/Omics/RiboSeq/ClinoNAA_Experiment/Riboseq/Arabidopsis_Genes_Metadata.tsv", header = TRUE)
genes.sig.meta <- merge(genes.sig.meta, Gene.Metadata, all.x = TRUE, all.y = FALSE)

#Build plot of Auxin-induced gene families
genes.sig.auxin <- data.frame("Family" = factor(rep(c("SAUR", "Aux/IAA", "GH3", "ARF"), each = 3), levels = c("SAUR", "Aux/IAA", "GH3", "ARF")),
                              "Number" = c(
                                length(grep("SAUR", Gene.Metadata$Symbol, value = TRUE))-length(grep("SAUR", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE))-length(grep("SAUR", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE)),
                                length(grep("SAUR", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE)),
                                length(grep("SAUR", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE)),
                                length(grep("^IAA| IAA", Gene.Metadata$Symbol, value = TRUE))-length(grep("^IAA| IAA", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE))-length(grep("^IAA| IAA", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE)),
                                length(grep("^IAA| IAA", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE)),
                                length(grep("^IAA| IAA", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE)),
                                length(grep("GH3\\.", Gene.Metadata$Symbol, value = TRUE))- length(grep("GH3\\.", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE))-length(grep("GH3\\.", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE)),
                                length(grep("GH3\\.", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE)),
                                length(grep("GH3\\.", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE)),
                                length(grep("^ARF[0-9]| ARF[0-9]", Gene.Metadata$Symbol, value = TRUE))-length(grep("^ARF[0-9]| ARF[0-9]", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE))-length(grep("^ARF[0-9]| ARF[0-9]", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE)),
                                length(grep("^ARF[0-9]| ARF[0-9]", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE)),
                                length(grep("^ARF[0-9]| ARF[0-9]", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE))
                              ),
                              "Regulation" = factor(rep(c("None", "Up", "Down"), 4), levels = c("None", "Up", "Down"))
)


ggplot(genes.sig.auxin, aes(Family, Number, fill = Regulation))+
  geom_col()+
  scale_fill_manual(values = c("grey80", viridis::inferno(11)[7], viridis::inferno(11)[5]))+
  ylab("Number of Genes")+
  xlab("Gene Family")+
  theme_light(base_size = 9)+
  theme(axis.text = element_text(color = "black"), legend.position = "right")


#Save plot to file
ggsave("Y:/Omics/RiboSeq/PolysomeEnrichment/Plots/NAA/AuxGenes.pdf", width = 12, height = 12, units = "cm")


#Build Volcano plots

#Define plot colors
col_down <- viridis::inferno(11)[5]
col_up <- viridis::inferno(11)[7]
col_neutral <- "grey80"

#Select data to plot
plotdata <- res_NAA_df

#Modify dataframe for plotting by adding coloring information
plotdata$ID <- row.names(plotdata)
plotdata$color[abs(plotdata$log2FoldChange) > fc.limit] <- 0
plotdata$color[which(plotdata$padj > p.val)] <- 1
plotdata$color[which(abs(plotdata$log2FoldChange) < fc.limit)] <- 1

#Plotting volcano plot
ggplot(plotdata, aes(log2FoldChange, -log10(padj), label = ID, color = log2FoldChange))+
  geom_point_rast(size = 1)+
  geom_point_rast(color = col_neutral, alpha = plotdata$color, size = 0.5)+
  scale_x_continuous(limits = c(-7,7), breaks = seq(-6,6,2))+
  scale_y_continuous(limits = c(0,25), breaks = seq(0,25,5))+
  scale_color_viridis_c(option = "inferno", limits = c(-7,7), begin = 0.1, end = 0.9, guide = NULL)+
  geom_hline(yintercept = -log10(p.val), linetype = 2)+
  geom_vline(xintercept = c(-fc.limit,fc.limit), linetype = 2)+
  ggtitle("20 µM NAA 90 min")+
  xlab("Log2 Foldchange")+
  ylab("-Log10 adjusted p-value")+
  annotate("text", x = -5, y = 23, label = length(subset(plotdata, padj <= p.val & log2FoldChange <= -fc.limit)$ID), size = 8, color = col_down)+
  annotate("text", x = 5, y = 23, label = length(subset(plotdata, padj <= p.val & log2FoldChange >= fc.limit)$ID), size = 8, color = col_up)+
  theme_light(base_size = 9)+
  theme(axis.text = element_text(color = "black"))

#Save volcano plot to file
ggsave("Y:/Omics/RiboSeq/PolysomeEnrichment/Plots/NAA/Volcano.pdf", width = 12, height = 12, units = "cm")

#//Build MA plots

#Shrink foldchanges of results for better plotting (check lfcshrink help for other types of shrinkage)
res.shrink <- lfcShrink(dds.NAA, res = res_NAA, type="ashr")
#Plot MA data
plotMA(res.shrink, ylim = c(-8,8))
#Print out plot to file
dev.print(pdf, "Y:/Omics/RiboSeq/PolysomeEnrichment/Plots/NAA/NAA_MAplot.pdf")


#Individual analysis for each condition depending on the treatment ####################################################################################################################################
# BUILD METADATA FOR CTRL SAMPLES ______________________________________________________________________________________

#Generate sample table containing experiment information 
Columns_with_Ctrl_tmp <- grep("Ctrl", colnames(STAR.counts),value = TRUE)
sample.table_Ctrl <- data.frame(sample = Columns_with_Ctrl_tmp)
rownames(sample.table_Ctrl) <- sample.table_Ctrl$sample 
sample.table_Ctrl$treatment <- rep(c("Ctrl"), each = 8)
sample.table_Ctrl$condition <- rep(c("Polysomes", "Total_RNA"), each = 4)
sample.table_Ctrl$replicate <- gsub("^.*_Ctrl_", "", sample.table_Ctrl$sample)

#Generate subset of the STAR.counts table containing just informations about control samples
STAR.counts_Ctrl_tmp <- grep("Ctrl", colnames(STAR.counts), value = TRUE)
STAR.counts_Ctrl <- STAR.counts[, STAR.counts_Ctrl_tmp]

#optional: write sample table to file for export
#write.table(sample.table, "./sample_table.txt", quote = FALSE, row.names = FALSE)

dds.Ctrl <- DESeqDataSetFromMatrix(STAR.counts_Ctrl,
                                  colData = sample.table_Ctrl,
                                  design = ~ replicate+condition) 

# DIFFERENTIAL GENE EXPRESSION FOR CTRL SAMPLES - comparison between Polysome fractions and Total RNA ####################################################################

#Run DESeq2 on the dataset
dds.Ctrl <- DESeq(dds.Ctrl)

#Normalize counts
dds.Ctrl <- estimateSizeFactors(dds.Ctrl)

#Write counts to file
write.csv(STAR.counts_Ctrl, "counts_Ctrl.csv", quote = FALSE)
#Write normalized counts to file
write.csv(counts(dds.Ctrl, normalized = TRUE), "counts_Ctrl_normalized.csv", quote = FALSE)

#Extract results with applied filters for p-value and log2 foldchange threshold
#Set p-value (set the variable here, so the results function and any manual filtering later rely on the same value)
p.val <- 0.05
fc.limit <- 1

#Extract results (contrast needs 3 values: Which independent variable to use, Numerator=Treatment, Denominator=Control)
res_conditions_Ctrl <- results(dds.Ctrl, alpha = p.val, contrast = c("condition", "Polysomes", "Total_RNA"))

#Print summary
summary(res_conditions_Ctrl)

#Convert to dataframe
res_conditions_Ctrl_df <- data.frame(res_conditions_Ctrl)

#Write results to file
write.csv(res_conditions_Ctrl_df, "Foldchanges_Ctrl.csv", quote = FALSE)

#Write only significantly regulated genes to file (log2FC > 1 $ p.adj < 0.05)
genes.sig_Ctrl <- subset(res_conditions_Ctrl_df, abs(log2FoldChange) >= fc.limit & padj <= p.val)
write.csv(genes.sig_Ctrl, "Foldchanges_sig_Ctrl.csv", quote = FALSE)
#write.table(genes.sig, "Y:/Omics/RiboSeq/PolysomeEnrichment/Significantly_regulated_genes.txt", quote = TRUE, row.names = FALSE, sep = "\t" )

#Plot log2foldchange in a histogram
ggplot(genes.sig_Ctrl, aes(x = log2FoldChange))+
  geom_histogram(binwidth = 0.5, fill = "blue", color = "blue")+
  labs(title = "Histogram of log2FoldChanges in Ctrl samples", x = "log2FoldChange", y = "Frequency")+
  theme_minimal()

# BUILD METADATA FOR NAA SAMPLES ______________________________________________________________________________________

#Generate sample table containing experiment information 
Columns_with_NAA_tmp <- grep("NAA", colnames(STAR.counts),value = TRUE)
sample.table_NAA <- data.frame(sample = Columns_with_NAA_tmp)
rownames(sample.table_NAA) <- sample.table_NAA$sample 
sample.table_NAA$treatment <- rep(c("NAA"), each = 8)
sample.table_NAA$condition <- rep(c("Polysomes", "Total_RNA"), each = 4)
sample.table_NAA$replicate <- gsub("^.*_NAA_", "", sample.table_NAA$sample)

#Generate subset of the STAR.counts table containing just informations about control samples
STAR.counts_NAA_tmp <- grep("NAA", colnames(STAR.counts), value = TRUE)
STAR.counts_NAA <- STAR.counts[, STAR.counts_NAA_tmp]

#optional: write sample table to file for export
#write.table(sample.table, "./sample_table.txt", quote = FALSE, row.names = FALSE)

dds.NAA_samples <- DESeqDataSetFromMatrix(STAR.counts_NAA,
                                   colData = sample.table_NAA,
                                   design = ~ replicate+condition) 

# DIFFERENTIAL GENE EXPRESSION FOR NAA SAMPLES - comparison between Polysome fractions and Total RNA ####################################################################

#Run DESeq2 on the dataset
dds.NAA_samples <- DESeq(dds.NAA_samples)

#Normalize counts
dds.NAA_samples <- estimateSizeFactors(dds.NAA_samples)

#Write counts to file
write.csv(STAR.counts_NAA, "counts_NAA_samples.csv", quote = FALSE)
#Write normalized counts to file
write.csv(counts(dds.NAA_samples, normalized = TRUE), "counts_NAA_samples_normalized.csv", quote = FALSE)

#Extract results with applied filters for p-value and log2 foldchange threshold
#Set p-value (set the variable here, so the results function and any manual filtering later rely on the same value)
p.val <- 0.05
fc.limit <- 1

#Extract results (contrast needs 3 values: Which independent variable to use, Numerator=Treatment, Denominator=Control)
res_conditions_NAA <- results(dds.NAA_samples, alpha = p.val, contrast = c("condition", "Polysomes", "Total_RNA"))

#Print summary
summary(res_conditions_NAA)

#Convert to dataframe
res_conditions_NAA_df <- data.frame(res_conditions_NAA)

#Write results to file
write.csv(res_conditions_NAA_df, "Foldchanges_NAA.csv", quote = FALSE)

#Write only significantly regulated genes to file (log2FC > 1 $ p.adj < 0.05)
genes.sig_NAA <- subset(res_conditions_NAA_df, abs(log2FoldChange) >= fc.limit & padj <= p.val)
write.csv(genes.sig_NAA, "Foldchanges_sig_NAA.csv", quote = FALSE)

#write.table(genes.sig, "Y:/Omics/RiboSeq/PolysomeEnrichment/Significantly_regulated_genes.txt", quote = TRUE, row.names = FALSE, sep = "\t" )

#Plot log2foldchange in a histogram
ggplot(res_conditions_NAA_df, aes(x = log2FoldChange))+
  geom_histogram(binwidth = 0.5, fill = "blue", color = "blue")+
  labs(title = "Histogram of log2FoldChanges in NAA samples", x = "log2FoldChange", y = "Frequency")+
  theme_minimal()


# BUILD METADATA FOR Total RNA SAMPLES ______________________________________________________________________________________

#Generate sample table containing experiment information 
Columns_with_Total_RNA_tmp <- grep("Total_RNA", colnames(STAR.counts),value = TRUE)
sample.table_Total_RNA <- data.frame(sample = Columns_with_Total_RNA_tmp)
rownames(sample.table_Total_RNA) <- sample.table_Total_RNA$sample 
sample.table_Total_RNA$treatment <- rep(c("Ctrl", "NAA"), each = 4)
sample.table_Total_RNA$condition <- rep(c("Total_RNA"), each = 8)
sample.table_Total_RNA$replicate <- gsub("^.*_Ctrl_", "", sample.table_Total_RNA$sample)
sample.table_Total_RNA$replicate <- gsub("^.*_NAA_", "", sample.table_Total_RNA$replicate)

#Generate subset of the STAR.counts table containing just informations about control samples
STAR.counts_Total_RNA_tmp <- grep("Total_RNA", colnames(STAR.counts), value = TRUE)
STAR.counts_Total_RNA <- STAR.counts[, STAR.counts_Total_RNA_tmp]

#optional: write sample table to file for export
#write.table(sample.table, "./sample_table.txt", quote = FALSE, row.names = FALSE)

dds.Total_RNA <- DESeqDataSetFromMatrix(STAR.counts_Total_RNA,
                                   colData = sample.table_Total_RNA,
                                   design = ~ replicate+treatment)

# DIFFERENTIAL GENE EXPRESSION FOR Total RNA SAMPLES - comparison between Ctrl and NAA samples ####################################################################

#Run DESeq2 on the dataset
dds.Total_RNA <- DESeq(dds.Total_RNA)

#Normalize counts
dds.Total_RNA <- estimateSizeFactors(dds.Total_RNA)

#Write counts to file
write.csv(STAR.counts_Total_RNA, "counts_Total_RNA.csv", quote = FALSE)
#Write normalized counts to file
write.csv(counts(dds.Total_RNA, normalized = TRUE), "counts_Total_RNA_normalized.csv", quote = FALSE)

#Extract results with applied filters for p-value and log2 foldchange threshold
#Set p-value (set the variable here, so the results function and any manual filtering later rely on the same value)
p.val <- 0.05
fc.limit <- 1

#Extract results (contrast needs 3 values: Which independent variable to use, Numerator=Treatment, Denominator=Control)
res_Total_RNA <- results(dds.Total_RNA, alpha = p.val, contrast = c("treatment", "Ctrl", "NAA"))

#Print summary
summary(res_Total_RNA)

#Convert to dataframe
res_Total_RNA_df <- data.frame(res_Total_RNA)

#Write results to file
write.csv(res_Total_RNA_df, "Foldchanges_Total_RNA.csv", quote = FALSE)

#Write only significantly regulated genes to file (log2FC > 1 $ p.adj < 0.05)
genes.sig_Total_RNA <- subset(res_Total_RNA_df, abs(log2FoldChange) >= fc.limit & padj <= p.val)
write.csv(genes.sig_Total_RNA, "Foldchanges_sig_Total_RNA.csv", quote = FALSE)

#write.table(genes.sig, "Y:/Omics/RiboSeq/PolysomeEnrichment/Significantly_regulated_genes.txt", quote = TRUE, row.names = FALSE, sep = "\t" )

#Plot log2foldchange in a histogram
ggplot(res_Total_RNA_df, aes(x = log2FoldChange))+
  geom_histogram(binwidth = 0.5, fill = "blue", color = "blue")+
  labs(title = "Histogram of log2FoldChanges in Total_RNA samples", x = "log2FoldChange", y = "Frequency")+
  theme_minimal()

# BUILD METADATA FOR POLYSOME ENRICHED SAMPLES ______________________________________________________________________________________

#Generate sample table containing experiment information 
Columns_with_Polysome_Fractions_tmp <- grep("Polysome_Fractions", colnames(STAR.counts),value = TRUE)
sample.table_Polysome_Fractions <- data.frame(sample = Columns_with_Polysome_Fractions_tmp)
rownames(sample.table_Polysome_Fractions) <- sample.table_Polysome_Fractions$sample 
sample.table_Polysome_Fractions$treatment <- rep(c("Ctrl", "NAA"), each = 4)
sample.table_Polysome_Fractions$condition <- rep(c("Polysome_Fractions"), each = 8)
sample.table_Polysome_Fractions$replicate <- gsub("^.*_Ctrl_", "", sample.table_Polysome_Fractions$sample)
sample.table_Polysome_Fractions$replicate <- gsub("^.*_NAA_", "", sample.table_Polysome_Fractions$replicate)

#Generate subset of the STAR.counts table containing just informations about control samples
STAR.counts_Polysome_Fractions_tmp <- grep("Polysome_Fractions", colnames(STAR.counts), value = TRUE)
STAR.counts_Polysome_Fractions <- STAR.counts[, STAR.counts_Polysome_Fractions_tmp]

#optional: write sample table to file for export
#write.table(sample.table, "./sample_table.txt", quote = FALSE, row.names = FALSE)

dds.Polysome_Fractions <- DESeqDataSetFromMatrix(STAR.counts_Polysome_Fractions,
                                        colData = sample.table_Polysome_Fractions,
                                        design = ~ replicate+treatment)

# DIFFERENTIAL GENE EXPRESSION FOR POLYSOME ENRICHED SAMPLES - comparison between Ctrl and NAA samples ####################################################################

#Run DESeq2 on the dataset
dds.Polysome_Fractions <- DESeq(dds.Polysome_Fractions)

#Normalize counts
dds.Polysome_Fractions <- estimateSizeFactors(dds.Polysome_Fractions)

#Write counts to file
write.csv(STAR.counts_Polysome_Fractions, "counts_Polysome_Fractions.csv", quote = FALSE)
#Write normalized counts to file
write.csv(counts(dds.Polysome_Fractions, normalized = TRUE), "counts_Polysome_Fractions_normalized.csv", quote = FALSE)

#Extract results with applied filters for p-value and log2 foldchange threshold
#Set p-value (set the variable here, so the results function and any manual filtering later rely on the same value)
p.val <- 0.05
fc.limit <- 1

#Extract results (contrast needs 3 values: Which independent variable to use, Numerator=Treatment, Denominator=Control)
res_Polysome_Fractions <- results(dds.Polysome_Fractions, alpha = p.val, contrast = c("treatment", "Ctrl", "NAA"))

#Print summary
summary(res_Polysome_Fractions)

#Convert to dataframe
res_Polysome_Fractions_df <- data.frame(res_Polysome_Fractions)

#Write results to file
write.csv(res_Polysome_Fractions_df, "Foldchanges_Polysome_Fractions.csv", quote = FALSE)

#Write only significantly regulated genes to file (log2FC > 1 $ p.adj < 0.05)
genes.sig_Polysome_Fractions <- subset(res_Polysome_Fractions_df, abs(log2FoldChange) >= fc.limit & padj <= p.val)
write.csv(genes.sig_Polysome_Fractions, "Foldchanges_sig_Polysome_Fractions.csv", quote = FALSE)

#write.table(genes.sig, "Y:/Omics/RiboSeq/PolysomeEnrichment/Significantly_regulated_genes.txt", quote = TRUE, row.names = FALSE, sep = "\t" )

#Plot log2foldchange in a histogram
ggplot(res_Polysome_Fractions_df, aes(x = log2FoldChange))+
  geom_histogram(binwidth = 0.5, fill = "blue", color = "blue")+
  labs(title = "Histogram of log2FoldChanges in Polysome enriched samples", x = "log2FoldChange", y = "Frequency")+
  theme_minimal()
