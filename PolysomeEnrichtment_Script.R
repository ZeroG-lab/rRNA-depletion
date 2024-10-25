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
library("MetBrewer")
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

#filter STAR.counts containing transcripts with 2 or more counts per million reads  
STAR.counts.cpm <- STAR.counts
for (i in 1:length(colnames(STAR.counts))) {
  
  STAR.counts.cpm[,i] <- STAR.counts[,i]/sum(STAR.counts[,i])*1000000
 
  print(i)
 
}
#Create STAR.counts subset that just contains transcripts with more than 2 or more counts per million reads 
STAR.counts.subset <- STAR.counts[apply(STAR.counts.cpm >= 2, 1, all),]  
  
 

# BUILD METADATA ______________________________________________________________________________________

#Generate sample table containing experiment information 
sample.table <- data.frame(sample = colnames(STAR.counts.subset))
rownames(sample.table) <- sample.table$sample 
sample.table$treatment <- rep(c("Ctrl", "NAA"), each = 4)
sample.table$condition <- rep(c("Polysomes", "Total_RNA"), each = 8)
sample.table$replicate <- gsub("^.*_NAA_", "", sample.table$sample)
sample.table$replicate <- gsub("^.*_Ctrl_", "", sample.table$replicate)

#optional: write sample table to file for export
#write.table(sample.table, "./sample_table.txt", quote = FALSE, row.names = FALSE)

dds.overall <- DESeqDataSetFromMatrix(STAR.counts.subset,
                                  colData = sample.table,
                                  design = ~ replicate+condition+treatment)

# DATA EXPLORATION ##########################################################################################

#Perform Variance Stabilized Transformation on the count data to equalize variances across means

#Use vst for larger datasets, fast
dds.overall.vst <- vst(dds.overall)
#Use rlog for smaller datasets or large sequencing depth differences between samples, slower
dds.overall.rlog <- rlog(dds.overall)

#Visualize mean standard deviation to decide which transformation to use (flatter = better)
#Each plot is printed to a file for reference after drawing
#Combined Data
meanSdPlot(assay(dds.overall), ranks = FALSE)
dev.print(pdf, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/meanSd_raw.pdf")
meanSdPlot(assay(dds.overall.vst), ranks = FALSE)
dev.print(pdf, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/meanSd_vst.pdf")
meanSdPlot(assay(dds.overall.rlog), ranks = FALSE)
dev.print(pdf, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/meanSd_rlog.pdf")

#Use variance stabilized data to visualize sample distances in a heatmap (= how similar are the samples?)

#Calculate euclidian sample distance matrix from preferred stabilization
#overall, unselected on a special condition
euc.dist.overall <- dist(t(assay(dds.overall.vst)))
sampleDistMatrix.overall <- as.matrix(euc.dist.overall)

#Draw heatmap
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
#overall, unselected on a special condition
pheatmap(sampleDistMatrix.overall,
         clustering_distance_rows = euc.dist.overall,
         clustering_distance_cols = euc.dist.overall,
         col = colors)
#Print out plot for reference
#dev.print(pdf, "./Plots/distheatmap.pdf") 

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
pcaData <- plotPCA.ext(dds.overall.vst, intgroup = "treatment", returnData = TRUE) 


pcaData$condition <- rep(c("Polysome", "Total_RNA"), each = 8)

#Select colors 
selected_colors <- MetBrewer::met.brewer("Morgenstern", n = 8)[c(1,7)]

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
    scale_color_manual(values = selected_colors)+
    ggtitle("Principal Component Analysis")+
    theme_light(base_size = 9)+
    theme(axis.text = element_text(color = "black"), legend.position = "right")
  
  #Save PCA plot
  ggsave(filename = paste0("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/PCA_", l+(l-1), "+", l+l, ".pdf"), width = 12, height = 12, units = "cm")
}

# DIFFERENTIAL GENE EXPRESSION ####################################################################

#Run DESeq2 on the dataset
dds.overall <- DESeq(dds.overall)

#Normalize counts
dds.overall <- estimateSizeFactors(dds.overall)

#Write counts to file
write.csv(STAR.counts, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/STAR_counts_overall.csv", quote = FALSE)
#Write normalized counts to file
write.csv(counts(dds.overall, normalized = TRUE), "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/STAR_counts_overall_normalized.csv", quote = FALSE)

#Extract results with applied filters for p-value and log2 foldchange threshold
#Set p-value (set the variable here, so the results function and any manual filtering later rely on the same value)
p.val <- 0.05
fc.limit <- 1

#Extract results (contrast needs 3 values: Which independent variable to use, Numerator=Treatment, Denominator=Control)
res_overall <- results(dds.overall, alpha = p.val, contrast = c("treatment", "NAA", "Ctrl"))

#Print summary
summary(res_overall)

#Convert to dataframe
res_overall_df <- data.frame(res_overall)

#Write results to file
write.csv(res_overall_df, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_overall.csv", quote = FALSE)

#Write only significantly regulated genes to file (log2FC > 1 $ p.adj < 0.05)
genes.sig <- subset(res_overall_df, abs(log2FoldChange) >= fc.limit & padj <= p.val)
write.csv(genes.sig, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_sig_overall.csv", quote = FALSE)
write.table(genes.sig, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Significantly_regulated_genes_overall.txt", quote = TRUE, row.names = FALSE, sep = "\t" )

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
  scale_fill_manual(values = c("grey80", met.brewer("Morgenstern", 8)[7], met.brewer("Morgenstern", 8)[2]))+
  ylab("Number of Genes")+
  xlab("Gene Family")+
  theme_light(base_size = 9)+
  theme(axis.text = element_text(color = "black"), legend.position = "right")
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/NAA/AuxGenes.pdf", width = 12, height = 12, units = "cm")


#Build Volcano plots

#Define plot colors
derain_colors <- met.brewer("Derain", n = 7)
col_down <- derain_colors[4]
col_up <- derain_colors[6]
col_neutral <- "grey80"

#Select data to plot
plotdata_overall <- res_overall_df

#Modify dataframe for plotting by adding coloring information
plotdata_overall$ID <- row.names(plotdata_overall)
plotdata_overall$Regulation <- "insignificant"  # Initialize all as neutral
plotdata_overall$Regulation[plotdata_overall$log2FoldChange >= fc.limit & plotdata_overall$padj <= p.val] <- "up"  
plotdata_overall$Regulation[plotdata_overall$log2FoldChange <= -fc.limit & plotdata_overall$padj <= p.val] <- "down"  


#Plotting volcano plot
ggplot(plotdata_overall, aes(log2FoldChange, -log10(padj), label = ID, color = Regulation))+
  geom_point_rast(size = 1)+
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,2))+
  scale_y_continuous(limits = c(0,25), breaks = seq(0,25,5))+
  scale_color_manual(values = c("down" = col_down, "up" = col_up, "insignificant" = col_neutral))+
  geom_hline(yintercept = -log10(p.val), linetype = 2)+
  geom_vline(xintercept = c(-fc.limit,fc.limit), linetype = 2)+
  ggtitle("20 µM NAA for 90 min")+
  xlab("Log2 Foldchange")+
  ylab("-Log10 adjusted p-value")+
  annotate("text", x = -4, y = 23, label = length(subset(plotdata_overall, padj <= p.val & log2FoldChange <= -fc.limit)$ID), size = 8, color = col_down)+
  annotate("text", x = 4, y = 23, label = length(subset(plotdata_overall, padj <= p.val & log2FoldChange >= fc.limit)$ID), size = 8, color = col_up)+
  theme_light(base_size = 9)+
  theme(axis.text = element_text(color = "black"))

#Save volcano plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/NAA/Volcano.pdf", width = 12, height = 12, units = "cm")


#//Build MA plots

#Shrink foldchanges of results for better plotting (check lfcshrink help for other types of shrinkage)
#res.shrink <- lfcShrink(dds.NAA, res = res_NAA, type="ashr")
#Plot MA data
#plotMA(res.shrink, ylim = c(-8,8))
#Print out plot to file
#dev.print(pdf, "Y:/Omics/RiboSeq/PolysomeEnrichment/Plots/NAA/NAA_MAplot.pdf")


# DIFFERENTIAL GENE EXPRESSION WITH INTERACTION FUNCTION ####################################################################

dds.interaction <- DESeqDataSetFromMatrix(STAR.counts.subset,
                                      colData = sample.table,
                                      design = ~ replicate+condition+treatment+condition:treatment)
#Define the reference level
dds.interaction$condition <- relevel(dds.interaction$condition, ref = "Total_RNA")
#Run DESeq2 on the dataset
dds.interaction <- DESeq(dds.interaction)

#Normalize counts
dds.interaction <- estimateSizeFactors(dds.interaction)

#Write counts to file
#write.csv(STAR.counts, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/STAR_counts_overall.csv", quote = FALSE)
#Write normalized counts to file
#write.csv(counts(dds.overall, normalized = TRUE), "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/STAR_counts_overall_normalized.csv", quote = FALSE)

#Extract results with applied filters for p-value and log2 foldchange threshold
#Set p-value (set the variable here, so the results function and any manual filtering later rely on the same value)
p.val <- 0.05
fc.limit <- 1

#write out possible results to choose what should be compared with eachother
resultsNames(dds.interaction)
#Extract results (contrast needs 3 values: Which independent variable to use, Numerator=Treatment, Denominator=Control)
res_interaction <- results(dds.interaction, alpha = p.val, name = "conditionPolysomes.treatmentNAA")

#Print summary
summary(res_interaction)

#Convert to dataframe
res_interaction_df <- data.frame(res_interaction)

#Write results to file
write.csv(res_interaction_df, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_interaction.csv", quote = FALSE)

#Write only significantly regulated genes to file (log2FC > 1 $ p.adj < 0.05)
genes.sig_interaction <- subset(res_interaction_df, abs(log2FoldChange) >= fc.limit & padj <= p.val)
write.csv(genes.sig_interaction, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_sig_interaction.csv", quote = FALSE)
write.table(genes.sig_interaction, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Significantly_regulated_genes_interaction.txt", quote = TRUE, row.names = FALSE, sep = "\t" )

#Add Metadata to significantly regulated genes and write out
genes.sig.meta_interaction <- genes.sig_interaction
genes.sig.meta_interaction$ID <- rownames(genes.sig.meta_interaction)
Gene.Metadata <- read.table("Y:/Omics/RiboSeq/ClinoNAA_Experiment/Riboseq/Arabidopsis_Genes_Metadata.tsv", header = TRUE)
genes.sig.meta_interaction <- merge(genes.sig.meta_interaction, Gene.Metadata, all.x = TRUE, all.y = FALSE)

#Build Volcano plots

#Define plot colors
#derain_colors <- met.brewer("Derain", n = 7)
#col_down <- derain_colors[4]
#col_up <- derain_colors[6]
#col_neutral <- "grey80"

#Select data to plot
#plotdata_interaction <- res_interaction_df

#Modify dataframe for plotting by adding coloring information
#plotdata_interaction$ID <- row.names(plotdata_interaction)
#plotdata_interaction$Regulation <- "insignificant"  # Initialize all as neutral
#plotdata_interaction$Regulation[plotdata_interaction$log2FoldChange >= fc.limit & plotdata_interaction$padj <= p.val] <- "up"  
#plotdata_interaction$Regulation[plotdata_interaction$log2FoldChange <= -fc.limit & plotdata_interaction$padj <= p.val] <- "down"  


#Plotting volcano plot
#ggplot(plotdata_interaction, aes(log2FoldChange, -log10(padj), label = ID, color = Regulation))+
  geom_point_rast(size = 1)+
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,2))+
  #scale_y_continuous(limits = c(0,25), breaks = seq(0,25,5))+
  scale_color_manual(values = c("down" = col_down, "up" = col_up, "insignificant" = col_neutral))+
  geom_hline(yintercept = -log10(p.val), linetype = 2)+
  geom_vline(xintercept = c(-fc.limit,fc.limit), linetype = 2)+
  #ggtitle("20 µM NAA for 90 min")+
  xlab("Log2 Foldchange")+
  ylab("-Log10 adjusted p-value")+
  annotate("text", x = -4, y = 23, label = length(subset(plotdata_interaction, padj <= p.val & log2FoldChange <= -fc.limit)$ID), size = 8, color = col_down)+
  annotate("text", x = 4, y = 23, label = length(subset(plotdata_interaction, padj <= p.val & log2FoldChange >= fc.limit)$ID), size = 8, color = col_up)+
  theme_light(base_size = 9)+
  theme(axis.text = element_text(color = "black"))

#Save volcano plot to file
#ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/NAA/Volcano_interaction.pdf", width = 12, height = 12, units = "cm")

  
#Individual analysis for each condition depending on the treatment ####################################################################################################################################
# BUILD METADATA FOR CTRL SAMPLES ______________________________________________________________________________________

#Generate sample table containing experiment information 
Columns_with_Ctrl_tmp <- grep("Ctrl", colnames(STAR.counts.subset),value = TRUE)
sample.table_Ctrl <- data.frame(sample = Columns_with_Ctrl_tmp)
rownames(sample.table_Ctrl) <- sample.table_Ctrl$sample 
sample.table_Ctrl$treatment <- rep(c("Ctrl"), each = 8)
sample.table_Ctrl$condition <- rep(c("Polysomes", "Total_RNA"), each = 4)
sample.table_Ctrl$replicate <- gsub("^.*_Ctrl_", "", sample.table_Ctrl$sample)

#Generate subset of the STAR.counts table containing just informations about control samples
STAR.counts_Ctrl_tmp <- grep("Ctrl", colnames(STAR.counts.subset), value = TRUE)
STAR.counts_Ctrl <- STAR.counts.subset[, STAR.counts_Ctrl_tmp]

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
write.csv(STAR.counts_Ctrl, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/counts_Ctrl.csv", quote = FALSE)
#Write normalized counts to file
write.csv(counts(dds.Ctrl, normalized = TRUE), "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/counts_Ctrl_normalized.csv", quote = FALSE)

#Extract results with applied filters for p-value and log2 foldchange threshold
#Set p-value (set the variable here, so the results function and any manual filtering later rely on the same value)
p.val <- 0.05
fc.limit <- 1

#Extract results (contrast needs 3 values: Which independent variable to use, Numerator=Treatment, Denominator=Control)
res_condition_Ctrl <- results(dds.Ctrl, alpha = p.val, contrast = c("condition", "Polysomes", "Total_RNA"))

#Print summary
summary(res_condition_Ctrl)

#Convert to dataframe
res_condition_Ctrl_df <- data.frame(res_condition_Ctrl)

#Write results to file
write.csv(res_condition_Ctrl_df, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_Ctrl.csv", quote = FALSE)

#Write only significantly regulated genes to file (log2FC > 1 $ p.adj < 0.05)
genes.sig_Ctrl <- subset(res_condition_Ctrl_df, abs(log2FoldChange) >= fc.limit & padj <= p.val)
write.csv(genes.sig_Ctrl, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_sig_Ctrl.csv", quote = FALSE)
#write.table(genes.sig_Ctrl, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Significantly_regulated_genes_Ctrl.txt", quote = TRUE, row.names = FALSE, sep = "\t" )

#Add Metadata to significantly regulated genes and write out
genes.sig.meta_Ctrl <- genes.sig_Ctrl
genes.sig.meta_Ctrl$ID <- rownames(genes.sig.meta_Ctrl)
Gene.Metadata <- read.table("Y:/Omics/RiboSeq/ClinoNAA_Experiment/Riboseq/Arabidopsis_Genes_Metadata.tsv", header = TRUE)
genes.sig.meta_Ctrl <- merge(genes.sig.meta_Ctrl, Gene.Metadata, all.x = TRUE, all.y = FALSE)

#Build Volcano plots

#Define plot colors
morgenstern_colors <- met.brewer("Morgenstern", n = 8)
col_down <- morgenstern_colors[1]
col_up <- morgenstern_colors[7]
col_neutral <- "grey80"

#Select data to plot
plotdata_Ctrl <- res_condition_Ctrl_df

#Modify dataframe for plotting by adding coloring information
plotdata_Ctrl$ID <- row.names(plotdata_Ctrl)
plotdata_Ctrl$Regulation <- "insignificant"  # Initialize all as neutral
plotdata_Ctrl$Regulation[plotdata_Ctrl$log2FoldChange >= fc.limit & plotdata_Ctrl$padj <= p.val] <- "up"  
plotdata_Ctrl$Regulation[plotdata_Ctrl$log2FoldChange <= -fc.limit & plotdata_Ctrl$padj <= p.val] <- "down"  

#Plotting volcano plot
ggplot(plotdata_Ctrl, aes(log2FoldChange, -log10(padj), label = ID, color = Regulation))+
  geom_point_rast(size = 1)+
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,2))+
  #scale_y_continuous(limits = c(0,25), breaks = seq(0,25,5))+
  scale_color_manual(values = c("down" = col_down, "up" = col_up, "insignificant" = col_neutral))+
  geom_hline(yintercept = -log10(p.val), linetype = 2)+
  geom_vline(xintercept = c(-fc.limit,fc.limit), linetype = 2)+
  ggtitle("DTG Under Control Conditions")+
  xlab("Log2 Foldchange")+
  ylab("-Log10 adjusted p-value")+
  annotate("text", x = -4, y = 55, label = length(subset(plotdata_Ctrl, padj <= p.val & log2FoldChange <= -fc.limit)$ID), size = 8, color = col_down)+
  annotate("text", x = 4, y = 55, label = length(subset(plotdata_Ctrl, padj <= p.val & log2FoldChange >= fc.limit)$ID), size = 8, color = col_up)+
  theme_light(base_size = 9)+
  theme(axis.text = element_text(color = "black"))
#Save volcano plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/Volcano_Ctrl.pdf", width = 12, height = 12, units = "cm")


# BUILD METADATA FOR NAA TREATED SAMPLES ______________________________________________________________________________________

#Generate sample table containing experiment information 
Columns_with_NAA_tmp <- grep("NAA", colnames(STAR.counts.subset),value = TRUE)
sample.table_NAA <- data.frame(sample = Columns_with_NAA_tmp)
rownames(sample.table_NAA) <- sample.table_NAA$sample 
sample.table_NAA$treatment <- rep(c("NAA"), each = 8)
sample.table_NAA$condition <- rep(c("Polysomes", "Total_RNA"), each = 4)
sample.table_NAA$replicate <- gsub("^.*_NAA_", "", sample.table_NAA$sample)

#Generate subset of the STAR.counts table containing just information about control samples
STAR.counts_NAA_tmp <- grep("NAA", colnames(STAR.counts.subset), value = TRUE)
STAR.counts_NAA <- STAR.counts.subset[, STAR.counts_NAA_tmp]

#optional: write sample table to file for export
#write.table(sample.table, "./sample_table.txt", quote = FALSE, row.names = FALSE)

dds.NAA <- DESeqDataSetFromMatrix(STAR.counts_NAA,
                                   colData = sample.table_NAA,
                                   design = ~ replicate+condition) 

# DIFFERENTIAL GENE EXPRESSION FOR NAA SAMPLES - comparison between Polysome fractions and Total RNA ####################################################################

#Run DESeq2 on the dataset
dds.NAA <- DESeq(dds.NAA)

#Normalize counts
dds.NAA <- estimateSizeFactors(dds.NAA)

#Write counts to file
write.csv(STAR.counts_NAA, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/counts_NAA.csv", quote = FALSE)
#Write normalized counts to file
write.csv(counts(dds.NAA, normalized = TRUE), "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/counts_NAA_normalized.csv", quote = FALSE)

#Extract results with applied filters for p-value and log2 foldchange threshold
#Set p-value (set the variable here, so the results function and any manual filtering later rely on the same value)
p.val <- 0.05
fc.limit <- 1

#Extract results (contrast needs 3 values: Which independent variable to use, Numerator=Treatment, Denominator=Control)
res_condition_NAA <- results(dds.NAA, alpha = p.val, contrast = c("condition", "Polysomes", "Total_RNA"))

#Print summary
summary(res_condition_NAA)

#Convert to dataframe
res_condition_NAA_df <- data.frame(res_condition_NAA)

#Write results to file
write.csv(res_condition_NAA_df, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_NAA.csv", quote = FALSE)

#Add Metadata to significantly regulated genes and write out
res_condition_NAA.meta <- res_condition_NAA
genes.sig.meta_Ctrl$ID <- rownames(genes.sig.meta_Ctrl)
Gene.Metadata <- read.table("Y:/Omics/RiboSeq/ClinoNAA_Experiment/Riboseq/Arabidopsis_Genes_Metadata.tsv", header = TRUE)
genes.sig.meta_Ctrl <- merge(genes.sig.meta_Ctrl, Gene.Metadata, all.x = TRUE, all.y = FALSE)

#Write only significantly regulated genes to file (log2FC > 1 $ p.adj < 0.05)
genes.sig_NAA <- subset(res_condition_NAA_df, abs(log2FoldChange) >= fc.limit & padj <= p.val)
write.csv(genes.sig_NAA, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_sig_NAA.csv", quote = FALSE)

#write.table(genes.sig_NAA, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Significantly_regulated_genes_NAA.txt", quote = TRUE, row.names = FALSE, sep = "\t" )

#Add Metadata to significantly regulated genes and write out
genes.sig.meta_NAA <- genes.sig_NAA
genes.sig.meta_NAA$ID <- rownames(genes.sig.meta_NAA)
Gene.Metadata <- read.table("Y:/Omics/RiboSeq/ClinoNAA_Experiment/Riboseq/Arabidopsis_Genes_Metadata.tsv", header = TRUE)
genes.sig.meta_NAA <- merge(genes.sig.meta_NAA, Gene.Metadata, all.x = TRUE, all.y = FALSE)

#Build Volcano plots

#Define plot colors
morgenstern_colors <- met.brewer("Morgenstern", n = 8)
col_down <- morgenstern_colors[1]
col_up <- morgenstern_colors[7]
col_neutral <- "grey80"

#Select data to plot
plotdata_NAA <- res_condition_NAA_df

#Modify dataframe for plotting by adding coloring information
plotdata_NAA$ID <- row.names(plotdata_NAA)
plotdata_NAA$Regulation <- "insignificant"  # Initialize all as neutral
plotdata_NAA$Regulation[plotdata_NAA$log2FoldChange >= fc.limit & plotdata_NAA$padj <= p.val] <- "up"  
plotdata_NAA$Regulation[plotdata_NAA$log2FoldChange <= -fc.limit & plotdata_NAA$padj <= p.val] <- "down"  

#Plotting volcano plot
ggplot(plotdata_NAA, aes(log2FoldChange, -log10(padj), label = ID, color = Regulation))+
  geom_point_rast(size = 1)+
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,2))+
  #scale_y_continuous(limits = c(0,25), breaks = seq(0,25,5))+
  scale_color_manual(values = c("down" = col_down, "up" = col_up, "insignificant" = col_neutral))+
  geom_hline(yintercept = -log10(p.val), linetype = 2)+
  geom_vline(xintercept = c(-fc.limit,fc.limit), linetype = 2)+
  ggtitle("DTG After NAA-Treatment")+
  xlab("Log2 Foldchange")+
  ylab("-Log10 adjusted p-value")+
  annotate("text", x = -4, y = 95, label = length(subset(plotdata_NAA, padj <= p.val & log2FoldChange <= -fc.limit)$ID), size = 8, color = col_down)+
  annotate("text", x = 4, y = 95, label = length(subset(plotdata_NAA, padj <= p.val & log2FoldChange >= fc.limit)$ID), size = 8, color = col_up)+
  theme_light(base_size = 9)+
  theme(axis.text = element_text(color = "black"))

#Save volcano plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/Volcano_NAA.pdf", width = 12, height = 12, units = "cm")

# BUILD METADATA FOR Total RNA SAMPLES ______________________________________________________________________________________

#Generate sample table containing experiment information 
Columns_with_Total_RNA_tmp <- grep("Total_RNA", colnames(STAR.counts.subset),value = TRUE)
sample.table_Total_RNA <- data.frame(sample = Columns_with_Total_RNA_tmp)
rownames(sample.table_Total_RNA) <- sample.table_Total_RNA$sample 
sample.table_Total_RNA$treatment <- rep(c("Ctrl", "NAA"), each = 4)
sample.table_Total_RNA$condition <- rep(c("Total_RNA"), each = 8)
sample.table_Total_RNA$replicate <- gsub("^.*_Ctrl_", "", sample.table_Total_RNA$sample)
sample.table_Total_RNA$replicate <- gsub("^.*_NAA_", "", sample.table_Total_RNA$replicate)

#Generate subset of the STAR.counts table containing just informations about control samples
STAR.counts_Total_RNA_tmp <- grep("Total_RNA", colnames(STAR.counts.subset), value = TRUE)
STAR.counts_Total_RNA <- STAR.counts.subset[, STAR.counts_Total_RNA_tmp]

#optional: write sample table to file for export
#write.table(sample.table, "./sample_table.txt", quote = FALSE, row.names = FALSE)

dds.Total_RNA <- DESeqDataSetFromMatrix(STAR.counts_Total_RNA,
                                   colData = sample.table_Total_RNA,
                                   design = ~ replicate+treatment)

# DIFFERENTIAL GENE EXPRESSION FOR Total RNA SAMPLES - comparison between NAA and Ctrl samples ####################################################################

#Run DESeq2 on the dataset
dds.Total_RNA <- DESeq(dds.Total_RNA)

#Normalize counts
dds.Total_RNA <- estimateSizeFactors(dds.Total_RNA)

#Write counts to file
write.csv(STAR.counts_Total_RNA, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/counts_Total_RNA.csv", quote = FALSE)
#Write normalized counts to file
write.csv(counts(dds.Total_RNA, normalized = TRUE), "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/counts_Total_RNA_normalized.csv", quote = FALSE)
Total_RNA_counts <- counts(dds.Total_RNA)
write.table(Total_RNA_counts, "Y:/Omics/RiboSeq/PolysomeEnrichment/Counts_Total_RNA_normalized.txt", quote = TRUE, row.names = TRUE, sep = "\t" )

#Extract results with applied filters for p-value and log2 foldchange threshold
#Set p-value (set the variable here, so the results function and any manual filtering later rely on the same value)
p.val <- 0.05
fc.limit <- 1

#Extract results (contrast needs 3 values: Which independent variable to use, Numerator=Treatment, Denominator=Control)
res_Total_RNA <- results(dds.Total_RNA, alpha = p.val, contrast = c("treatment", "NAA", "Ctrl"))

#Print summary
summary(res_Total_RNA)

#Convert to dataframe
res_Total_RNA_df <- data.frame(res_Total_RNA)

#Write results to file
write.csv(res_Total_RNA_df, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_Total_RNA.csv", quote = FALSE)

#Write only significantly regulated genes to file (log2FC > 1 $ p.adj < 0.05)
genes.sig_Total_RNA <- subset(res_Total_RNA_df, abs(log2FoldChange) >= fc.limit & padj <= p.val)
write.csv(genes.sig_Total_RNA, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_sig_Total_RNA.csv", quote = FALSE)

write.table(genes.sig_Total_RNA, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Significantly_regulated_genes_Total_RNA.txt", quote = TRUE, row.names = FALSE, sep = "\t" )

#Add Metadata to significantly regulated genes and write out
genes.sig.meta_Total_RNA <- genes.sig_Total_RNA
genes.sig.meta_Total_RNA$ID <- rownames(genes.sig.meta_Total_RNA)
Gene.Metadata <- read.table("Y:/Omics/RiboSeq/ClinoNAA_Experiment/Riboseq/Arabidopsis_Genes_Metadata.tsv", header = TRUE)
genes.sig.meta_Total <- merge(genes.sig.meta_Total_RNA, Gene.Metadata, all.x = TRUE, all.y = FALSE)


#Build Volcano plots

#Define plot colors
morgenstern_colors <- met.brewer("Morgenstern", n = 8)
col_down <- morgenstern_colors[1]
col_up <- morgenstern_colors[7]
col_neutral <- "grey80"

#Select data to plot
plotdata_Total_RNA <- res_Total_RNA_df

#Modify dataframe for plotting by adding coloring information
plotdata_Total_RNA$ID <- row.names(plotdata_Total_RNA)
plotdata_Total_RNA$Regulation <- "insignificant"  # Initialize all as neutral
plotdata_Total_RNA$Regulation[plotdata_Total_RNA$log2FoldChange >= fc.limit & plotdata_Total_RNA$padj <= p.val] <- "up"  
plotdata_Total_RNA$Regulation[plotdata_Total_RNA$log2FoldChange <= -fc.limit & plotdata_Total_RNA$padj <= p.val] <- "down"  

#Plotting volcano plot
ggplot(plotdata_Total_RNA, aes(log2FoldChange, -log10(padj), label = ID, color = Regulation))+
  geom_point_rast(size = 1)+
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,2))+
  #scale_y_continuous(limits = c(0,25), breaks = seq(0,25,5))+
  scale_color_manual(values = c("down" = col_down, "up" = col_up, "insignificant" = col_neutral))+
  geom_hline(yintercept = -log10(p.val), linetype = 2)+
  geom_vline(xintercept = c(-fc.limit,fc.limit), linetype = 2)+
  ggtitle("DTG in Total RNA samples under NAA vs Ctrl conditions")+
  xlab("Log2 Foldchange")+
  ylab("-Log10 adjusted p-value")+
  annotate("text", x = -4, y = 160, label = length(subset(plotdata_Total_RNA, padj <= p.val & log2FoldChange <= -fc.limit)$ID), size = 8, color = col_down)+
  annotate("text", x = 4, y = 160, label = length(subset(plotdata_Total_RNA, padj <= p.val & log2FoldChange >= fc.limit)$ID), size = 8, color = col_up)+
  theme_light(base_size = 9)+
  theme(axis.text = element_text(color = "black"))

#Save volcano plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/Volcano_Total_RNA.pdf", width = 12, height = 12, units = "cm")


# BUILD METADATA FOR POLYSOME ENRICHED SAMPLES ______________________________________________________________________________________

#Generate sample table containing experiment information 
Columns_with_Polysome_Fractions_tmp <- grep("Polysome_Fractions", colnames(STAR.counts.subset),value = TRUE)
sample.table_Polysome_Fractions <- data.frame(sample = Columns_with_Polysome_Fractions_tmp)
rownames(sample.table_Polysome_Fractions) <- sample.table_Polysome_Fractions$sample 
sample.table_Polysome_Fractions$treatment <- rep(c("Ctrl", "NAA"), each = 4)
sample.table_Polysome_Fractions$condition <- rep(c("Polysome_Fractions"), each = 8)
sample.table_Polysome_Fractions$replicate <- gsub("^.*_Ctrl_", "", sample.table_Polysome_Fractions$sample)
sample.table_Polysome_Fractions$replicate <- gsub("^.*_NAA_", "", sample.table_Polysome_Fractions$replicate)

#Generate subset of the STAR.counts table containing just informations about control samples
STAR.counts_Polysome_Fractions_tmp <- grep("Polysome_Fractions", colnames(STAR.counts.subset), value = TRUE)
STAR.counts_Polysome_Fractions <- STAR.counts.subset[, STAR.counts_Polysome_Fractions_tmp]

#optional: write sample table to file for export
#write.table(sample.table, "./sample_table.txt", quote = FALSE, row.names = FALSE)

dds.Polysome_Fractions <- DESeqDataSetFromMatrix(STAR.counts_Polysome_Fractions,
                                        colData = sample.table_Polysome_Fractions,
                                        design = ~ replicate+treatment)

# DIFFERENTIAL GENE EXPRESSION FOR POLYSOME ENRICHED SAMPLES - comparison between NAA and Ctrl samples ####################################################################

#Run DESeq2 on the dataset
dds.Polysome_Fractions <- DESeq(dds.Polysome_Fractions)

#Normalize counts
dds.Polysome_Fractions <- estimateSizeFactors(dds.Polysome_Fractions)

#Write counts to file
write.csv(STAR.counts_Polysome_Fractions, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/counts_Polysome_Fractions.csv", quote = FALSE)
#Write normalized counts to file
write.csv(counts(dds.Polysome_Fractions, normalized = TRUE), "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/counts_Polysome_Fractions_normalized.csv", quote = FALSE)

#Extract results with applied filters for p-value and log2 foldchange threshold
#Set p-value (set the variable here, so the results function and any manual filtering later rely on the same value)
p.val <- 0.05
fc.limit <- 1

#Extract results (contrast needs 3 values: Which independent variable to use, Numerator=Treatment, Denominator=Control)
res_Polysome_Fractions <- results(dds.Polysome_Fractions, alpha = p.val, contrast = c("treatment", "NAA", "Ctrl"))

#Print summary
summary(res_Polysome_Fractions)

#Convert to dataframe
res_Polysome_Fractions_df <- data.frame(res_Polysome_Fractions)

#Write results to file
write.csv(res_Polysome_Fractions_df, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_Polysome_Fractions.csv", quote = FALSE)

#Write only significantly regulated genes to file (log2FC > 1 $ p.adj < 0.05)
genes.sig_Polysome_Fractions <- subset(res_Polysome_Fractions_df, abs(log2FoldChange) >= fc.limit & padj <= p.val)
write.csv(genes.sig_Polysome_Fractions, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Foldchanges_sig_Polysome_Fractions.csv", quote = FALSE)

write.table(genes.sig, "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Counts/FoldChanges/Significantly_regulated_genes_Polysome_Fractions.txt", quote = TRUE, row.names = FALSE, sep = "\t" )

#Add Metadata to significantly regulated genes and write out
genes.sig.meta_Polysome_Fractions <- genes.sig_Polysome_Fractions
genes.sig.meta_Polysome_Fractions$ID <- rownames(genes.sig.meta_Polysome_Fractions)
Gene.Metadata <- read.table("Y:/Omics/RiboSeq/ClinoNAA_Experiment/Riboseq/Arabidopsis_Genes_Metadata.tsv", header = TRUE)
genes.sig.meta_Polysome_Fractions <- merge(genes.sig.meta_Polysome_Fractions, Gene.Metadata, all.x = TRUE, all.y = FALSE)


#Build Volcano plots

#Define plot colors
morgenstern_colors <- met.brewer("Morgenstern", n = 8)
col_down <- morgenstern_colors[1]
col_up <- morgenstern_colors[7]
col_neutral <- "grey80"


#Select data to plot
plotdata_Polysome_Fractions <- res_Polysome_Fractions_df

#Modify dataframe for plotting by adding coloring information
plotdata_Polysome_Fractions$ID <- row.names(plotdata_Polysome_Fractions)
plotdata_Polysome_Fractions$Regulation <- "insignificant"  # Initialize all as neutral
plotdata_Polysome_Fractions$Regulation[plotdata_Polysome_Fractions$log2FoldChange >= fc.limit & plotdata_Polysome_Fractions$padj <= p.val] <- "up"  
plotdata_Polysome_Fractions$Regulation[plotdata_Polysome_Fractions$log2FoldChange <= -fc.limit & plotdata_Polysome_Fractions$padj <= p.val] <- "down"  

#Plotting volcano plot
ggplot(plotdata_Polysome_Fractions, aes(log2FoldChange, -log10(padj), label = ID, color = Regulation))+
  geom_point_rast(size = 1)+
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,2))+
  #scale_y_continuous(limits = c(0,25), breaks = seq(0,25,5))+
  scale_color_manual(values = c("down" = col_down, "up" = col_up, "insignificant" = col_neutral))+
  geom_hline(yintercept = -log10(p.val), linetype = 2)+
  geom_vline(xintercept = c(-fc.limit,fc.limit), linetype = 2)+
  ggtitle("DTG in Polysome-enriched samples under NAA vs Ctrl conditions")+
  xlab("Log2 Foldchange")+
  ylab("-Log10 adjusted p-value")+
  annotate("text", x = -4, y = 115, label = length(subset(plotdata_Polysome_Fractions, padj <= p.val & log2FoldChange <= -fc.limit)$ID), size = 8, color = col_down)+
  annotate("text", x = 4, y = 115, label = length(subset(plotdata_Polysome_Fractions, padj <= p.val & log2FoldChange >= fc.limit)$ID), size = 8, color = col_up)+
  theme_light(base_size = 9)+
  theme(axis.text = element_text(color = "black"))

#Save volcano plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/Volcano_Polysome_Fractions.pdf", width = 12, height = 12, units = "cm")


#Plot log2FoldChanges from Polysome enriched samples against Total_RNA samples in a scatterplot ________________________________________________________________________

#Create new column in genes.sig data.frame called "ID" and contains rownames of previous data.frame
genes.sig_Total_RNA$ID <- rownames(genes.sig_Total_RNA)
genes.sig_Polysome_Fractions$ID <- rownames(genes.sig_Polysome_Fractions)
#merge the two data.frames on that common "ID" column
log2foldchanges_combined <- merge(genes.sig_Total_RNA, genes.sig_Polysome_Fractions, by = "ID")
#Plot data in a scatterplot 
ggplot(data = log2foldchanges_combined, aes(x = log2foldchanges_combined$log2FoldChange.x, y = log2foldchanges_combined$log2FoldChange.y))+
  geom_point()+
  labs(titel = "Scatterplot of log2FoldChanges from Polysome enriched and Total_RNA samples", x ="Log2FoldChanges Total_RNA samples", y = "Log2FoldChanges Polysome enriched samples")+
  theme_minimal()
#Save volcano plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/Scatterplot_Total_RNA_vs_Polysomes_sigDEGs.pdf", width = 12, height = 12, units = "cm")

#Plot log2FoldChanges from Polysome enriched samples against Total_RNA samples in a scatterplot but for all genes, not just significant DEGs 
unfiltered_res_Total_RNA <- res_Total_RNA_df
unfiltered_res_Total_RNA$ID <- rownames(unfiltered_res_Total_RNA)
unfiltered_res_Polysome_fractions <- res_Polysome_Fractions_df
unfiltered_res_Polysome_fractions$ID <- rownames(unfiltered_res_Polysome_fractions)
unfiltered_res_combined <- merge(unfiltered_res_Total_RNA, unfiltered_res_Polysome_fractions, by = "ID")
#Plot data in a scatterplot 
ggplot(data = unfiltered_res_combined, aes(x = unfiltered_res_combined$log2FoldChange.x, y = unfiltered_res_combined$log2FoldChange.y))+
  geom_point()+
  labs(titel = "Scatterplot of log2FoldChanges from Polysome enriched and Total_RNA samples", x ="Log2FoldChanges Total_RNA samples", y = "Log2FoldChanges Polysome enriched samples")+
  theme_minimal()
#Save volcano plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/Scatterplot_Total_RNA_vs_Polysomes_unfiltered.pdf", width = 12, height = 12, units = "cm")

#Building Venn diagrams comparing NAA treatment vs Ctrl in Total_RNA and Polysome enriched samples 

#Create subsets of results table after DESeq for up- or downregulated genes 
Venn.Total_RNA.up <- rownames(subset(res_Total_RNA_df, res_Total_RNA_df$log2FoldChange >= fc.limit & res_Total_RNA_df$padj <= p.val))
Venn.Total_RNA.down <- rownames(subset(res_Total_RNA_df, res_Total_RNA_df$log2FoldChange <= -fc.limit & res_Total_RNA_df$padj <= p.val)) 
Venn.Polysome.up <- rownames(subset(res_Polysome_Fractions_df, res_Polysome_Fractions_df$log2FoldChange >= fc.limit & res_Polysome_Fractions_df$padj <= p.val))
Venn.Polysome.down <- rownames(subset(res_Polysome_Fractions_df, res_Polysome_Fractions_df$log2FoldChange <= -fc.limit & res_Polysome_Fractions_df$padj <= p.val))
Venn.NAA.up <- rownames(subset(res_condition_NAA_df, res_condition_NAA_df$log2FoldChange >= fc.limit & res_condition_NAA_df$padj <= p.val))
Venn.NAA.down <- rownames(subset(res_condition_NAA_df, res_condition_NAA_df$log2FoldChange <= -fc.limit & res_condition_NAA_df$padj <= p.val))
Venn.Ctrl.up <- rownames(subset(res_condition_Ctrl_df, res_condition_Ctrl_df$log2FoldChange >= fc.limit & res_condition_Ctrl_df$padj <= p.val))
Venn.Ctrl.down <- rownames(subset(res_condition_Ctrl_df, res_condition_Ctrl_df$log2FoldChange <= -fc.limit & res_condition_Ctrl_df$padj <= p.val))


#Between Polysome Fractions and Total_RNA samples 
venn.diagram(list("Total_RNA_up" = Venn.Total_RNA.up, "Total_RNA_down" = Venn.Total_RNA.down, "Polysome_up" = Venn.Polysome.up, "Polysome_down" = Venn.Polysome.down),
             filename = "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/Venn_Polysomes_vs_Total_RNA.png",
             fill = c("#FFB178FF", "#A56457FF", "#C7A2B6FF", "#98768EFF"),
             cat.pos = c(-10,9,0,0),
             disable.logging = TRUE,
             imagetype = "png")

venn.diagram(list("Total_RNA_up" = Venn.Total_RNA.up, "Polysome_up" = Venn.Polysome.up),
             filename = "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/Venn_Polysomes_vs_Total_RNA_up.png",
             fill = c("#FFB178FF", "#C7A2B6FF"),
             cat.pos = c(-10,20),
             disable.logging = TRUE,
             imagetype = "png")


venn.diagram(list("Total_RNA_down" = Venn.Total_RNA.down, "Polysome_down" = Venn.Polysome.down),
             filename = "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/Venn_Polysomes_vs_Total_RNA_down.png",
             fill = c("#A56457FF", "#98768EFF"),
             cat.pos = c(-10,30),
             disable.logging = TRUE,
             imagetype = "png")

#Between NAA vs. Ctrl conditions 
venn.diagram(list("Ctrl_up" = Venn.Ctrl.up, "Ctrl_down" = Venn.Ctrl.down, "NAA_up" = Venn.NAA.up, "NAA_down" = Venn.NAA.down),
             filename = "Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/Venn_NAA_vs_Ctrl.png",
             fill = c("#6F9969FF", "#EFC86EFF", "#454A74FF", "#97C684FF"),
             cat.pos = c(-5,5,-2,2),
             disable.logging = TRUE,
             imagetype = "png")

#Gene Ontology (GO) Enrichment Analysis ############################################################################################################

#Function to prepare named and sorted gene lists for GO enrichment from DESeq2 results object
prep.geneList <- function(x, fc.limit = 1, padj = 0.05) {
  list.tmp <- subset(data.frame(x), abs(log2FoldChange) >= fc.limit & padj <= p.val)$log2FoldChange
  names(list.tmp) <- rownames(subset(x, abs(log2FoldChange) >= fc.limit & padj <= p.val))
  list.tmp <- sort(list.tmp, decreasing = TRUE)
  return(list.tmp)
}

#Prepare gene list on the fly with the custom function above - structure example 
#geneList <- prep.geneList(res_X_df)

#Prepare gene list for data set: Polysomes NAA vs. Ctrl ______________________________________________________________________________________ 
geneList_P <- prep.geneList(res_Polysome_Fractions_df) 

#GO ENRICHMENT of DEGs in Polysome enriched samples (NAA vs. Ctrl) ___________________________________________________________________________ 
Morgenstern <- met.brewer("Morgenstern")

gse_P_all <- gseGO(geneList=geneList_P, 
             ont ="ALL", 
             keyType = "TAIR", 
             #nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.At.tair.db, 
             pAdjustMethod = "none")

#Plot data for Polysome fractions (NAA vs. Ctrl) in a Splitplots _________________________________________________________________________________
dotplot(gse_P_all, showCategory=12, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 9) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/SplitplotGOall_Polysome_Fractions_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

gse_P_BP <- gseGO(geneList=geneList_P, 
                   ont ="BP", 
                   keyType = "TAIR", 
                   #nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = org.At.tair.db, 
                   pAdjustMethod = "none")

#Plot data
dotplot(gse_P_BP, showCategory=12, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 9) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/SplitplotGOBP_Polysome_Fractions_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

gse_P_MF <- gseGO(geneList=geneList_P, 
                  ont ="MF", 
                  keyType = "TAIR", 
                  #nPerm = 10000, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.At.tair.db, 
                  pAdjustMethod = "none")
#Plot data
dotplot(gse_P_MF, showCategory=12, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 9) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/SplitplotGOMF_Polysome_Fractions_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

gse_P_CC <- gseGO(geneList=geneList_P, 
                  ont ="CC", 
                  keyType = "TAIR", 
                  #nPerm = 10000, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.At.tair.db, 
                  pAdjustMethod = "none")
#Plot data
dotplot(gse_P_CC, showCategory=12, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 9) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/SplitplotGOCC_Polysome_Fractions_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

#Plot data in category network plot 
cnetplot(gse_P_all, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_T),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 10)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")+
  labs(x = NULL, y = NULL)+
  theme(legend.position = "right", axis.text = element_text(color = "black")) 
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/cnetGOall_Polysomes_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

#Plot data 
cnetplot(gse_P_BP, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_T),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 10)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")+
  labs(x = NULL, y = NULL)+
  theme(legend.position = "right", axis.text = element_text(color = "black")) 
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/cnetGOBP_Polysomes_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

#Plot data 
cnetplot(gse_P_MF, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_T),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 10)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")+
  labs(x = NULL, y = NULL)+
  theme(legend.position = "right", axis.text = element_text(color = "black")) 
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/cnetGOMF_Polysomes_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

#Plot data 
cnetplot(gse_P_CC, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_T),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 10)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")+
  labs(x = NULL, y = NULL)+
  theme(legend.position = "right", axis.text = element_text(color = "black")) 
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/cnetGOCC_Polysomes_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

#Plot data in Ridgeplot
ridgeplot(gse_P_all)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/RidgeplotGOall_Polysomes_NAA_vs_Ctrl.pdf", width = 15, height = 17, units = "cm")

#Plot data 
ridgeplot(gse_P_BP)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/RidgeplotGOBP_Polysomes_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

#Plot data 
ridgeplot(gse_P_MF)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/RidgeplotGOMF_Polysomes_NAA_vs_Ctrl.pdf", width = 10, height = 10, units = "cm")

#Plot data 
ridgeplot(gse_P_CC)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/RidgeplotGOCC_Polysomes_NAA_vs_Ctrl.pdf", width = 10, height = 10, units = "cm")


#Prepare gene list for data set: Total_RNA NAA vs. Ctrl 
geneList_T <- prep.geneList(res_Total_RNA_df) 

gse_T_all <- gseGO(geneList=geneList_T, 
             ont ="ALL", 
             keyType = "TAIR", 
             #nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.At.tair.db, 
             pAdjustMethod = "none")

#Plot data for Total_RNA samples (comparison NAA vs. Ctrl) in a Splitplot
dotplot(gse_T_all, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 8) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/SplitplotGOall_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

gse_T_BP <- gseGO(geneList=geneList_T, 
                   ont ="BP", 
                   keyType = "TAIR", 
                   #nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = org.At.tair.db, 
                   pAdjustMethod = "none")

#Plot data 
dotplot(gse_T_BP, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 9) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/SplitplotGOBP_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

gse_T_MF <- gseGO(geneList=geneList_T, 
                   ont ="MF", 
                   keyType = "TAIR", 
                   #nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = org.At.tair.db, 
                   pAdjustMethod = "none")

#Plot data 
dotplot(gse_T_MF, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 9) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/SplitplotGOMF_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

gse_T_CC <- gseGO(geneList=geneList_T, 
                   ont ="CC", 
                   keyType = "TAIR", 
                   #nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = org.At.tair.db, 
                   pAdjustMethod = "none")

#Plot data 
dotplot(gse_T_CC, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 9) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/SplitplotGOCC_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

#Plot data for Total_RNA samples (comparison NAA vs. Ctrl) in a category network 
Morgenstern <- met.brewer("Morgenstern")

cnetplot(gse_T_all, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_T),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 10)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")+
  labs(x = NULL, y = NULL)+
  theme(legend.position = "right", axis.text = element_text(color = "black")) 
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/cnetGOall_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

#Plot data 
cnetplot(gse_T_BP, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_T),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 10)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")+
  labs(x = NULL, y = NULL)+
  theme(legend.position = "right", axis.text = element_text(color = "black")) 
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/cnetGOBP_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

#Plot data 
cnetplot(gse_T_MF, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_T),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 10)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")+
  labs(x = NULL, y = NULL)+
  theme(legend.position = "right", axis.text = element_text(color = "black")) 
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/cnetGOMF_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

#Plot data 
cnetplot(gse_T_CC, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_T),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 10)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")+
  labs(x = NULL, y = NULL)+
  theme(legend.position = "right", axis.text = element_text(color = "black")) 
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/cnetGOCC_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

#Plot data in Ridgeplot
ridgeplot(gse_T_all)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/RidgeplotGOall_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 20, units = "cm")

#Plot data 
ridgeplot(gse_T_BP)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/RidgeplotGOBP_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 20, units = "cm")

#Plot data 
ridgeplot(gse_T_MF)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/RidgeplotGOMF_Total_RNA_NAA_vs_Ctrl.pdf", width = 10, height = 10, units = "cm")

#Plot data 
ridgeplot(gse_T_CC)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/RidgeplotGOCC_Total_RNA_NAA_vs_Ctrl.pdf", width = 10, height = 10, units = "cm")

#Prepare gene list for data set: Polysome Fractions Ctrl vs. Total_RNA Ctrl  
geneList_Ctrls <- prep.geneList(res_condition_Ctrl_df)

gse_Ctrl_all <- gseGO(geneList=geneList_Ctrls, 
             ont ="ALL", 
             keyType = "TAIR", 
             #nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.At.tair.db, 
             pAdjustMethod = "none")

#Plot data for Ctrl samples in a Splitplot
dotplot(gse_Ctrl_all, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 8) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/SplitplotGOall_Ctrls.pdf", width = 15, height = 15, units = "cm")

gse_Ctrl_BP <- gseGO(geneList=geneList_Ctrls, 
                      ont ="BP", 
                      keyType = "TAIR", 
                      #nPerm = 10000, 
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = 0.05, 
                      verbose = TRUE, 
                      OrgDb = org.At.tair.db, 
                      pAdjustMethod = "none")

#Plot data 
dotplot(gse_Ctrl_BP, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 8) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/SplitplotGOBP_Ctrls.pdf", width = 15, height = 15, units = "cm")

gse_Ctrl_MF <- gseGO(geneList=geneList_Ctrls, 
                      ont ="MF", 
                      keyType = "TAIR", 
                      #nPerm = 10000, 
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = 0.05, 
                      verbose = TRUE, 
                      OrgDb = org.At.tair.db, 
                      pAdjustMethod = "none")

#Plot data 
dotplot(gse_Ctrl_MF, showCategory=7, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 8) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/SplitplotGOMF_Ctrls.pdf", width = 17, height = 15, units = "cm")

gse_Ctrl_CC <- gseGO(geneList=geneList_Ctrls, 
                      ont ="CC", 
                      keyType = "TAIR", 
                      #nPerm = 10000, 
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = 0.05, 
                      verbose = TRUE, 
                      OrgDb = org.At.tair.db, 
                      pAdjustMethod = "none")

#Plot data for Ctrl samples in a Splitplot
dotplot(gse_Ctrl_CC, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 8) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/SplitplotGOCC_Ctrls.pdf", width = 15, height = 15, units = "cm")

#Plot data for Ctrl samples in a category network 
Morgenstern <- met.brewer("Morgenstern")

cnetplot(gse_Ctrl_all, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_Ctrls),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 5)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/cnetGOall_Ctrls.pdf", width = 15, height = 15, units = "cm")

#Plot data
cnetplot(gse_Ctrl_BP, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_Ctrls),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 5)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/cnetGOBP_Ctrls.pdf", width = 15, height = 15, units = "cm")

#Plot data
cnetplot(gse_Ctrl_MF, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_Ctrls),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 5)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/cnetGOMF_Ctrls.pdf", width = 15, height = 15, units = "cm")

#Plot data 
cnetplot(gse_Ctrl_CC, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_Ctrls),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 5)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/cnetGOCC_Ctrls.pdf", width = 15, height = 15, units = "cm")

#Plot data in Ridgeplot
ridgeplot(gse_Ctrl_all)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/RidgeplotGOall_Ctrls.pdf", width = 15, height = 25, units = "cm")

#Plot data 
ridgeplot(gse_Ctrl_BP)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/RidgeplotGOBP_Ctrls.pdf", width = 15, height = 30, units = "cm")

#Plot data 
ridgeplot(gse_Ctrl_MF)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/RidgeplotGOMF_Ctrls.pdf", width = 15, height = 30, units = "cm")

#Plot data 
ridgeplot(gse_Ctrl_CC)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/RidgeplotGOCC_Ctrls.pdf", width = 15, height = 30, units = "cm")

