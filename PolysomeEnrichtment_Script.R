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
install_github(repo = "lcalviell/Ribo-seQC")

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
# RNA-Seq of polysome enriched and total RNA samples ##########################################################
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

#Write counts to file
#write.csv(STAR.counts, "Y:/Omics/RiboSeq/PolysomeEnrichment/STAR_Output/STAR.counts.csv", quote = FALSE, row.names = TRUE, sep = "\t")


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
#write.table(sample.table, "./sample_table.txt", quote = FALSE, row.names = TRUE, sep = ",")

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
genes.sig.auxin <- data.frame("Family" = factor(rep(c("SAUR", "Aux/IAA", "GH3", "ARF"), each = 4), levels = c("SAUR", "Aux/IAA", "GH3", "ARF")),
                              "Number" = c(
                                length(which(rowSums(STAR.counts[Gene.Metadata$ID[grep("SAUR", Gene.Metadata$Symbol)],]) < 17)),
                                length(grep("SAUR", Gene.Metadata$Symbol, value = TRUE))-length(grep("SAUR", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE))-length(grep("SAUR", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE))-length(which(rowSums(STAR.counts[Gene.Metadata$ID[grep("SAUR", Gene.Metadata$Symbol)],]) < 17)),
                                length(grep("SAUR", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE)),
                                length(grep("SAUR", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE)),
                                length(which(rowSums(STAR.counts[Gene.Metadata$ID[grep("^IAA| IAA", Gene.Metadata$Symbol)],]) < 17)),
                                length(grep("^IAA| IAA", Gene.Metadata$Symbol, value = TRUE))-length(grep("^IAA| IAA", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE))-length(grep("^IAA| IAA", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE))-length(which(rowSums(STAR.counts[Gene.Metadata$ID[grep("^IAA| IAA", Gene.Metadata$Symbol)],]) < 17)),
                                length(grep("^IAA| IAA", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE)),
                                length(grep("^IAA| IAA", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE)),
                                length(which(rowSums(STAR.counts[Gene.Metadata$ID[grep("GH3\\.", Gene.Metadata$Symbol)],]) < 17)),
                                length(grep("GH3\\.", Gene.Metadata$Symbol, value = TRUE))- length(grep("GH3\\.", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE))-length(grep("GH3\\.", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE))-length(which(rowSums(STAR.counts[Gene.Metadata$ID[grep("GH3\\.", Gene.Metadata$Symbol)],]) < 17)),
                                length(grep("GH3\\.", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE)),
                                length(grep("GH3\\.", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE)),
                                length(which(rowSums(STAR.counts[Gene.Metadata$ID[grep("^ARF[0-9]| ARF[0-9]", Gene.Metadata$Symbol)],]) < 17)),
                                length(grep("^ARF[0-9]| ARF[0-9]", Gene.Metadata$Symbol, value = TRUE))-length(grep("^ARF[0-9]| ARF[0-9]", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE))-length(grep("^ARF[0-9]| ARF[0-9]", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE))-length(which(rowSums(STAR.counts[Gene.Metadata$ID[grep("^ARF[0-9]| ARF[0-9]", Gene.Metadata$Symbol)],]) < 17)),
                                length(grep("^ARF[0-9]| ARF[0-9]", subset(genes.sig.meta, log2FoldChange > 0)$Symbol, value = TRUE)),
                                length(grep("^ARF[0-9]| ARF[0-9]", subset(genes.sig.meta, log2FoldChange < 0)$Symbol, value = TRUE))
                              ),
                              "Regulation" = factor(rep(c("Unexpressed", "None", "Up", "Down"), 4), levels = c("Unexpressed", "None", "Up", "Down"))
)


ggplot(genes.sig.auxin, aes(Family, Number, fill = Regulation))+
  geom_col()+
  scale_fill_manual(values = c("grey90", "grey80", met.brewer("Morgenstern", 8)[7], met.brewer("Morgenstern", 8)[2]))+
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
p.val <- 0.9
fc.limit <- 0

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

#Plot normalized counts for specific genes depending on treatment or control conditions 
#Highest up-regulated transcript under both conditions   
AT2G23170_O <- plotCounts(dds.overall, gene = "AT2G23170", intgroup = c("treatment", "condition"), returnData = TRUE)
ggplot(AT2G23170_O, aes(x = treatment, y = count))+
  geom_point(position = position_jitter(w=0.1, h=0))+
  scale_y_log10()+
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "red")+
  facet_grid(cols = vars(condition))
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/NAA/AT2G23170_interaction.pdf", width = 20, height = 10, units = "cm")

#significant up-regulated transcript exclusively in polysome fractions, but difference to expression in total RNA samples was insignificant 
AT1G14060_O <- plotCounts(dds.overall, gene = "AT1G14060", intgroup = c("treatment", "condition"), returnData = TRUE)
ggplot(AT1G14060_O, aes(x = treatment, y = count))+
  geom_point(position = position_jitter(w=0.1, h=0))+
  scale_y_log10()+
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "red")+
  facet_grid(cols = vars(condition))
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/NAA/AT1G14060_interaction.pdf", width = 20, height = 10, units = "cm")

  
AT4G10910_T <- plotCounts(dds.Total_RNA, gene = "AT4G10910", intgroup = "treatment", returnData = TRUE)

ggplot(AT4G10910_T, aes(x = treatment, y = count))+
  geom_point(position = position_jitter(w=0.1, h=0))
  #scale_y_log10()
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/NAA/AT4G10910_Total_RNA_interaction.pdf", width = 12, height = 12, units = "cm")

AT4G10910_P <- plotCounts(dds.Total_RNA, gene = "AT4G10910", intgroup = "treatment", returnData = TRUE)

ggplot(AT4G10910_P, aes(x = treatment, y = count))+
  geom_point(position = position_jitter(w=0.1, h=0))+
scale_y_log10()
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/NAA/AT4G10910_Total_RNA_interaction.pdf", width = 12, height = 12, units = "cm")  
  
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

#Build data frame containing foldchanges and p-values for all transcripts in total RNA samples and polysome fractions
total_DEGs_Total_RNA_df <- as.data.frame(results(dds.Total_RNA)) 
total_DEGs_Total_RNA_df$ID <- rownames(total_DEGs_Total_RNA_df)
total_DEGs_Polysome_Fractions_df <- as.data.frame(results(dds.Polysome_Fractions))
total_DEGs_Polysome_Fractions_df$ID <- rownames(total_DEGs_Polysome_Fractions_df)
combined_df <- merge(total_DEGs_Total_RNA_df, total_DEGs_Polysome_Fractions_df, by = "ID", all = TRUE)

#Modify dataframe for plotting by adding coloring information based on significance 
combined_df$Significance <- "grey80"  # Initialize all as neutral
combined_df$alpha <- 0
combined_df$Significance[(abs(combined_df$log2FoldChange.x) >= fc.limit & combined_df$padj.x <= p.val)] <- "#FFB178FF"
combined_df$Significance[combined_df$ID %in% rownames(genes.sig_interaction)]<- "magenta"  
combined_df$Significance[(abs(combined_df$log2FoldChange.y) >= fc.limit & combined_df$padj.y <= p.val)] <- "#A56457FF"
combined_df$Significance[
  (abs(combined_df$log2FoldChange.x) >= fc.limit & combined_df$padj.x <= p.val) &
    (abs(combined_df$log2FoldChange.y) >= fc.limit & combined_df$padj.y <= p.val)
] <- "#DFBBC8FF"
combined_df$alpha[combined_df$Significance != "grey80"] <- 1

ggplot(data = combined_df, aes(x = combined_df$log2FoldChange.x, y = combined_df$log2FoldChange.y))+
  geom_point(color = combined_df$Significance, alpha = combined_df$alpha)+
  geom_smooth(method = "lm", se = FALSE, color = "grey")+
  geom_abline(intercept = 0, slope = 1)+
  labs(title = "log2FoldChanges from polysome enriched and total RNA samples", x = "Log2FoldChanges Total_RNA samples", y = "Log2FoldChanges Polysome enriched samples")+
  theme_minimal()+
  theme_light(base_size = 9)
#Save volcano plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/Scatterplot_Total_RNA_vs_Polysomes_unfiltered_interaction.pdf", width = 15, height = 15, units = "cm")

library(gridExtra)
library(grid)

#merging and filtering data frame
genes.sig.meta_Polysome_Fractions <- merge(genes.sig.meta_Polysome_Fractions, Gene.Metadata, all.x = TRUE, all.y = FALSE)
combined_df_meta <- merge(combined_df, Gene.Metadata, all.x = TRUE, all.y = FALSE)

table_excl_in_P <- combined_df_meta[grep("#A56457FF", combined_df_meta$Significance), ]
table_excl_in_P <- table_excl_in_P[order(-table_excl_in_P$log2FoldChange.y), ]
table_excl_in_P <- head(table_excl_in_P, 10) 

genes.sig.grouped <- table_excl_in_P[, c(1,17,9,13,3,7)] 

table_excl_in_T <- combined_df_meta[grep("#FFB178FF", combined_df_meta$Significance), ]
table_excl_in_T <- table_excl_in_T[order(-table_excl_in_T$log2FoldChange.x), ]
table_excl_in_T <- head(table_excl_in_T, 10)

genes.sig.grouped <- rbind(genes.sig.grouped, table_excl_in_T[, c(1,17,9,13,3,7)])

combined_df_meta <- combined_df_meta[order(-combined_df_meta$log2FoldChange.x), ]
combined_df_meta <- head(combined_df_meta, 10)

genes.sig.grouped <- rbind(genes.sig.grouped, combined_df_meta[, c(1,17,9,13,3,7)]) 

#Renaming columns 
colnames(genes.sig.grouped)[colnames(genes.sig.grouped) == "log2FoldChange.y"] <- "log2FoldChange"
colnames(genes.sig.grouped)[colnames(genes.sig.grouped) == "log2FoldChange.x"] <- "total RNA log2FoldChange"
colnames(genes.sig.grouped)[colnames(genes.sig.grouped) == "padj.y"] <- "padj"
colnames(genes.sig.grouped)[colnames(genes.sig.grouped) == "padj.x"] <- "total RNA padj"

# Round the relevant columns to two decimal places
genes.sig.grouped$log2FoldChange <- round(genes.sig.grouped$log2FoldChange, 2)
genes.sig.grouped$padj <- round(genes.sig.grouped$padj, 3)
genes.sig.grouped$`total RNA log2FoldChange`<- round(genes.sig.grouped$`total RNA log2FoldChange`, 2) 
genes.sig.grouped$`total RNA padj` <- round(genes.sig.grouped$`total RNA padj`, 3)

#Renaming columns
colnames(genes.sig.grouped)[colnames(genes.sig.grouped) == "total RNA log2FoldChange"] <- "log2FoldChange"
colnames(genes.sig.grouped)[colnames(genes.sig.grouped) == "total RNA padj"] <- "padj"

#adjust symbol column 
genes.sig.grouped$Symbol <- gsub(",.*$", "", genes.sig.grouped$Symbol)

# Define row colors with transparency 
row_colors <- c(adjustcolor("#A56457FF", alpha.f = 0.8),
                adjustcolor("#FFB178FF", alpha.f = 0.8),
                adjustcolor("#DFBBC8FF", alpha.f = 0.8))

# Repeat the colors for the number of rows needed
row_colors <- rep(row_colors, each = 10)

# Create a table grob with background colors for rows
table_grob <- tableGrob(
  genes.sig.grouped,
  theme = ttheme_default(
    core = list(
      bg_params = list(fill = row_colors),  # Set background colors
      fg_params = list(col = "black")  # Text color
    )
  )
)

#Safe table as PDF file
pdf("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Genes.sig.grouped.pdf", width = 20, height = 20)
grid.draw(table_grob)
dev.off()

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

gse_P_all_simp <- simplify(gse_P_all)

#Plot data for Polysome fractions (NAA vs. Ctrl) in a Splitplots _________________________________________________________________________________
dotplot(gse_P_all_simp, showCategory=12, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 9) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/simpSplitplotGOall_Polysome_Fractions_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

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

gse_P_BP_simp <- simplify(gse_P_BP)

#Plot data
dotplot(gse_P_BP_simp, showCategory=12, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 9) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/PolysomeFractions/simpSplitplotGOBP_Polysome_Fractions_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

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
cnetplot(gse_P_BP_simp, 
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

gse_T_all_simp <- simplify(gse_T_all)

#Plot data for Total_RNA samples (comparison NAA vs. Ctrl) in a Splitplot
dotplot(gse_T_all_simp, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 8) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/simpSplitplotGOall_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

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

gse_T_BP_simp <- simplify(gse_T_BP)

#Plot data 
dotplot(gse_T_BP_simp, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 9) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/simpSplitplotGOBP_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 15, units = "cm")

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
ridgeplot(gse_T_all_simp)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/simpRidgeplotGOall_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 20, units = "cm")

#Plot data 
ridgeplot(gse_T_BP_simp)+
  labs(x = "Enrichment Distribution")+
  theme_light(base_size = 9)+
  scale_fill_gradientn(colors = Morgenstern)+ 
  theme(legend.position = "right", axis.text = element_text(color = "black"))
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/Total_RNASamples/simpRidgeplotGOBP_Total_RNA_NAA_vs_Ctrl.pdf", width = 15, height = 20, units = "cm")

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

gse_Ctrl_all_simp <- simplify(gse_Ctrl_all)

#Plot data for Ctrl samples in a Splitplot
dotplot(gse_Ctrl_all_simp, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 8) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/simpSplitplotGOall_Ctrls.pdf", width = 15, height = 15, units = "cm")

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

gse_Ctrl_BP_simp <- simplify(gse_Ctrl_BP)

#Plot data 
dotplot(gse_Ctrl_BP_simp, showCategory=8, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 8) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/simpSplitplotGOBP_Ctrls.pdf", width = 15, height = 15, units = "cm")

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

gse_Ctrl_MF_simp <- simplify(gse_Ctrl_MF)

#Plot data 
dotplot(gse_Ctrl_MF_simp, showCategory=7, split=".sign") + facet_grid(.~.sign) + theme_light(base_size = 8) + scale_fill_gradientn(colors = Morgenstern)
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/simpSplitplotGOMF_Ctrls.pdf", width = 17, height = 15, units = "cm")

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

cnetplot(gse_Ctrl_all_simp, 
         categorySize="pvalue",
         color.params = list(foldChange=geneList_Ctrls),
         cex.params = list(gene_label = 0.6, category_label = 0.8),
         showCategorie = 5)+
  theme_light(base_size = 8)+
  scale_color_gradientn(colors = Morgenstern, name = "log2FoldChange")
#Save plot to file
ggsave("Y:/Omics/_GitHub_Repositories/Franzi/PolysomeEnrichment/Plots/GO/CtrlSamples/simpcnetGOall_Ctrls.pdf", width = 15, height = 15, units = "cm")

#Plot data
cnetplot(gse_Ctrl_BP_simp, 
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
ridgeplot(gse_Ctrl_BP_simp)+
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

#Build table with significant top 10 transcripts for total RNA samples, polysome fractions and in both 
sorted_genes.sig.meta_Polysome_Fractions <- genes.sig.meta_Polysome_Fractions[order(-genes.sig.meta_Polysome_Fractions$log2FoldChange), ] 
sorted_genes.sig.meta_Polysome_Fractions <- head(sorted_genes.sig.meta_Polysome_Fractions, 10)
table <- sorted_genes.sig.meta_Polysome_Fractions[, c(1,3,7,9)]

sorted_genes.sig.meta_Total_RNA <- genes.sig.meta_Total[order(-genes.sig.meta_Total$log2FoldChange), ]
sorted_genes.sig.meta_Total_RNA <- head(sorted_genes.sig.meta_Total_RNA, 10)
table <- rbind(table, sorted_genes.sig.meta_Total_RNA[, c(1,3,7,9)])  

common_genes <- intersect(genes.sig.meta_Polysome_Fractions$ID, genes.sig.meta_Total$ID)
common_df <- genes.sig.meta_Polysome_Fractions[genes.sig.meta_Polysome_Fractions$ID %in% common_genes, ]

exclusive.sig.genes.meta_Polysome_Fractions <- genes.sig.meta_Polysome_Fractions[!(genes.sig.meta_Polysome_Fractions$ID %in% genes.sig.meta_Total$ID)]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Ribo-Seq analyses to compare occurence of rRNA fragments under different growth conditions  #################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# PREPARE RAW DATA ############################################################################################

#Concatenate files from different lanes into single file
#Rename files to a uniform format if not already done 
#setwd("Y:/Omics/RiboSeq/Growth_Conditions_rRNA/RawData/cleaned") 
#file.rename(dir(), gsub("-SMALLRNA.*Seq", "",dir()))
#Only work with gzipped files for more efficient storage
#Keep all files in a single folder

# QUALITY CONTROL _____________________________________________________________________________________________

#Run FastQC on the read files
setwd("Y:/Omics/RiboSeq/Growth_Conditions_rRNA/RawData/cleaned/") #switch to directory containing the reads
system(paste("wsl fastqc *.fastq.gz --outdir ../FastQC/")) #Run FastQC over all files and save results to specified folder

# SET DIRECTORIES _____________________________________________________________________________________________

#Set linux path to genome directory
genomeDir <- "/mnt/y/Omics/RiboSeq/Growth_Conditions_rRNA/Annotations/STAR_Index"
#Set linux path to genome fasta file
genomeFastaFiles <- "/mnt/y/Omics/RiboSeq/Growth_Conditions_rRNA/Annotations/TAIR10_chr_all.fas"
#Set linux path to annotation GTF file
GTFfile <- "/mnt/y/Omics/RiboSeq/Growth_Conditions_rRNA/Annotations/Araport11_GTF_genes_transposons.20241001.gtf"

# STAR MAPPING ################################################################################################

#Generate Reference Genome 
system(paste("wsl ~/STAR-2.7.11b/source/STAR",
             "--runThreadN 30", #number of threads to use for computation
             "--runMode genomeGenerate", #run mode to generate index
             "--genomeDir", genomeDir, #directory in which to store the index (see above)
             "--genomeFastaFiles", genomeFastaFiles, #directory of the genome fasta (see above)
             "--sjdbGTFfile", GTFfile, #directory of the annotation (see above)
             "--sjdbOverhang 63", # max read length - 1
             "--genomeSAindexNbases 12")) #length of SA pre-indexing string, scaled to Arabidopsis genome according to manual: min(14, log2(GenomeLength)/2 - 1)

#Run STAR mapping in single-end mode 
setwd("Y:/Omics/RiboSeq/Growth_Conditions_rRNA/RawData/cleaned/") #switch to directory containing the reads

for (k in list.files(pattern = ".fastq.gz$")) { 
  system(paste0("wsl echo \"Processing ", k, "\"; ~/STAR-2.7.11b/source/STAR ",
                "--readFilesCommand zcat ",
                "--genomeDir ", genomeDir, " ",
                "--readFilesIn ", k, " ",
                "--outFileNamePrefix ~/STAR_Output/", gsub("fastq\\.gz$", "", k), " ",
                "--quantMode GeneCounts ",
                "--runThreadN 12 ", 
                "--alignSJoverhangMin 8 ", 
                "--alignSJDBoverhangMin 2 ", 
                "--outSAMtype BAM SortedByCoordinate ", 
                "--outSAMmultNmax 1 ", 
                "--outMultimapperOrder Random ", 
                "--outFilterMismatchNmax 1 ", 
                "--outFilterMultimapNmax 20 ", 
                "--outFilterType BySJout ", 
                "--outReadsUnmapped Fastx"))
  print("Moving STAR Output from Linux filesystem to Windows...")
  system("wsl mv ~/STAR_Output/* /mnt/y/Omics/RiboSeq/Growth_Conditions_rRNA/STAR_Output/")
  print("Done")
}


setwd("Y:/Omics/RiboSeq/Growth_Conditions_rRNA") #switch back to original working directory

#MultiQC #####################################################################################################

#run MultiQC to get a summarised report
#run in Ubuntu shell because it's not working with the system("wsl ...") command,  
# cd /mnt/y/Omics/RiboSeq/Growth_Conditions_rRNA/ #set working directory
# multiqc . -o MultiQC/ #performing MultiQC analysis with data from working directory, MultiQC report will be saved in the defined output directory 

# STAR DATA IMPORT ###########################################################################################

#Set working directory of data to be analyzed
setwd("Y:/Omics/RiboSeq/Growth_Conditions_rRNA/")

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

## RIBOSEQC ____________________________________________________________________

### If not already done, get twobit programs and place them in F:\Annotations_RNASeq
### Download faToTwoBit program from: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
### Download twoBitInfo program from: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
### Tutorial: https://genome.ucsc.edu/goldenpath/help/twoBit.html

### PREPARE TWOBIT FILE 
system("wsl /mnt/y/Omics/Annotations_RNASeq/faToTwoBit /mnt/y/Omics/Annotations_RNASeq/TAIR10_chr_all.fas /mnt/y/Omics/Annotations_RNASeq/genome.2bit")
### CHECK TWOBIT FILE INTEGRITY 
system("wsl /mnt/y/Omics/Annotations_RNASeq/twoBitInfo /mnt/y/Omics/Annotations_RNASeq/genome.2bit stdout | sort -k2rn")

### MODIFY PREPARE_ANNOTATION_FILES FUNCTION TO PROPERLY EXTRACT EXON INFORMATION 
prepare_annotation_files <- function (annotation_directory, twobit_file, gtf_file, scientific_name = "Homo.sapiens", 
                                      annotation_name = "genc25", export_bed_tables_TxDb = T, 
                                      forge_BSgenome = T, create_TxDb = T) 
{
  DEFAULT_CIRC_SEQS <- unique(c("chrM", "MT", "MtDNA", "mit", 
                                "Mito", "mitochondrion", "dmel_mitochondrion_genome", 
                                "Pltd", "ChrC", "Pt", "chloroplast", "Chloro", "2micron", 
                                "2-micron", "2uM", "Mt", "NC_001879.2", "NC_006581.1", 
                                "ChrM", "mitochondrion_genome"))
  annotation_name <- gsub(annotation_name, pattern = "_", 
                          replacement = "")
  annotation_name <- gsub(annotation_name, pattern = "-", 
                          replacement = "")
  if (!dir.exists(annotation_directory)) {
    dir.create(path = annotation_directory, recursive = T)
  }
  annotation_directory <- normalizePath(annotation_directory)
  twobit_file <- normalizePath(twobit_file)
  gtf_file <- normalizePath(gtf_file)
  for (f in c(twobit_file, gtf_file)) {
    if (file.access(f, 0) == -1) {
      stop("\n                 The following files don't exist:\n", 
           f, "\n")
    }
  }
  scientific_name_spl <- strsplit(scientific_name, "[.]")[[1]]
  ok <- length(scientific_name_spl) == 2
  if (!ok) {
    stop("\"scientific_name\" must be two words separated by a \".\", like \"Homo.sapiens\"")
  }
  seqinfotwob <- seqinfo(TwoBitFile(twobit_file))
  circss <- seqnames(seqinfotwob)[which(seqnames(seqinfotwob) %in% 
                                          DEFAULT_CIRC_SEQS)]
  seqinfotwob@is_circular[which(seqnames(seqinfotwob) %in% 
                                  DEFAULT_CIRC_SEQS)] <- TRUE
  pkgnm <- paste("BSgenome", scientific_name, annotation_name, 
                 sep = ".")
  circseed <- circss
  if (length(circseed) == 0) {
    circseed <- NULL
  }
  if (forge_BSgenome) {
    cat(paste("Creating the BSgenome package ... ", date(), 
              "\n", sep = ""))
    seed_text <- paste("Package: BSgenome.", scientific_name, 
                       ".", annotation_name, "\n", "Title: Full genome sequences for ", 
                       scientific_name, ", ", annotation_name, "\n", "Description: Full genome sequences for ", 
                       scientific_name, ", ", annotation_name, "\n", "Version: 1.0", 
                       "\n", "organism: ", scientific_name, "\n", "common_name: ", 
                       scientific_name, "\n", "provider: NA", "\n", "provider_version: ", 
                       annotation_name, "\n", "release_date: NA", "\n", 
                       "release_name: NA", "\n", "source_url: NA", "\n", 
                       "organism_biocview: ", scientific_name, "\n", "BSgenomeObjname: ", 
                       scientific_name, "\n", "seqs_srcdir: ", dirname(twobit_file), 
                       "\n", "seqfile_name: ", basename(twobit_file), sep = "")
    seed_dest <- paste(annotation_directory, "/", basename(twobit_file), 
                       "_", scientific_name, "_seed", sep = "")
    if (length(circseed) == 0) {
      writeLines(text = seed_text, con = seed_dest)
    }
    if (length(circseed) == 1) {
      seed_text <- paste(seed_text, "\n", "circ_seqs: \"", 
                         circseed, "\"", sep = "")
      writeLines(text = seed_text, con = seed_dest)
    }
    if (length(circseed) > 1) {
      circseed <- paste("c(\"", paste(circseed, collapse = ","), 
                        "\")", sep = "")
      circseed <- gsub(circseed, pattern = ",", replacement = "\",\"")
      cat(seed_text, "\n", "circ_seqs: ", circseed, "\n", 
          sep = "", file = seed_dest)
    }
    unlink(paste(annotation_directory, pkgnm, sep = "/"), 
           recursive = T)
    forgeBSgenomeDataPkg(x = seed_dest, destdir = annotation_directory, 
                         seqs_srcdir = dirname(twobit_file))
    cat(paste("Creating the BSgenome package --- Done! ", 
              date(), "\n", sep = ""))
    cat(paste("Installing the BSgenome package ... ", date(), 
              "\n", sep = ""))
    install(paste(annotation_directory, pkgnm, sep = "/"), 
            upgrade = F)
    cat(paste("Installing the BSgenome package --- Done! ", 
              date(), "\n", sep = ""))
  }
  if (create_TxDb) {
    cat(paste("Creating the TxDb object ... ", date(), "\n", 
              sep = ""))
    annotation <- makeTxDbFromGFF(file = gtf_file, format = "gtf", 
                                  chrominfo = seqinfotwob)
    saveDb(annotation, file = paste(annotation_directory, 
                                    "/", basename(gtf_file), "_TxDb", sep = ""))
    cat(paste("Creating the TxDb object --- Done! ", date(), 
              "\n", sep = ""))
    cat(paste("Extracting genomic regions ... ", date(), 
              "\n", sep = ""))
    genes <- genes(annotation)
    exons_ge <- exonsBy(annotation, by = "gene")
    exons_ge <- reduce(exons_ge)
    cds_gen <- cdsBy(annotation, "gene")
    cds_ge <- reduce(cds_gen)
    threeutrs <- reduce(GenomicRanges::setdiff(unlist(threeUTRsByTranscript(annotation)), 
                                               unlist(cds_ge), ignore.strand = FALSE))
    fiveutrs <- reduce(GenomicRanges::setdiff(unlist(fiveUTRsByTranscript(annotation)), 
                                              unlist(cds_ge), ignore.strand = FALSE))
    introns <- reduce(GenomicRanges::setdiff(unlist(intronsByTranscript(annotation)), 
                                             unlist(exons_ge), ignore.strand = FALSE))
    nc_exons <- reduce(GenomicRanges::setdiff(unlist(exons_ge), 
                                              reduce(c(unlist(cds_ge), fiveutrs, threeutrs)), 
                                              ignore.strand = FALSE))
    ov <- findOverlaps(threeutrs, genes)
    ov <- split(subjectHits(ov), queryHits(ov))
    threeutrs$gene_id <- CharacterList(lapply(ov, FUN = function(x) {
      names(genes)[x]
    }))
    ov <- findOverlaps(fiveutrs, genes)
    ov <- split(subjectHits(ov), queryHits(ov))
    fiveutrs$gene_id <- CharacterList(lapply(ov, FUN = function(x) {
      names(genes)[x]
    }))
    ov <- findOverlaps(introns, genes)
    ov <- split(subjectHits(ov), queryHits(ov))
    introns$gene_id <- CharacterList(lapply(ov, FUN = function(x) {
      names(genes)[x]
    }))
    ov <- findOverlaps(nc_exons, genes)
    ov <- split(subjectHits(ov), queryHits(ov))
    nc_exons$gene_id <- CharacterList(lapply(ov, FUN = function(x) {
      names(genes)[x]
    }))
    intergenicRegions <- genes
    strand(intergenicRegions) <- "*"
    intergenicRegions <- gaps(reduce(intergenicRegions))
    intergenicRegions <- intergenicRegions[strand(intergenicRegions) == 
                                             "*"]
    cds_tx <- cdsBy(annotation, "tx", use.names = T)
    txs_gene <- transcriptsBy(annotation, by = "gene")
    genes_red <- reduce(sort(genes(annotation)))
    exons_tx <- exonsBy(annotation, "tx", use.names = T)
    transcripts_db <- transcripts(annotation)
    intron_names_tx <- intronsByTranscript(annotation, use.names = T)
    nsns <- exonicParts(annotation)
    exsss_cds <- exons_tx[names(cds_tx)]
    chunks <- seq(1, length(cds_tx), by = 20000)
    if (chunks[length(chunks)] < length(cds_tx)) {
      chunks <- c(chunks, length(cds_tx))
    }
    mapp <- GRangesList()
    for (i in 1:(length(chunks) - 1)) {
      if (i != (length(chunks) - 1)) {
        mapp <- suppressWarnings(c(mapp, pmapToTranscripts(cds_tx[chunks[i]:(chunks[i + 
                                                                                      1] - 1)], transcripts = exsss_cds[chunks[i]:(chunks[i + 
                                                                                                                                            1] - 1)])))
      }
      if (i == (length(chunks) - 1)) {
        mapp <- suppressWarnings(c(mapp, pmapToTranscripts(cds_tx[chunks[i]:(chunks[i + 
                                                                                      1])], transcripts = exsss_cds[chunks[i]:(chunks[i + 
                                                                                                                                        1])])))
      }
    }
    cds_txscoords <- unlist(mapp)
    cat(paste("Extracting ids and biotypes ... ", date(), 
              "\n", sep = ""))
    trann <- unique(mcols(import.gff2(gtf_file, colnames = c("gene_id", 
                                                             "gene_biotype", "gene_type", "gene_name", "gene_symbol", 
                                                             "transcript_id", "transcript_biotype", "transcript_type"))))
    trann <- trann[!is.na(trann$transcript_id), ]
    trann <- data.frame(unique(trann), stringsAsFactors = F)
    if (sum(!is.na(trann$transcript_biotype)) == 0 & sum(!is.na(trann$transcript_type)) == 
        0) {
      trann$transcript_biotype <- "no_type"
    }
    if (sum(!is.na(trann$transcript_biotype)) == 0) {
      trann$transcript_biotype <- NULL
    }
    if (sum(!is.na(trann$transcript_type)) == 0) {
      trann$transcript_type <- NULL
    }
    if (sum(!is.na(trann$gene_biotype)) == 0 & sum(!is.na(trann$gene_type)) == 
        0) {
      trann$gene_type <- "no_type"
    }
    if (sum(!is.na(trann$gene_name)) == 0 & sum(!is.na(trann$gene_symbol)) == 
        0) {
      trann$gene_name <- "no_name"
    }
    if (sum(!is.na(trann$gene_biotype)) == 0) {
      trann$gene_biotype <- NULL
    }
    if (sum(!is.na(trann$gene_type)) == 0) {
      trann$gene_type <- NULL
    }
    if (sum(!is.na(trann$gene_name)) == 0) {
      trann$gene_name <- NULL
    }
    if (sum(!is.na(trann$gene_symbol)) == 0) {
      trann$gene_symbol <- NULL
    }
    colnames(trann) <- c("gene_id", "gene_biotype", "gene_name", 
                         "transcript_id", "transcript_biotype")
    trann <- DataFrame(trann)
    unq_intr <- sort(unique(unlist(intron_names_tx)))
    names(unq_intr) <- NULL
    all_intr <- unlist(intron_names_tx)
    ov <- findOverlaps(unq_intr, all_intr, type = "equal")
    ov <- split(subjectHits(ov), queryHits(ov))
    a_nam <- CharacterList(lapply(ov, FUN = function(x) {
      unique(names(all_intr)[x])
    }))
    unq_intr$type = "J"
    unq_intr$tx_name <- a_nam
    mat_genes <- match(unq_intr$tx_name, trann$transcript_id)
    g <- unlist(apply(cbind(1:length(mat_genes), Y = elementNROWS(mat_genes)), 
                      FUN = function(x) rep(x[1], x[2]), MARGIN = 1))
    g2 <- split(trann[unlist(mat_genes), "gene_id"], g)
    unq_intr$gene_id <- CharacterList(lapply(g2, unique))
    ncrnas <- nc_exons[!nc_exons %over% genes[trann$gene_id[trann$gene_biotype == 
                                                              "protein_coding"]]]
    ncisof <- nc_exons[nc_exons %over% genes[trann$gene_id[trann$gene_biotype == 
                                                             "protein_coding"]]]
    ifs <- seqinfo(annotation)
    translations <- as.data.frame(ifs)
    translations$genetic_code <- "1"
    translations$genetic_code[rownames(translations) %in% 
                                c("chrM", "MT", "MtDNA", "mit", "mitochondrion")] <- "2"
    translations$genetic_code[rownames(translations) %in% 
                                c("Mito")] <- "3"
    translations$genetic_code[rownames(translations) %in% 
                                c("dmel_mitochondrion_genome")] <- "5"
    circs <- ifs@seqnames[which(ifs@is_circular)]
    suppressPackageStartupMessages(library(pkgnm, character.only = TRUE))
    genome <- get(pkgnm)
    tocheck <- as.character(runValue(seqnames(cds_tx)))
    tocheck <- cds_tx[!tocheck %in% circs]
    seqcds <- extractTranscriptSeqs(genome, transcripts = tocheck)
    cd <- unique(translations$genetic_code[!rownames(translations) %in% 
                                             circs])
    trsl <- suppressWarnings(translate(seqcds, genetic.code = getGeneticCode(cd), 
                                       if.fuzzy.codon = "solve"))
    trslend <- as.character(narrow(trsl, end = width(trsl), 
                                   width = 1))
    stop_inannot <- NA
    if (names(sort(table(trslend), decreasing = T)[1]) == 
        "*") {
      stop_inannot <- "*"
    }
    cds_txscoords$gene_id <- trann$gene_id[match(as.vector(seqnames(cds_txscoords)), 
                                                 trann$transcript_id)]
    cds_cc <- cds_txscoords
    strand(cds_cc) <- "*"
    sta_cc <- resize(cds_cc, width = 1, "start")
    sta_cc <- unlist(pmapFromTranscripts(sta_cc, exons_tx[seqnames(sta_cc)], 
                                         ignore.strand = F))
    sta_cc$gene_id <- trann$gene_id[match(names(sta_cc), 
                                          trann$transcript_id)]
    sta_cc <- sta_cc[sta_cc$hit]
    strand(sta_cc) <- structure(as.vector(strand(transcripts_db)), 
                                names = transcripts_db$tx_name)[names(sta_cc)]
    sta_cc$type <- "start_codon"
    mcols(sta_cc) <- mcols(sta_cc)[, c("exon_rank", "type", 
                                       "gene_id")]
    sto_cc <- resize(cds_cc, width = 1, "end")
    sto_cc <- shift(sto_cc, -2)
    if (is.na(stop_inannot)) {
      sto_cc <- resize(trim(shift(sto_cc, 3)), width = 1, 
                       fix = "end")
    }
    sto_cc <- unlist(pmapFromTranscripts(sto_cc, exons_tx[seqnames(sto_cc)], 
                                         ignore.strand = F))
    sto_cc <- sto_cc[sto_cc$hit]
    sto_cc$gene_id <- trann$gene_id[match(names(sto_cc), 
                                          trann$transcript_id)]
    strand(sto_cc) <- structure(as.vector(strand(transcripts_db)), 
                                names = transcripts_db$tx_name)[names(sto_cc)]
    sto_cc$type <- "stop_codon"
    mcols(sto_cc) <- mcols(sto_cc)[, c("exon_rank", "type", 
                                       "gene_id")]
    cat(paste("Defining most common start/stop codons ... ", 
              date(), "\n", sep = ""))
    start_stop_cc <- sort(c(sta_cc, sto_cc))
    start_stop_cc$transcript_id <- names(start_stop_cc)
    start_stop_cc$most_up_downstream <- FALSE
    start_stop_cc$most_frequent <- FALSE
    df <- cbind.DataFrame(start(start_stop_cc), start_stop_cc$type, 
                          start_stop_cc$gene_id)
    colnames(df) <- c("start_pos", "type", "gene_id")
    upst <- by(df$start_pos, INDICES = df$gene_id, function(x) {
      x == min(x) | x == max(x)
    })
    start_stop_cc$most_up_downstream <- unlist(upst[unique(df$gene_id)])
    mostfr <- by(df[, c("start_pos", "type")], INDICES = df$gene_id, 
                 function(x) {
                   mfreq <- table(x)
                   x$start_pos %in% as.numeric(names(which(mfreq[, 
                                                                 1] == max(mfreq[, 1])))) | x$start_pos %in% 
                     as.numeric(names(which(mfreq[, 2] == max(mfreq[, 
                                                                    2]))))
                 })
    start_stop_cc$most_frequent <- unlist(mostfr[unique(df$gene_id)])
    names(start_stop_cc) <- NULL
    mostupstr_tx <- sum(LogicalList(split(start_stop_cc$most_up_downstream, 
                                          start_stop_cc$transcript_id)))[as.character(seqnames(cds_txscoords))]
    cds_txscoords$upstr_stasto <- mostupstr_tx
    mostfreq_tx <- sum(LogicalList(split(start_stop_cc$most_frequent, 
                                         start_stop_cc$transcript_id)))[as.character(seqnames(cds_txscoords))]
    cds_txscoords$mostfreq_stasto <- mostfreq_tx
    cds_txscoords$lentx <- sum(width(exons_tx[as.character(seqnames(cds_txscoords))]))
    df <- cbind.DataFrame(as.character(seqnames(cds_txscoords)), 
                          width(cds_txscoords), start(cds_txscoords), cds_txscoords$mostfreq_stasto, 
                          cds_txscoords$gene_id)
    colnames(df) <- c("txid", "cdslen", "utr5len", "var", 
                      "gene_id")
    repres_freq <- by(df[, c("txid", "cdslen", "utr5len", 
                             "var")], df$gene_id, function(x) {
                               x <- x[order(x$var, x$utr5len, x$cdslen, decreasing = T), 
                               ]
                               x <- x[x$var == max(x$var), ]
                               ok <- x$txid[which(x$cdslen == max(x$cdslen) & x$utr5len == 
                                                    max(x$utr5len) & x$var == max(x$var))][1]
                               if (length(ok) == 0 | is.na(ok[1])) {
                                 ok <- x$txid[1]
                               }
                               ok
                             })
    df <- cbind.DataFrame(as.character(seqnames(cds_txscoords)), 
                          width(cds_txscoords), start(cds_txscoords), cds_txscoords$upstr_stasto, 
                          cds_txscoords$gene_id)
    colnames(df) <- c("txid", "cdslen", "utr5len", "var", 
                      "gene_id")
    repres_upstr <- by(df[, c("txid", "cdslen", "utr5len", 
                              "var")], df$gene_id, function(x) {
                                x <- x[order(x$var, x$utr5len, x$utr5len, decreasing = T), 
                                ]
                                x <- x[x$var == max(x$var), ]
                                ok <- x$txid[which(x$cdslen == max(x$cdslen) & x$utr5len == 
                                                     max(x$utr5len) & x$var == max(x$var))][1]
                                if (length(ok) == 0 | is.na(ok[1])) {
                                  ok <- x$txid[1]
                                }
                                ok
                              })
    df <- cbind.DataFrame(as.character(seqnames(cds_txscoords)), 
                          width(cds_txscoords), start(cds_txscoords), cds_txscoords$upstr_stasto, 
                          cds_txscoords$gene_id)
    colnames(df) <- c("txid", "cdslen", "utr5len", "var", 
                      "gene_id")
    repres_len5 <- by(df[, c("txid", "cdslen", "utr5len", 
                             "var")], df$gene_id, function(x) {
                               x <- x[order(x$utr5len, x$var, x$cdslen, decreasing = T), 
                               ]
                               ok <- x$txid[which(x$utr5len == max(x$utr5len) & 
                                                    x$var == max(x$var))][1]
                               if (length(ok) == 0 | is.na(ok[1])) {
                                 ok <- x$txid[1]
                               }
                               ok
                             })
    cds_txscoords$reprentative_mostcommon <- as.character(seqnames(cds_txscoords)) %in% 
      unlist(repres_freq)
    cds_txscoords$reprentative_boundaries <- as.character(seqnames(cds_txscoords)) %in% 
      unlist(repres_upstr)
    cds_txscoords$reprentative_5len <- as.character(seqnames(cds_txscoords)) %in% 
      unlist(repres_len5)
    unq_stst <- start_stop_cc
    mcols(unq_stst) <- NULL
    unq_stst <- sort(unique(unq_stst))
    ov <- findOverlaps(unq_stst, start_stop_cc, type = "equal")
    ov <- split(subjectHits(ov), queryHits(ov))
    unq_stst$type <- CharacterList(lapply(ov, FUN = function(x) {
      unique(start_stop_cc$type[x])
    }))
    unq_stst$transcript_id <- CharacterList(lapply(ov, FUN = function(x) {
      start_stop_cc$transcript_id[x]
    }))
    unq_stst$gene_id <- CharacterList(lapply(ov, FUN = function(x) {
      unique(start_stop_cc$gene_id[x])
    }))
    unq_stst$reprentative_mostcommon <- sum(!is.na(match(unq_stst$transcript_id, 
                                                         unlist(as(repres_freq, "CharacterList"))))) > 0
    unq_stst$reprentative_boundaries <- sum(!is.na(match(unq_stst$transcript_id, 
                                                         unlist(as(repres_upstr, "CharacterList"))))) > 0
    unq_stst$reprentative_5len <- sum(!is.na(match(unq_stst$transcript_id, 
                                                   unlist(as(repres_len5, "CharacterList"))))) > 0
    GTF_annotation <- list(transcripts_db, txs_gene, ifs, 
                           unq_stst, cds_tx, intron_names_tx, cds_gen, exons_tx, 
                           nsns, unq_intr, genes, threeutrs, fiveutrs, ncisof, 
                           ncrnas, introns, intergenicRegions, trann, cds_txscoords, 
                           translations, pkgnm, stop_inannot)
    names(GTF_annotation) <- c("txs", "txs_gene", "seqinfo", 
                               "start_stop_codons", "cds_txs", "introns_txs", "cds_genes", 
                               "exons_txs", "exons_bins", "junctions", "genes", 
                               "threeutrs", "fiveutrs", "ncIsof", "ncRNAs", "introns", 
                               "intergenicRegions", "trann", "cds_txs_coords", 
                               "genetic_codes", "genome_package", "stop_in_gtf")
    txs_all <- unique(GTF_annotation$trann$transcript_id)
    txs_exss <- unique(names(GTF_annotation$exons_txs))
    txs_notok <- txs_all[!txs_all %in% txs_exss]
    if (length(txs_notok) > 0) {
      set.seed(666)
      cat(paste("Warning: ", length(txs_notok), " txs with incorrect/unspecified exon boundaries - e.g. trans-splicing events, examples: ", 
                paste(txs_notok[sample(1:length(txs_notok), 
                                       size = min(3, length(txs_notok)), replace = F)], 
                      collapse = ", "), " - ", date(), "\n", sep = ""))
    }
    save(GTF_annotation, file = paste(annotation_directory, 
                                      "/", basename(gtf_file), "_Rannot", sep = ""))
    cat(paste("Rannot object created!   ", date(), "\n", 
              sep = ""))
    if (export_bed_tables_TxDb == T) {
      cat(paste("Exporting annotation tables ... ", date(), 
                "\n", sep = ""))
      for (bed_file in c("fiveutrs", "threeutrs", "ncIsof", 
                         "ncRNAs", "introns", "cds_txs_coords")) {
        bf <- GTF_annotation[[bed_file]]
        bf_t <- bf
        if (length(bf) > 0) {
          bf_t <- data.frame(chromosome = seqnames(bf), 
                             start = start(bf), end = end(bf), name = ".", 
                             score = width(bf), strand = strand(bf))
          meccole <- mcols(bf)
          for (mecc in names(meccole)) {
            if (is(meccole[, mecc], "CharacterList") | 
                is(meccole[, mecc], "NumericList") | is(meccole[, 
                                                                mecc], "IntegerList")) {
              meccole[, mecc] <- paste(meccole[, mecc], 
                                       collapse = ";")
            }
          }
          bf_t <- cbind.data.frame(bf_t, meccole)
        }
        write.table(bf_t, file = paste(annotation_directory, 
                                       "/", bed_file, "_similbed.bed", sep = ""), 
                    sep = "\t", quote = FALSE, row.names = FALSE, 
                    col.names = F)
      }
      write.table(GTF_annotation$trann, file = paste(annotation_directory, 
                                                     "/table_gene_tx_IDs", sep = ""), sep = "\t", 
                  quote = FALSE, row.names = FALSE)
      seqi <- as.data.frame(GTF_annotation$seqinfo)
      seqi$chromosome <- rownames(seqi)
      write.table(seqi, file = paste(annotation_directory, 
                                     "/seqinfo", sep = ""), sep = "\t", quote = FALSE, 
                  row.names = FALSE)
      gen_cod <- as.data.frame(GTF_annotation$genetic_codes)
      gen_cod$chromosome <- rownames(gen_cod)
      write.table(gen_cod, file = paste(annotation_directory, 
                                        "/genetic_codes", sep = ""), sep = "\t", quote = FALSE, 
                  row.names = FALSE)
      cat(paste("Exporting annotation tables --- Done! ", 
                date(), "\n", sep = ""))
    }
  }
}

### SET DIRECTORY 
setwd("Y:/Omics/RiboSeq/Growth_Conditions_rRNA/") 

### PREPARE ANNOTATION FILE 
prepare_annotation_files(annotation_directory = "Y:/Omics/RiboSeq/Growth_Conditions_rRNA/Annotations",
                         twobit_file = "Y:/Omics/Annotations_RNASeq/genome.2bit",
                         gtf_file = "Y:/Omics/RiboSeq/Growth_Conditions_rRNA/Annotations/Araport11_GTF_genes_transposons.20241001.gtf",
                         scientific_name = "Arabidopsis.thaliana",
                         annotation_name = "Araport11")

### RUN RIBOSEQC 
### BUILD PATH TO BAM FILES 
bamFiles <- dir(path = "Y:/Omics/RiboSeq/Growth_Conditions_rRNA/STAR_Output/",
                pattern = ".*Coord.out.bam",
                full.names = TRUE)

RiboseQC_analysis(annotation_file = "Y:/Omics/RiboSeq/Growth_Conditions_rRNA/Annotations/Araport11_GTF_genes_transposons.20241001.gtf_Rannot",
                  bam_files = bamFiles,
                  write_tmp_files = FALSE,
                  sample_names = gsub("\\.Aligned\\..*bam$","", bamFiles),
                  extended_report = FALSE,
                  stranded = FALSE)

#create HTML report 
setwd("Y:/Omics/RiboSeq/Growth_Conditions_rRNA/STAR_Output/")
create_html_report(input_files = list.files(pattern = "results", path = "."),
                   input_sample_names = c("3wLR1","3wLR2", "2wSR1", "2wSR2", "8dLR1", "8dLR2", "8dSR1","8dSR2", "ER1", "ER2", "LmR1", "LmR2", "LpR1", "LpR2", "SmR1", "SmR2", "SpR1", "SpR2"),
                   output_file = "Y:/Omics/RiboSeq/Growth_Conditions_rRNA/STAR_Output/report.html",
                   extended = FALSE)

#filter raw data based on sequence length 
setwd("Y:/Omics/RiboSeq/Growth_Conditions_rRNA/RawData/filtered/")

for (i in list.files(pattern = ".fastq.gz")) {
  system(paste0("wsl zcat ", i, 
                " | seqkit seq -m 20 -M 30 ", 
                "-o Filtered_Output/", i))
}

# Description:
# This R script identifies the 100 most abundant sequences from a set of BAM files,
# counts their occurrences, and outputs the results in a table format. 
# Additionally, it generates a heat map to visualize the relative abundance of these sequences 
# across all samples.
#
# Instructions:
# 1. Set the working directory to the folder containing the BAM files.
# 2. Execute the script to obtain the output table and heat map.
#
# Inputs:
# - BAM files located in the specified working directory.
#
# Outputs:
# - A CSV file containing the counts of the 100 most abundant sequences.
# - A heat map in PNG format showing the abundance of the top sequences across samples.
#
# Requirements:
# - R version 4.0 or higher.
#
# Author: Maik Böhmer
# Date: 16.11.2024


# Load necessary libraries
library(Rsamtools)
library(dplyr)
library(tidyr)
library(GenomicFeatures)
library(rtracklayer)
library(txdbmaker)
library(Biostrings)
library(ggplot2)
library(viridis)     # best. color. palette. evar.
library(ggthemes)    # has a clean theme for ggplot2


# set working directory
setwd("~/RiboSeq_cotaminants")

# Function to count sequences and map to genes in a BAM file
count_sequences <- function(bam_file, gene_ranges) {
  # Open the BAM file
  bam <- scanBam(BamFile(bam_file))
  
  # Extract sequences and mapping positions
  sequences <- as.character(bam[[1]]$seq)  # Convert DNAStringSet to character
  chrom <- bam[[1]]$rname
  start <- bam[[1]]$pos
  end <- start + width(DNAStringSet(sequences)) - 1
  
  # Create a GRanges object for the sequences
  sequences_gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end), strand = "*")
  
  # Map sequences to genes
  hits <- findOverlaps(sequences_gr, gene_ranges, ignore.strand = TRUE)
  gene_ids <- as.character(mcols(gene_ranges)[subjectHits(hits), "gene_id"])
  
  # Create a tibble with sequences and gene IDs
  sequence_counts_tbl <- tibble(Sequence = sequences[queryHits(hits)], Gene = gene_ids)
  
  # Count the appearances of each sequence
  sequence_counts_tbl <- sequence_counts_tbl %>%
    group_by(Sequence, Gene) %>%
    summarise(Count = n(), .groups = "drop")
  
  return(sequence_counts_tbl)
}

# Path to the annotation file (e.g., GTF or GFF)
annotation_file <- "Araport11_GTF_genes_transposons.20241001.gtf"

# Create a TxDb object from the annotation file using txdbmaker
txdb <- txdbmaker::makeTxDbFromGFF(annotation_file, format = "gtf")

# Extract gene ranges
gene_ranges <- genes(txdb)

# Get list of BAM files in the current working directory
bam_files <- list.files(pattern = "\\.bam$")

# Initialize an empty list to store sequence counts from each file
sequence_counts_list <- list()

# Loop over each BAM file and get sequence counts
for (bam_file in bam_files) {
  sequence_counts <- count_sequences(bam_file, gene_ranges)
  sequence_counts <- sequence_counts %>%
    mutate(File = bam_file)
  sequence_counts_list[[bam_file]] <- sequence_counts
}

# Combine sequence counts from all files into a single tibble
combined_sequence_counts <- bind_rows(sequence_counts_list)

# Get the top 100 most abundant sequences
top_sequences <- combined_sequence_counts %>%
  group_by(Sequence, Gene) %>%
  summarise(TotalCount = sum(Count), .groups = "drop") %>%
  arrange(desc(TotalCount)) %>%
  slice_head(n = 100)

# Create a comparison table
comparison_table <- top_sequences %>%
  left_join(combined_sequence_counts, by = c("Gene" , "Sequence")) %>%
  ungroup() %>%
  dplyr::select(Gene, Sequence, File, Count)

test <- comparison_table[order(-comparison_table[,2], comparison_table[,1]), ]

# Group by 'sequence' and summarize
comparison_table <- comparison_table %>%
  group_by(Sequence, File) %>%
  summarize(
    Gene = paste(unique(Gene), collapse = ", "),
    Count = first(Count),  # Take the first count (assuming they are identical)
    .groups = 'drop'
  )

# Pivot the table to a wider format for comparison
comparison_table_wide <- comparison_table %>%
  pivot_wider(names_from = File, values_from = Count, values_fill = list(Count = 0))

# Print the comparison table
print(comparison_table_wide)

# Optionally, write the comparison table to a CSV file
write.csv(comparison_table_wide, "comparison_table.csv", row.names = FALSE)

# Create a heat map visualization by sequence
heatmap_data <- comparison_table_wide %>%
  gather(key = "File", value = "Count", -Sequence, -Gene)

ggplot(heatmap_data, aes(x = File, y = Sequence, fill = Count)) +
  geom_tile() +
  facet_grid(Gene~.,scales="free_y", space="free_y") +
  theme(strip.text.y.right = element_text(angle = 0),
        axis.text=element_text(size=6)) +
  scale_fill_viridis(name="Count") +
  labs(title = "Read count distribution of top contaminants among samples", x = "Sample", y = "Fragment")