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
# Authors: Maik Böhmer and Dario Ricciardi
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
setwd("F:/RiboSeq/ClinoNAA_Experiment/Riboseq/STAR_Output/")

# Function to count sequences and map to genes in a BAM file
count_sequences <- function(bam_file, gene_ranges) {
  
  print("Extracting sequences")
  
  # Open the BAM file
  bam <- scanBam(BamFile(bam_file))
  
  # Extract sequences and mapping positions
  seqs <- data.frame("sequences" = as.character(bam[[1]]$seq),
                     "chrom" = bam[[1]]$rname,
                     "start" = bam[[1]]$pos,
                     "end" = bam[[1]]$pos + width(DNAStringSet(as.character(bam[[1]]$seq))) - 1)
  
  print("Counting sequences and fetching gene IDs")
  
  # Count all sequences and only display uniques with their counts
  seqs.count <- count(group_by_all(seqs))
  
  # Create a GRanges object for the sequences
  sequences_gr <- GRanges(seqnames = seqs.count$chrom, ranges = IRanges(start = seqs.count$start, end = seqs.count$end), strand = "*")
  
  # Map sequences to genes
  hits <- findOverlaps(sequences_gr, gene_ranges, ignore.strand = TRUE, select = "first")
  gene_ids <- as.character(mcols(gene_ranges)[hits, "gene_id"])
  
  # Add gene ids to counted sequences
  seqs.count$gene_id <- gene_ids
  seqs.count <- data.frame(seqs.count)
  
  
  # Return data
  return(seqs.count)
}

# Path to the annotation file (e.g., GTF or GFF)
annotation_file <- "X:/Omics/RiboSeq/Growth_Conditions_rRNA/Annotations/Araport11_GTF_genes_transposons.20241001.gtf"

# Create a TxDb object from the annotation file using txdbmaker
txdb <- txdbmaker::makeTxDbFromGFF(annotation_file, format = "gtf")

# Extract gene ranges
gene_ranges <- genes(txdb)

# Get list of BAM files in the current working directory
bam_files <- list.files(pattern = "\\.bam$")[25:32]

# Initialize an empty list to store sequence counts from each file
sequence_counts_list <- list()

# Loop over each BAM file and get sequence counts
for (bam_file in bam_files) {
  print(paste("Working on", bam_file))
  sequence_counts <- count_sequences(bam_file, gene_ranges)
  samplename <- gsub("\\.Aligned.*$","", bam_file)
  sequence_counts$sample <- samplename
  sequence_counts <- sequence_counts[order(sequence_counts$n, decreasing = TRUE),]
  sequence_counts_list[[samplename]] <- sequence_counts
  print("Done")
}


# Build comparison table over all samples

# Combine sequence counts from all files into a single tibble
combined_sequence_counts <- bind_rows(lapply(sequence_counts_list, head, n = 2000))

# Cut down columns
combined_sequence_counts <- combined_sequence_counts[, c(1,5,6,7)]
colnames(combined_sequence_counts) <- c("Sequence", "Count", "Gene_ID", "Sample")
# Reshape to wide format
combined_sequence_counts <- reshape(combined_sequence_counts, idvar = c("Sequence", "Gene_ID"), timevar = "Sample", direction = "wide")

# OPTIONAL: Convert absolute counts to percentage of counts in relation to total counts per sample
for (i in 1:length(sequence_counts_list)){
  combined_sequence_counts[,2+i] <- combined_sequence_counts[,2+i] / sum(sequence_counts_list[[i]]$n) * 100
}

# Replace NA with 0
combined_sequence_counts[is.na(combined_sequence_counts)] <- 0
combined_sequence_counts$Gene_ID[combined_sequence_counts$Gene_ID == 0] <- NA
# Sort by row-wise sums
combined_sequence_counts <- combined_sequence_counts[order(rowSums(combined_sequence_counts[,c(3:length(colnames(combined_sequence_counts)))]), decreasing = TRUE),]
# Beautify colnames
colnames(combined_sequence_counts) <- gsub("^Count\\.","", colnames(combined_sequence_counts))


# Function to summarize the same sequences even if they are mapped to different loci
accumulateFragments <- function(x, accumulate){
  contams <- data.frame()
  # Loop through all sequences from short to long
  for(k in x[order(nchar(x$Sequence)), 1]){
    # Pick out set of sequences that match the shortest one of interest (fuzzy matching if desired)
    if(accumulate == "fuzzy"){
      temp <- x[grepl(substr(k, 2, nchar(k)-1), x$Sequence),]
    }else if(accumulate == "shortest"){
      temp <- x[grepl(k, x$Sequence),]
    }else if (accumulate == "sequence"){
      temp <- x[x$Sequence == k,]
    }else{print("accumulate must be either 'sequence', 'shortest' or 'fuzzy'")}
    # Sum up counts and collapse all mappings into a single character string
    temp <- temp[order(nchar(temp$Sequence)),]
    temp[1,3:length(colnames(temp))] <- colSums(temp[,3:length(colnames(temp))])
    if(length(paste(unique(na.omit(temp$Gene_ID)), sep = ",")) > 0){
      temp$Gene_ID[1] <- paste(sort(unique(na.omit(temp$Gene_ID))), collapse = ",")
    }else temp$Gene_ID[1] <- NA
    
    #Bind temp to contams
    contams <- na.omit(rbind(contams, temp[1,]))
    #Eliminate matches from fragment list
      x <- x[!x$Sequence %in% unique(temp$Sequence),]

  }
  return(contams[order(-rowSums(contams[, 3:length(colnames(contams))])),])
}

#Summarize sequences across loci
# Accumulate counts of the same sequences
# "sequence" to summarize by exact sequence
# "shortest" to summarize all sequences that share a shortest common sequence
# "fuzzy" is the same as shortest but allows for overhangs of 1 nt at each end
accumulated_sequences <- accumulateFragments(combined_sequence_counts, accumulate = "fuzzy")

# OPTIONAL: Add biotype metadata to sequences
Gene.Metadata <- read.table("X:/Omics/RiboSeq/Growth_Conditions_rRNA/Annotations/Arabidopsis_Genes_Metadata.tsv", header = TRUE)
accumulated_sequences$ID <- gsub(",.*$","", accumulated_sequences$Gene_ID)
accumulated_sequences <- merge(accumulated_sequences, Gene.Metadata[,1:2], sort = FALSE, all.x = TRUE)
accumulated_sequences <- accumulated_sequences[,c("Sequence", "Gene_ID", "type", grep("ID|Sequence|Gene_ID|type",colnames(accumulated_sequences), invert = TRUE, value = TRUE))]
accumulated_sequences <- accumulated_sequences[order(rowSums(accumulated_sequences[,grep("ID|Sequence|Gene_ID|type",colnames(accumulated_sequences), invert = TRUE, value = TRUE)]), decreasing = TRUE),]

# Optionally, write the comparison table to a CSV file
#write.csv(accumulated_sequences, "Contaminating_Sequences.csv", row.names = FALSE)

# Create a heat map visualization by sequence

# Transform data to long format for plotting and use top 20 contaminant sequences
heatmap_data <- accumulated_sequences[1:10,] %>%
  gather(key = "Sample", value = "Count", -Sequence, -Gene_ID, -type)
heatmap_data$Gene_ID <- factor(accumulated_sequences$Gene_ID[1:10], levels = unique(accumulated_sequences$Gene_ID[1:20]), ordered = TRUE)
# Add sequence length information
heatmap_data$Sequence <- paste0("...", heatmap_data$Sequence, "...", " (", nchar(heatmap_data$Sequence), ")")
# Calculate percentages to add to plot
heatmap_data$perc[heatmap_data$Sequence == heatmap_data$Sequence[1]] <- paste0(round(colSums(accumulated_sequences[1:5, 4:length(colnames(accumulated_sequences))])), "%")

ggplot(heatmap_data, aes(x = Sample, y = Sequence, fill = Count)) +
  geom_tile() +
  facet_grid(Gene_ID~.,scales="free_y", space="free_y") +
  theme(strip.text.y.right = element_text(angle = 0),
        axis.text=element_text(size=8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(margin = margin(b = 20))) +
  scale_fill_viridis(name="Percent", option = "inferno") +
  geom_text(label = heatmap_data$perc, vjust = -2, size = 3.5) +
  labs(title = "Read count distribution of top contaminants among samples", x = "Sample", y = "Fragment")+
  coord_cartesian(clip = "off")

ggsave("X:/Omics/RiboSeq/Growth_Conditions_rRNA/Contaminant_Heatmap_byShortestFuzzy_Trimmed.pdf", width = 30, height = 15, units = "cm")

