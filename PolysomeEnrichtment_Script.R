#RNA-Seq Pipeline for total RNA and polysome enriched samples of NAA treated and control plants 

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
    colnames(STAR.counts) <- gsub("ReadsPerGene.*$","",list.files("./STAR_Output/", pattern = "ReadsPerGene")[1:length(colnames(STAR.counts))])
  }
  rm(temp,i)
}

# BUILD METADATA ______________________________________________________________________________________

#Generate sample table containing experiment information 
sample.table <- data.frame(sample = colnames(STAR.counts))
rownames(sample.table) <- sample.table$sample
sample.table$treatment <- rep(c("Ctrl", "NAA"), each = 4)
sample.table$condition <- rep(c("Polysomes", "Total RNA"), each = 8)
sample.table$replicate <- gsub("^.*_NAA_", "", sample.table$sample)
sample.table$replicate <- gsub("^.*_Ctrl_", "", sample.table$replicate)
sample.table$replicate <- gsub("\\.$", "", sample.table$replicate)
