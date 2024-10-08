#RNA-Seq Pipeline for total RNA and polysome enriched samples of NAA treated and control plants 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RNA-Seq #####################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# PREPARE RAW DATA ############################################################################################

#Concatenate files from different lanes into single file
#Rename files to a uniform format if not already done (eg: Experiment_Sample_Treatment_Replicate_Mate -> NAA_90min_Ctrl_1_1.fastq.gz)
#setwd("Y:/Omics Data/RiboSeq/PolysomeEnrichment/RawData/cleaned)
#file.rename(dir(), gsub("-MACE_D.*Seq", "",dir()))
#Only work with gzipped files for more efficient storage
#Keep all files in a single folder


# QUALITY CONTROL _____________________________________________________________________________________________

#Run FastQC on the read files
setwd("Y:/Omics Data/RiboSeq/PolysomeEnrichment/RawData/cleaned/") #switch to directory containing the reads
system(paste("wsl fastqc *.fastq.gz --outdir ../FastQC/")) #Run FastQC over all files and save results to specified folder

