#############################################################
#
# Ref to the ARTICLE
# 
# Davide Bulgarelli
# 
# 
# d.bulgarelli@dundee.ac.uk
# 
# Revison March 2019
# 
# script to reproduce calculations associated to Controls in Additional files 2
# 
# Disclaimer: the manuscript is currently submitted, the script might be subjected to changes derived from the revision process
#
#############################################################

#############################################################
# Libraries and functions required
#############################################################

rm(list=ls())

#set working directory (from DB's cpu)
setwd("C:/Users/DB42008/Box Sync/Davide_lab_manuscripts/Rodrigo_NT_2018/R_data_no16S/")

#############################################################
#Additional file 1 Table S3
#############################################################

#Import the classification using 16S rRNA genes
datafile_phylum_16S <- read.delim("Additional_file2_ws3.txt", skip=1, sep = "\t", header=TRUE, row.names=1 )
datafile_phylum_16S

#asses the normality of the distributions
shapiro.test(datafile_phylum_16S$Soil.0)
shapiro.test(datafile_phylum_16S$Soil.25)
shapiro.test(datafile_phylum_16S$Soil.100)
shapiro.test(datafile_phylum_16S$NA.)
#Stats
cor.test(datafile_phylum_16S$Soil.0, datafile_phylum_16S$NA.,
          method = "spearman")
cor.test(datafile_phylum_16S$Soil.25, datafile_phylum_16S$NA.,
         method = "spearman")
cor.test(datafile_phylum_16S$Soil.100, datafile_phylum_16S$NA.,
         method = "spearman")

#####################################
#Figure 4
####################################
#Import the list of OTUs enriched in the rhizosphere microhabitat
#############
#N0%
############
Desert.0 <- read.delim("Additional_file2_ws4.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
North.0 <- read.delim("Additional_file2_ws5.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Modern.0 <- read.delim("Additional_file2_ws6.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
#Subset for adjusted p value <0.05
Desert.0_FDR005 <-  Desert.0[(rownames(Desert.0)[which(Desert.0$padj <0.05)]), ] 
North.0_FDR005 <-  North.0[(rownames(North.0)[which(North.0$padj <0.05)]), ] 
Modern.0_FDR005 <-  Modern.0[(rownames(Modern.0)[which(Modern.0$padj <0.05)]), ] 
#Subset for cultivar enriched
Desert.0_FDR005_enriched <- Desert.0_FDR005[(rownames(Desert.0_FDR005)[which(Desert.0_FDR005$log2FoldChange < 0)]), ]
North.0_FDR005_enriched <- North.0_FDR005[(rownames(North.0_FDR005)[which(North.0_FDR005$log2FoldChange < 0)]), ]
Modern.0_FDR005_enriched <- Modern.0_FDR005[(rownames(Modern.0_FDR005)[which(Modern.0_FDR005$log2FoldChange < 0)]), ]
#Inspect the file
dim(Desert.0_FDR005_enriched)
dim(North.0_FDR005_enriched)
dim(Modern.0_FDR005_enriched)
#remove contaminant OTUs (if present)
#import the list of contaminant OTUs DB lab from PMID:30116224
contaminants <- read.delim("contaminant_OTUs_30116224.txt", header = T, row.names = 1)
rownames(contaminants)
length(setdiff(rownames(Desert.0_FDR005_enriched), rownames(contaminants)))
#656
length(setdiff(rownames(North.0_FDR005_enriched), rownames(contaminants)))
#386
length(setdiff(rownames(Modern.0_FDR005_enriched), rownames(contaminants)))
#716

#############
#N25%
############
Desert.25 <- read.delim("Additional_file2_ws7.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
North.25 <- read.delim("Additional_file2_ws8.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Modern.25 <- read.delim("Additional_file2_ws9.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
#Subset for adjusted p value <0.05
Desert.25_FDR005 <-  Desert.25[(rownames(Desert.25)[which(Desert.25$padj <0.05)]), ] 
North.25_FDR005 <-  North.25[(rownames(North.25)[which(North.25$padj <0.05)]), ] 
Modern.25_FDR005 <-  Modern.25[(rownames(Modern.25)[which(Modern.25$padj <0.05)]), ] 
#Subset for cultivar enriched
Desert.25_FDR005_enriched <- Desert.25_FDR005[(rownames(Desert.25_FDR005)[which(Desert.25_FDR005$log2FoldChange < 0)]), ]
North.25_FDR005_enriched <- North.25_FDR005[(rownames(North.25_FDR005)[which(North.25_FDR005$log2FoldChange < 0)]), ]
Modern.25_FDR005_enriched <- Modern.25_FDR005[(rownames(Modern.25_FDR005)[which(Modern.25_FDR005$log2FoldChange < 0)]), ]
#Inspect the file
dim(Desert.25_FDR005_enriched)
dim(North.25_FDR005_enriched)
dim(Modern.25_FDR005_enriched)
#remove contaminant OTUs (if present)
length(setdiff(rownames(Desert.25_FDR005_enriched), rownames(contaminants)))
#499
length(setdiff(rownames(North.25_FDR005_enriched), rownames(contaminants)))
#372
length(setdiff(rownames(Modern.25_FDR005_enriched), rownames(contaminants)))
#501

#############
#N100
############
Desert.100 <- read.delim("Additional_file2_ws10.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
North.100 <- read.delim("Additional_file2_ws11.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Modern.100 <- read.delim("Additional_file2_ws12.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
#Subset for adjusted p value <0.05
Desert.100_FDR005 <-  Desert.100[(rownames(Desert.100)[which(Desert.100$padj <0.05)]), ] 
North.100_FDR005 <-  North.100[(rownames(North.100)[which(North.100$padj <0.05)]), ] 
Modern.100_FDR005 <-  Modern.100[(rownames(Modern.100)[which(Modern.100$padj <0.05)]), ] 
#Subset for cultivar enriched
Desert.100_FDR005_enriched <- Desert.100_FDR005[(rownames(Desert.100_FDR005)[which(Desert.100_FDR005$log2FoldChange < 0)]), ]
North.100_FDR005_enriched <- North.100_FDR005[(rownames(North.100_FDR005)[which(North.100_FDR005$log2FoldChange < 0)]), ]
Modern.100_FDR005_enriched <- Modern.100_FDR005[(rownames(Modern.100_FDR005)[which(Modern.100_FDR005$log2FoldChange < 0)]), ]
#Inspect the file
dim(Desert.100_FDR005_enriched)
dim(North.100_FDR005_enriched)
dim(Modern.100_FDR005_enriched)
#remove contaminant OTUs (if present)
length(setdiff(rownames(Desert.100_FDR005_enriched), rownames(contaminants)))
#467
length(setdiff(rownames(North.100_FDR005_enriched), rownames(contaminants)))
#354
length(setdiff(rownames(Modern.100_FDR005_enriched), rownames(contaminants)))
#398

#Import the list of OTUs enriched in the genotype-effect comparison
#############
#N0%-Desert
############
Desert_Modern.0 <- read.delim("D0_M0_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Modern_Desert.0 <- read.delim("M0_D0_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Desert_North.0 <- read.delim("D0_N0_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
North_Desert.0 <- read.delim("N0_D0_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
##remove contaminant OTUs (if present)
length(setdiff(rownames(Desert_Modern.0), rownames(contaminants)))
#52
length(setdiff(rownames(Modern_Desert.0), rownames(contaminants)))
#43
length(setdiff(rownames(Desert_North.0), rownames(contaminants)))
#124
length(setdiff(rownames(North_Desert.0), rownames(contaminants)))
#8

#############
#N0%-North/Modern
############
North_Modern.0 <- read.delim("N0_M0_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Modern_North.0 <- read.delim("M0_N0_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
##remove contaminant OTUs (if present)
length(setdiff(rownames(North_Modern.0), rownames(contaminants)))
#9
length(setdiff(rownames(Modern_North.0), rownames(contaminants)))
#63

#########################################
#Compile taxonomy information for Additional file 1 Figure S5
#write.table(Desert_Modern.0[setdiff(rownames(Desert_Modern.0), rownames(contaminants)), ], file ="D0_M0_taxonomy_noCM.txt", sep = "\t")
#write.table(Modern_Desert.0[setdiff(rownames(Modern_Desert.0), rownames(contaminants)), ], file ="M0_D0_taxonomy_noCM.txt", sep = "\t")
#write.table(Desert_North.0[setdiff(rownames(Desert.North.0), rownames(contaminants)), ], file ="D0_N0_taxonomy_noCM.txt", sep = "\t")
#write.table(North_Desert.0[setdiff(rownames(North.Desert.0), rownames(contaminants)), ], file ="N0_D0_taxonomy_noCM.txt", sep = "\t")
#write.table(North_Modern.0[setdiff(rownames(North_Modern.0), rownames(contaminants)), ], file ="N0_M0_taxonomy_noCM.txt", sep = "\t")
#write.table(Modern_North.0[setdiff(rownames(Modern_North.0), rownames(contaminants)), ], file ="M0_N0_taxonomy_noCM.txt", sep = "\t")

#############
#N25%-Desert
############
Desert_Modern.25 <- read.delim("D25_M25_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Modern_Desert.25 <- read.delim("M25_D25_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Desert.North.25 <- read.delim("D25_N25_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
##remove contaminant OTUs (if present)
length(setdiff(rownames(Desert_Modern.25), rownames(contaminants)))
#9
length(setdiff(rownames(Modern_Desert.25), rownames(contaminants)))
#12
length(setdiff(rownames(Desert.North.25), rownames(contaminants)))
#2
#############
#N25%-North/Modern
############
North_Modern.25 <- read.delim("N25_M25_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Modern_North.25 <- read.delim("M25_N25_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
##remove contaminant OTUs (if present)
length(setdiff(rownames(North_Modern.25), rownames(contaminants)))
#12
length(setdiff(rownames(Modern_North.25), rownames(contaminants)))
#16

#############
#N100%-Desert
############
Desert_Modern.100 <- read.delim("D100_M100_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Desert_North.100 <- read.delim("D100_N100_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
North_Desert.100 <- read.delim("N100_D100_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
##remove contaminant OTUs (if present)
length(setdiff(rownames(Desert_Modern.100), rownames(contaminants)))
#6
length(setdiff(rownames(Desert_North.100), rownames(contaminants)))
#5
length(setdiff(rownames(North_Desert.100), rownames(contaminants)))
#1
#############
#N100%-North/Modern
############
North_Modern.100 <- read.delim("N100_M100_FDR005.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
##remove contaminant OTUs (if present)
length(setdiff(rownames(North_Modern.100), rownames(contaminants)))
#2
