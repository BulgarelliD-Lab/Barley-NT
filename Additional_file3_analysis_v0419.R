#############################################################
#
# Ref to the ARTICLE
# 
# Davide Bulgarelli
# 
# 
# d.bulgarelli@dundee.ac.uk
# 
# Revison April 2019
# 
# script to reproduce calculations and figures associated to Additional files 3
# 
# Disclaimer: the manuscript is currently submitted, the script might be subjected to changes derived from the revision process
#
#############################################################

#############################################################
# Libraries and functions required
#############################################################


#loand the package(s) needed for the analysis
#source("http://bioconductor.org/biocLite.R")
#biocLite("gplots")
#biocLite("marray")
#biocLite("ggpubr")
library(phyloseq)
library(vegan)
library(gplots)
library(marray)
library(ggpubr)

rm(list=ls())

#set working directory (from DB's cpu)
setwd("C:/Users/DB42008/Box Sync/Davide_lab_manuscripts/Rodrigo_NT_2018/R_data_no16S/")

#print the session info for reproducibility of the code
sessionInfo()

#############################################################

######################
#Additional file 3
######################

#####################
#Figure 5
#####################

#taxonomic characterisation: correlation at phylum level
datafile_phylum <- read.delim("Additional_file3_ws2.txt", sep = "\t", header=TRUE, row.names=1 )
datafile_phylum 

#data normalisation
datafile_phylum_RA <- t(t(datafile_phylum)/colSums(datafile_phylum)) * 100 
datafile_phylum_RA 
colSums(datafile_phylum_RA)

# determine the average % of reads for each phylum
datafile_phylum_RA_mean_sorted <- datafile_phylum_RA [(order(-rowSums(datafile_phylum_RA))), ] 

#Identifythe phyla whose average abundance in above 1% and their aggregated relative abundance
Phylum_mean_topRank <- datafile_phylum_RA_mean_sorted[rownames(datafile_phylum_RA_mean_sorted)[which(rowMeans(datafile_phylum_RA_mean_sorted) > 1)], ]
dim(Phylum_mean_topRank)
colSums(Phylum_mean_topRank)

#overall
mean(colSums(Phylum_mean_topRank))
#transform in matrix for plotting purposes
Phylum_mean_topRank <- as.matrix(Phylum_mean_topRank) 
rownames(Phylum_mean_topRank)
colnames(Phylum_mean_topRank)
#re-arrange the order samples
plotting_samples <- c("X2023_R1R2_PE", "X2024_R1R2_PE", "X2025_R1R2_PE","X2006_R1R2_PE", "X2007_R1R2_PE", "X2009_R1R2_PE", "X2011_R1R2_PE", "X2012_R1R2_PE", "X2013_R1R2_PE", "X2000_R1R2_PE", "X2001_R1R2_PE", "X2002_R1R2_PE")
Phylum_mean_topRank <- Phylum_mean_topRank[, plotting_samples]
Phylum_mean_topRank

#Import the classification using 16S rRNA genes
datafile_phylum_16S <- read.delim("Additional_file2_ws3.txt", skip=1, sep = "\t", header=TRUE, row.names=1 )
datafile_phylum_16S

# determine the average % of reads for each phylum
datafile_phylum_16S_mean_sorted <- datafile_phylum_16S [(order(-rowSums(datafile_phylum_16S))), ] 
datafile_phylum_16S_mean_sorted <- datafile_phylum_16S_mean_sorted * 100
colSums(datafile_phylum_16S_mean_sorted)

#Identifythe phyla whose average abundance in above 1% and their aggregated relative abundance
Phylum_mean_topRank_16S <- datafile_phylum_16S_mean_sorted[rownames(datafile_phylum_16S_mean_sorted)[which(rowMeans(datafile_phylum_16S_mean_sorted) > 1)], ]
dim(Phylum_mean_topRank_16S)
colSums(Phylum_mean_topRank_16S)

#overall
mean(colSums(Phylum_mean_topRank_16S))
#transform in matrix for plotting purposes
Phylum_mean_topRank_16S <- as.matrix(Phylum_mean_topRank_16S) 
rownames(Phylum_mean_topRank_16S)
colnames(Phylum_mean_topRank_16S)

#re-arrange the order samples and remove control samples
plotting_samples_16S <- c("Soil.0", "Desert.0", "North.0","Modern.0", "Soil.25", "Desert.25", "North.25","Modern.25", "Soil.100", "Desert.100", "North.100","Modern.100")
datafile_phylum_16S <- Phylum_mean_topRank_16S[, plotting_samples_16S]
datafile_phylum_16S <- as.matrix(datafile_phylum_16S)

#draw correlations plot for dominant phyla in both dataset
#rename the phyla in the 16S dataset
rownames_datafile_phylum_16S <- c("Proteobacteria", "Bacteroidetes", "Acidobacteria", "Actinobacteria", "Verrucomicrobia", "Crenarchaeota", "Gemmatimonadetes", "Planctomycetes", "Chloroflexi", "Firmicutes", "WS3", "Nitrospirae")
rownames(datafile_phylum_16S) <- rownames_datafile_phylum_16S 
rownames_datafile_phylum_16S

#identify the top phylum in common between the datasets and their relative contributions
common_phyla <- intersect(rownames(Phylum_mean_topRank), rownames(datafile_phylum_16S))
common_phyla
#metagenomics
mean(colSums(Phylum_mean_topRank[common_phyla, ]))
#16S rRNA gene
mean(colSums(datafile_phylum_16S[common_phyla, ]))

#Correlation plots 
Common_phyla_Metagenomics <- as.data.frame(rowMeans(Phylum_mean_topRank[common_phyla, ]))
colnames(Common_phyla_Metagenomics) <- c("Metagenomics_reads")
Common_phyla_16S <- as.data.frame(rowMeans(datafile_phylum_16S[common_phyla, ]))
colnames(Common_phyla_16S) <- c("rRNA_reads")

#merge the datasets
Correlation_phyla <- cbind(Common_phyla_Metagenomics, Common_phyla_16S)
Correlation_phyla <- log2(Correlation_phyla)
#plotting
dev.off()
ggscatter(Correlation_phyla, x = "rRNA_reads", y = "Metagenomics_reads", add = "reg.line", conf.int = TRUE,                                  # Add confidence interval
                  add.params = list(color = "blue", fill = "lightgray"), label = rownames(Correlation_phyla))
#asses the normality of the distributions
shapiro.test(Correlation_phyla$rRNA_reads)
shapiro.test(Correlation_phyla$Metagenomics_reads)
#Stats
cor.test( ~ Metagenomics_reads + rRNA_reads,
          data=Correlation_phyla,
          method = "pearson")

#####################################
#Additional file 1: Figure S7
####################################

#import the data frame with relative abundance of seed function level 1
datafile_SEED1 <- read.delim("Additional_file3_ws3.txt", sep = "\t", header=TRUE, row.names=1 )
datafile_SEED1_RA <- t(t(datafile_SEED1)/colSums(datafile_SEED1)) * 1000
colSums(datafile_SEED1_RA)
colnames(datafile_SEED1_RA)
dim(datafile_SEED1_RA)

#import the design file
design_metagenomics <- read.delim("Additional_file3_ws1.txt", sep = "\t", header=TRUE, row.names=1 )
design_metagenomics

#create a phyloseq object to generate an ordination on SEED function
SEED_counts <- otu_table(datafile_SEED1_RA, taxa_are_rows=TRUE)
SEED_design <- sample_data(design_metagenomics)
SEED_phyloseq <- merge_phyloseq(SEED_counts, SEED_design)

#PCoA bray distance
SEED_prop_bray <- ordinate(SEED_phyloseq, "PCoA", "bray")
plot_ordination(SEED_phyloseq, SEED_prop_bray, color = "Genotype")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(SEED_phyloseq, SEED_prop_bray, shape ="Microhabitat", color = "Genotype")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("brown","yellow", "magenta", "green"))
p + ggtitle("PCoA SEED data, Bray distance")

#BC distance adonis
BC <- phyloseq::distance(SEED_phyloseq, "bray")
#Microhabitat effect
adonis(BC ~ Microhabitat, data= design_metagenomics, permutations = 5000)
#Genotype effect
#Subset for rhizosphere samples
SEED_phyloseq_rhizo <- subset_samples(SEED_phyloseq, Microhabitat == "Rhizosphere")
design_metagenomics_rhizo <- design_metagenomics[colnames(otu_table(SEED_phyloseq_rhizo)), ]
BC_rhizo <- phyloseq::distance(SEED_phyloseq_rhizo, "bray")
adonis(BC_rhizo ~ Genotype, data= design_metagenomics_rhizo, permutations = 5000)

#end
