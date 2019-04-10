#############################################################
#
# Ref to the ARTICLE
# 
# Rodrigo Alegria Terrazas & Davide Bulgarelli
# r.alegriaterrazas@dundee.ac.uk
# d.bulgarelli@dundee.ac.uk
# 
# Revison March 2019
# 
# script to reproduce calculations and figures presented in the manuscript
# 
# Disclaimer: the manuscript is currently submitted, the script might be subjected to changes derived from the revision process
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#1st time installation
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
#biocLite("DESeq2")
#biocLite("PMCMR")
#biocLite("vegan")

#biocLite("agricolae")#HSD test~ I dont think we use it

#required packages 
library("phyloseq")
library("ggplot2")
library("vegan")
library ("ape")
library("PMCMR")
library("plyr")
library("car")#


#from davides figure 5 and S7 analysis
#biocLite("gplots")
#biocLite("marray")
#biocLite("ggpubr")
library(gplots)
library(marray)
library(ggpubr)



#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()


#set the working directory
setwd("C:/Users/ra42320/Box Sync/Davide_lab_manuscripts/Rodrigo_NT_2018/R_data/Rodrigo")
#############################################################
#import the count matrix and the desing file
#############################################################


#OTU table generated using QIIME 1.9.0. 
#additional file 2, supplementary worksheet 2
dat_info <- read.delim("JH02_NT_otu_table_non_chloro_non_mito.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)


#extract the OTUids from the GreenGenes OTU table to filter silva otu table in QIIME
rownames(dat_info)
#write(rownames(dat_info), "JH02_NT_noPlant_OTUs_id_gg135.txt")

#inspect the file
dim(dat_info)
colnames(dat_info)
rownames(dat_info)

#extract the total number of reads clustered at OTU 97% identiy (the number underneath the Hv identifier represents the total number of reads clustered for that sample) 
#min number of reads Hv965 - 25810
OTU_97_reads <- sort(colSums(dat_info[, 1:50]))
OTU_97_reads

#total reads  5262776
OTU_97_reads_sum <- sum(colSums(dat_info[, 1:50]))
OTU_97_reads_sum

#design file (Additional file 2, supplementary worksheet1)
design <- read.delim("Map_JH02_NT_noCoast_treat.txt", sep = "\t", header=TRUE, row.names=1)
design

#extract count information
dat_count_noPlants <- dat_info[, 1:50]
dim(dat_count_noPlants)
colnames(dat_count_noPlants)

#and a new taxa table
dat_tax_noPlants <- as.data.frame(dat_info[rownames(dat_count_noPlants), 51])
rownames(dat_tax_noPlants) <- rownames(dat_count_noPlants)
dat_tax_noPlants[1:5, ]
#save the above file and in excel we will create a new tax table where each column represents a taxonomic rank
#write.table(dat_tax_noPlants, file="JH02_NT_dat_tax_noPlants.txt", sep="\t")
#file "JH02_NT_dat_tax_noPlants.txt" was used to create JH02_NT_dat_tax_noPlant_ordered for downstream analyses

#############################################################
#Figure 3: High taxonomic ranks distribution
#Data required: Additional_file2_ws3.txt
#############################################################

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

# Stacked Bar Plot with Colors and Legend
barplot(datafile_phylum_16S, main="Phylum Distribution, 16S rRNA gene dataset",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkblue","red", "green", "cyan", "orange", "brown", "black", "magenta", "grey", "white", "Purple", "Pink"), beside=FALSE,   legend = rownames(datafile_phylum_16S))

#Due to size limits the legend covers part of the graph, save the graph as .eps file and in illustrator uou can adjust this (and other)  graphical issues
#barplot_no legend: you can save this as image and then reconstruct the leged using the previous figure. Note: the order of samples can be inferred from the command on line 256
barplot(datafile_phylum_16S, main="Phylum Distribution, 16S rRNA gene dataset",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkblue","red", "green", "cyan", "orange", "brown", "black", "magenta", "grey", "white", "Purple", "Pink"), beside=FALSE)



##########################################################################

#############################################################
#Genererate the phyloseq object
#Data required: dat_count_noPlants; design, JH02_NT_dat_tax_noPlants_ordered.txt, and 97_otus.tree.gz
#############################################################

#The OTU Table counts
JH02_NT_OTU <- otu_table(dat_count_noPlants, taxa_are_rows=TRUE)

#The taxonomy information
#Note that the file JH02_NT_dat_tax_noPlants_ordered.txt has been generated from the output of lines  
#it is a tab-delimited file with 8 columns, the column headers are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
JH02_NT_taxa_ordered <- read.delim ("JH02_NT_dat_tax_noPlants_ordered.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = FALSE)
JH02_NT_taxa <- tax_table(as.matrix(JH02_NT_taxa_ordered))
dim(JH02_NT_taxa)

#The mapping file 
JH02_NT_map <- sample_data(design)
JH02_NT_map
#The phylogenetic tree: the OTU table has been generated using a closed reference approach agains the greengenes 13_05 database use the corresponding phylogenetic tree 
JH02_NT_tree <- read_tree_greengenes("97_otus.tree.gz")

#check whether the tree is rooted
is.rooted(JH02_NT_tree)

#merge the files and create the phyloseq object
JH02_NT_data_phyloseq <- merge_phyloseq(JH02_NT_OTU, JH02_NT_taxa, JH02_NT_map,  JH02_NT_tree)

#inspect the generated data 8874 OTUs, 5262776 reads
JH02_NT_data_phyloseq
sum(colSums(otu_table(JH02_NT_data_phyloseq)))
dim(dat_count_noPlants)
sum(colSums(dat_count_noPlants))


##########################################################################

#remove the control samples from the dataset

JH02_NT_data_phyloseq_2<- subset_samples(JH02_NT_data_phyloseq, Microhabitat=="Rhizosphere")
JH02_NT_data_phyloseq_3<- subset_samples(JH02_NT_data_phyloseq, Microhabitat=="Bulk")
subset_combined<- merge_phyloseq(JH02_NT_data_phyloseq_2,JH02_NT_data_phyloseq_3)
JH02_NT_data_phyloseq_4 <-subset_combined

design_2 <- design[colnames(otu_table(JH02_NT_data_phyloseq_4)), ]

##########################################################################
#Dry weight analysis
##########################################################################

#Data visualisation#re-arrange the order for plotting purposes

rhizo_samples <- rownames(design_2)[which(design_2$Microhabitat!= "Bulk")]
design_2_rhizo <- design_2[rhizo_samples, ]

design_2_rhizo_plot <-as.data.frame(design_2_rhizo)

#re-order the factors
design_2_rhizo_plot$Description <- ordered(design_2_rhizo$Description, levels=c("North.0","North.25", "North.100","Desert.0","Desert.25","Desert.100","Modern.0", "Modern.25", "Modern.100"))

p <- ggplot(design_2_rhizo, aes(x=Description, y=Dryweight, fill=B1K)) + geom_boxplot ()
p + scale_fill_manual(values = c("B1K.31.01" = "green","B1K.12.10" ="yellow","Morex" ="magenta"))


#test the normality of the observed data (data not normally distributed!) 
shapiro.test(design_2_rhizo$Dryweight)

qqnorm(design_2_rhizo$Dryweight)
qqline(design_2_rhizo$Dryweight, col = 2)

kruskal.test(Dryweight ~ B1K, data =design_2_rhizo )
kruskal.test(Dryweight ~ Description, data =design_2_rhizo )
posthoc.kruskal.dunn.test (x= design_2_rhizo$Dryweight, g= design_2_rhizo$Description, p.adjust.method="BH")


#################################################################################################################
#Plant N content

design_2_rhizo$Description <- ordered(design_2_rhizo$Description, levels=c("North.0","North.25", "North.100","Desert.0","Desert.25","Desert.100","Modern.0", "Modern.25", "Modern.100"))

p <- ggplot(design_2_rhizo, aes(x=Description, y=N_content, fill=B1K)) + geom_boxplot ()
p + scale_fill_manual(values = c("B1K.31.01" = "green","B1K.12.10" ="yellow","Morex" ="magenta"))


#test the normality of the observed data
######data not normally distributed! 

shapiro.test(design_2_rhizo$N_content)

qqnorm(design_2_rhizo$N_content)
qqline(design_2_rhizo$N_content, col = 2)

kruskal.test(N_content ~ Treatment, data =design_2_rhizo)
kruskal.test(N_content ~ Description, data =design_2_rhizo)
posthoc.kruskal.dunn.test (x= design_2_rhizo$N_content, g= design_2_rhizo$Description, p.adjust.method="BH")


####################################################################################################################
#NO3 & NH4


design_2$Description <- ordered(design_2$Description, levels=c("Soil.0", "Soil.25", "Soil.100", "North.0","North.25", "North.100","Desert.0","Desert.25","Desert.100","Modern.0", "Modern.25", "Modern.100"))

p <- ggplot(design_2, aes(x=Description, y=NO3, fill=B1K)) + geom_boxplot ()
p+scale_fill_manual(values=c("Bulk"="brown","B1K.31.01" = "green","B1K.12.10" ="yellow","Morex" ="magenta" ))


p <- ggplot(design_2, aes(x=Description, y=NH4, fill=B1K)) + geom_boxplot ()
p+scale_fill_manual(values=c("Bulk"="brown","B1K.31.01" = "green","B1K.12.10" ="yellow","Morex" ="magenta" ))


#normality test
##data not normally distributed
#
shapiro.test(design_2$NO3)
shapiro.test(design_2$NH4)

qqnorm(design_2$NO3)
qqline(design_2$NO3, col = 2)

qqnorm(design_2$NH4)
qqline(design_2$NH4, col = 2)

#NO3stats


#not significant for NO3 and NH4
wilcox.test(NO3 ~ Microhabitat, data = design_2 )
wilcox.test(NH4 ~ Microhabitat, data = design_2 )

kruskal.test(NO3 ~ B1K, data =design_2)#Not significant
kruskal.test(NO3 ~ Treatment, data =design_2)#significant
kruskal.test(NO3 ~ Description, data =design_2)#significant
posthoc.kruskal.dunn.test (x= design_2$NO3, g= design_2$Description, p.adjust.method="BH")


kruskal.test(NH4 ~ B1K, data =design_2)#Not significant
kruskal.test(NH4 ~ Treatment, data =design_2)#significant
kruskal.test(NH4 ~ Description, data =design_2)#significant
posthoc.kruskal.dunn.test (x= design_2$NH4, g= design_2$Description, p.adjust.method="BH")


##################################################################################################################
#Figure 4: Alphadiversity calculations
#Data required: design_2; JH02_JH03_RHM_data_phyloseq_rare_table_counts_2.txt
##################################################################################################################

####################Remind to do abundance filtering first######################################################

#rarefy the dataset 1238880 reads, 7639 otus, 25810 reads/sample
#JH02_NT_data_phyloseq_rare <- rarefy_even_depth(JH02_NT_data_phyloseq_4, rngseed=TRUE)
#sum(colSums(otu_table(JH02_NT_data_phyloseq_rare)))

#extract and save the OTU table for reproducibility of the code
#JH02_NT_data_phyloseq_rare_table <- as.data.frame(otu_table(JH02_NT_data_phyloseq_rare))

#inspect the generated file
#class(JH02_NT_data_phyloseq_rare_table)
#dim(JH02_NT_data_phyloseq_rare_table)

#save the file for the reproducibility of the code
#write.table(JH02_NT_data_phyloseq_rare_table, file="JH02_NT_data_phyloseq_rare_table.txt", sep="\t")


#import the rarefied OTU counts (note file name counts2.txt??? this file has been generated in excel and includes the #OTU ID as header of the first column)
dat_count_rare <- read.delim("JH02_NT_data_phyloseq_rare_table.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#inspect the generated file
dim(dat_count_rare)
colSums(dat_count_rare)

#generate a new phyloseq object wich will contain only the rarefied counts and the design file (only these two pieces of information are required for alphadiversity calculation)
JH02_NT_OTU_rare <- otu_table(dat_count_rare, taxa_are_rows=TRUE)
JH02_NT_map_rare <- sample_data(design_2)
JH02_NT_data_rare_phyloseq <- merge_phyloseq(JH02_NT_OTU_rare, JH02_NT_map_rare)

#Inspect the generated file
JH02_NT_data_rare_phyloseq
sample_sums(JH02_NT_data_rare_phyloseq)

#Index calculations
JH02_NT_alpha_rare <-  estimate_richness(JH02_NT_data_rare_phyloseq, measures = c("Observed", "Shannon", "chao1"))

#generate a new dataframes for data visualisation

#Sample information

#Genotype
design_genotype <- as.data.frame(design_2[, 9])
rownames(design_genotype) <- rownames(design_2)
colnames(design_genotype) <- c("Genotype")

#Sample
design_description <- as.data.frame(design_2[, 19])
rownames(design_description) <- rownames(design_2)
colnames(design_description) <- c("Description")

#data frame Genotype_Description
design_GD <- cbind(design_genotype, design_description)

#Observed OTUs
JH02_NT_alpha_rare_Observed <- as.data.frame(JH02_NT_alpha_rare[ ,1])
rownames(JH02_NT_alpha_rare_Observed) <- rownames(JH02_NT_alpha_rare)
colnames(JH02_NT_alpha_rare_Observed) <- c("Observed")

#Combine the dataset sample description and Observed OTUs
JH02_NT_alpha_rare_Observed_D <- cbind(design_GD, JH02_NT_alpha_rare_Observed)

#Order the levels according to a defined order
#JH02_NT_alpha_rare_Observed_D$Description <- ordered(JH02_NT_alpha_rare_Observed_D$Description , levels=c("Soil.0","Soil.25","Soil.100","Desert.0","Desert.25","Desert.100","North.0","North.25", "North.100","Modern.0", "Modern.25", "Modern.100"))
JH02_NT_alpha_rare_Observed_D$Description <- ordered(JH02_NT_alpha_rare_Observed_D$Description , levels=c("Soil.0","North.0","Desert.0","Modern.0", "Soil.25","North.25","Desert.25","Modern.25","Soil.100","North.100","Desert.100", "Modern.100"))


#plotting
p <- ggplot(JH02_NT_alpha_rare_Observed_D, aes(x=Description, y=Observed, fill=Genotype)) + geom_boxplot()
p+scale_fill_manual(values=c("yellow", "green", "brown", "magenta"))

#Shannon
JH02_NT_alpha_rare_Shannon <- as.data.frame(JH02_NT_alpha_rare[ ,4])
rownames(JH02_NT_alpha_rare_Shannon) <- rownames(JH02_NT_alpha_rare)
colnames(JH02_NT_alpha_rare_Shannon) <- c("Shannon")

#Combine the dataset sample description and Shannon OTUs
JH02_NT_alpha_rare_Shannon_D <- cbind(design_GD, JH02_NT_alpha_rare_Shannon)

#Order the levels according to a defined order
#JH02_NT_alpha_rare_Shannon_D$Description <- ordered(JH02_NT_alpha_rare_Shannon_D$Description, levels=c("Soil.0","Soil.25","Soil.100","Desert.0","Desert.25","Desert.100","North.0","North.25", "North.100","Modern.0", "Modern.25", "Modern.100"))
JH02_NT_alpha_rare_Shannon_D$Description <- ordered(JH02_NT_alpha_rare_Shannon_D$Description, levels=c("Soil.0","North.0","Desert.0","Modern.0", "Soil.25","North.25","Desert.25","Modern.25","Soil.100","North.100","Desert.100", "Modern.100"))

#plotting
p <- ggplot(JH02_NT_alpha_rare_Shannon_D, aes(x=Description, y=Shannon, fill= Genotype)) + geom_boxplot()
p + scale_fill_manual(values = c("yellow", "green", "brown", "magenta"))

#Chao1
JH02_NT_alpha_rare_Chao1 <- as.data.frame(JH02_NT_alpha_rare[ ,2])
rownames(JH02_NT_alpha_rare_Chao1) <- rownames(JH02_NT_alpha_rare)
colnames(JH02_NT_alpha_rare_Chao1) <- c("Chao1")

#Combine the dataset sample description and Chao1 OTUs
JH02_NT_alpha_rare_Chao1_D <- cbind(design_GD, JH02_NT_alpha_rare_Chao1)

#Order the levels according to a defined order
#JH02_NT_alpha_rare_Chao1_D$Description <- ordered(JH02_NT_alpha_rare_Chao1_D$Description, levels=c("Soil.0","Soil.25","Soil.100","Desert.0","Desert.25","Desert.100","North.0","North.25", "North.100","Modern.0", "Modern.25", "Modern.100"))
JH02_NT_alpha_rare_Chao1_D$Description <- ordered(JH02_NT_alpha_rare_Chao1_D$Description, levels=c("Soil.0","North.0","Desert.0","Modern.0", "Soil.25","North.25","Desert.25","Modern.25","Soil.100","North.100","Desert.100", "Modern.100"))


#plotting
p <- ggplot(JH02_NT_alpha_rare_Chao1_D, aes(x=Description, y=Chao1, fill= Genotype)) + geom_boxplot()
p + scale_fill_manual(values = c("yellow", "green", "brown", "magenta"))


#generate a new dataframe for statistical analysis
JH02_NT_alpha_rare_info <- cbind(design_2, JH02_NT_alpha_rare)

#check the new dataset: it contains both the description of the samples and alpha 
JH02_NT_alpha_rare_info


#Test the soil effect on alphadiversity indices
#Soil effect
#generate a dataframe for statistical analysis

design_microhabitat <- as.data.frame(design_2[, 7])
rownames(design_microhabitat) <- rownames(design_2)
colnames(design_microhabitat) <- c("Microhabitat")

JH02_NT_alpha_rare_info_micro <- cbind(design_microhabitat, JH02_NT_alpha_rare)

#check normality
#(data normally distributed)
shapiro.test(JH02_NT_alpha_rare_info_micro$Observed)

#Observed- significant
wilcox.test(Observed ~ Microhabitat, data = JH02_NT_alpha_rare_info_micro)
#Shannon-significant
wilcox.test(Shannon ~ Microhabitat, data = JH02_NT_alpha_rare_info_micro)
#Chao1-significant
wilcox.test(Chao1 ~ Microhabitat, data = JH02_NT_alpha_rare_info_micro)

t.test(Chao1 ~ Microhabitat, data = JH02_NT_alpha_rare_info_micro)
t.test(Observed ~ Microhabitat, data = JH02_NT_alpha_rare_info_micro)
t.test(Shannon ~ Microhabitat, data = JH02_NT_alpha_rare_info_micro)

#Genotype not corrected for soil (use the parameter 'Description' in the dataset)
#Observed (not significant)
kruskal.test(Observed ~ Description, data = JH02_NT_alpha_rare_info)
posthoc.kruskal.dunn.test (x=JH02_NT_alpha_rare_info$Observed, g=JH02_NT_alpha_rare_info$Description, p.adjust.method="BH")

#Shannon (significant!)
kruskal.test(Shannon ~ Description, data = JH02_NT_alpha_rare_info)
posthoc.kruskal.dunn.test (x=JH02_NT_alpha_rare_info$Shannon, g=JH02_NT_alpha_rare_info$Description, p.adjust.method="BH")

#Chao1 (not significant)
kruskal.test(Chao1 ~ Description, data = JH02_NT_alpha_rare_info)
posthoc.kruskal.dunn.test (x=JH02_NT_alpha_rare_info$Chao1, g=JH02_NT_alpha_rare_info$Description, p.adjust.method="BH")

###################removing soil samples 

JH02_NT_alpha_rare_info_rhizo<- subset(JH02_NT_alpha_rare_info, Microhabitat=="Rhizosphere")

#Observed (not significant)
kruskal.test(Observed ~ Description, data = JH02_NT_alpha_rare_info_rhizo)
posthoc.kruskal.dunn.test (x=JH02_NT_alpha_rare_info_rhizo$Observed, g=JH02_NT_alpha_rare_info_rhizo$Description, p.adjust.method="BH")

#Shannon (not significant)
kruskal.test(Shannon ~ Description, data = JH02_NT_alpha_rare_info_rhizo)
posthoc.kruskal.dunn.test (x=JH02_NT_alpha_rare_info_rhizo$Shannon, g=JH02_NT_alpha_rare_info_rhizo$Description, p.adjust.method="BH")

#Chao1 (not significant)
kruskal.test(Chao1 ~ Description, data = JH02_NT_alpha_rare_info_rhizo)
posthoc.kruskal.dunn.test (x=JH02_NT_alpha_rare_info_rhizo$Chao1, g=JH02_NT_alpha_rare_info_rhizo$Description, p.adjust.method="BH")

#NH4stats
res.aov_obs <- aov(Shannon ~ B1K*Treatment, data =JH02_NT_alpha_rare_info_rhizo  )
summary(res.aov_obs)
HSD.test(res.aov_obs,"B1K", group=TRUE,console=TRUE)
HSD.test(res.aov_obs,"Treatment", group=TRUE,console=TRUE)

#generate model using the description for post hoc test
res.aov_obs_2 <- aov(Shannon ~ Description, data = JH02_NT_alpha_rare_info_rhizo)
summary(res.aov_obs_2)
TukeyHSD(res.aov_obs_2, which = "Description")
HSD.test(res.aov_obs_2,"Description", group=TRUE,console=TRUE)


#############################################################NOn RAREFIED###########
#Figure 5: Betadiversity calculations
#Data required: design_2; JH02_NT_data_phyloseq_4
#############################################################

#abundance filtering
#Remove OTUs not seen more than 5 times in at least 20% of the samples
JH02_NT_data_phyloseq_5 = filter_taxa(JH02_NT_data_phyloseq_4, function(x) sum(x > 5) > (0.2*length(x)), TRUE)
#8874 taxa 5087783 reads
JH02_NT_data_phyloseq_4
#2628 taxa 4868770 reads
JH02_NT_data_phyloseq_5
sum(sample_sums(JH02_NT_data_phyloseq_4))
sum(sample_sums(JH02_NT_data_phyloseq_5))

#ratio filtered reads/total reads 95.69532
ratio <- sum(sample_sums(JH02_NT_data_phyloseq_5))/sum(sample_sums(JH02_NT_data_phyloseq_4))*100


#Transform the count in relative abundance cpm
JH02_NT_data_phyloseq_prop <- transform_sample_counts(JH02_NT_data_phyloseq_5,  function(x) 1e+06 * x/sum(x))
sample_sums(JH02_NT_data_phyloseq_prop) 

#PCoA weighted unifrac distance
JH02_NT_data_phyloseq_prop_wunifrac <- ordinate(JH02_NT_data_phyloseq_prop, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(JH02_NT_data_phyloseq_prop, JH02_NT_data_phyloseq_prop_wunifrac , color = "Treatment")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH02_NT_data_phyloseq_prop, JH02_NT_data_phyloseq_prop_wunifrac , shape ="B1K", color = "Treatment")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("red","green","blue"))
p + ggtitle("PCoA 16S data, weighted Unifrac distance")

#PCoA bray distance
JH02_NT_data_phyloseq_prop_bray <- ordinate(JH02_NT_data_phyloseq_prop, "PCoA", "bray")
plot_ordination(JH02_NT_data_phyloseq_prop, JH02_NT_data_phyloseq_prop_bray , color = "Treatment")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH02_NT_data_phyloseq_prop, JH02_NT_data_phyloseq_prop_bray , shape ="B1K", color = "Treatment")
p = p + geom_point(size = 6, alpha = 0.75)
p = p + scale_colour_manual(values = c("red","green","blue"))
p + ggtitle("PCoA 16S data, Bray distance")

#adonis calculations

#rhizosphere effect

#WU distance
WU <- phyloseq::distance(JH02_NT_data_phyloseq_prop, "unifrac", weighted= TRUE)
adonis(WU ~ B1K * Microhabitat, data= design_2, permutations = 5000)

#BC distance
BC <- phyloseq::distance(JH02_NT_data_phyloseq_prop, "bray")
adonis(BC ~ B1K * Microhabitat, data= design_2, permutations = 5000)

#genotype effect (rhizosphere only)
#Subsetting
JH02_NT_data_phyloseq_prop_rhizo <- subset_samples(JH02_NT_data_phyloseq_prop, Microhabitat == "Rhizosphere")
design_rhizosphere <- design_2[colnames(otu_table(JH02_NT_data_phyloseq_prop_rhizo)), ]

#WU distance
WU <- phyloseq::distance(JH02_NT_data_phyloseq_prop_rhizo, "unifrac", weighted= TRUE)
adonis(WU ~ B1K * Treatment, data= design_rhizosphere, permutations = 5000)

#BC distance
BC <- phyloseq::distance(JH02_NT_data_phyloseq_prop_rhizo, "bray")
adonis(BC ~ B1K * Treatment, data= design_rhizosphere, permutations = 5000)

#treatment effect
#WU distance
WU <- phyloseq::distance(JH02_NT_data_phyloseq_prop, "unifrac", weighted= TRUE)
adonis(WU ~ B1K * Treatment, data= design_2, permutations = 5000)

#BC distance
BC <- phyloseq::distance(JH02_NT_data_phyloseq_prop, "bray")
adonis(BC ~ B1K * Treatment, data= design_2, permutations = 5000)

######try only soils samples

# treatment effect (soil only)
#Subsetting
JH02_NT_data_phyloseq_prop_bulk <- subset_samples(JH02_NT_data_phyloseq_prop, Microhabitat == "Bulk")
design_bulk <- design_2[colnames(otu_table(JH02_NT_data_phyloseq_prop_bulk)), ]

#WU distance
WU <- phyloseq::distance(JH02_NT_data_phyloseq_prop_bulk, "unifrac", weighted= TRUE)
adonis(WU ~ Treatment, data= design_bulk, permutations = 5000)

#BC distance
BC <- phyloseq::distance(JH02_NT_data_phyloseq_prop_bulk, "bray")
adonis(BC ~ Treatment, data= design_bulk, permutations = 5000)


#####################################################################
#BETADIVERSITY RARIFIED SAMPLES
##Data required: design_2; JH02_NT_data_phyloseq_rare OTU table need to be extracted B4

####################################################################
####create phyloseq object: JH02_NT_data_phyloseq_rare_2#############

#abundance filtering
#Remove OTUs not seen more than 5 times in at least 20% of the samples
JH02_NT_data_phyloseq_rare_2 = filter_taxa(JH02_NT_data_phyloseq_rare, function(x) sum(x > 5) > (0.2*length(x)), TRUE)
#8874 taxa 5087783 reads

JH02_NT_data_phyloseq_rare #16/01/18   7639 taxa/// 1238880 reads
#2628 taxa 4868770 reads
JH02_NT_data_phyloseq_rare_2 #16/01/18 1204 taxa///1064450 reads
sum(sample_sums(JH02_NT_data_phyloseq_rare))
sum(sample_sums(JH02_NT_data_phyloseq_rare_2))

#ratio filtered reads/total reads 95.69532
# #16/01/18  85.92035
ratio <- sum(sample_sums(JH02_NT_data_phyloseq_rare_2))/sum(sample_sums(JH02_NT_data_phyloseq_rare))*100
ratio


#PCoA weighted unifrac distance
JH02_NT_data_phyloseq_rare_wunifrac <- ordinate(JH02_NT_data_phyloseq_rare_2, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(JH02_NT_data_phyloseq_rare_2, JH02_NT_data_phyloseq_rare_wunifrac , color = "Treatment")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH02_NT_data_phyloseq_rare_2, JH02_NT_data_phyloseq_rare_wunifrac , shape ="B1K", color = "Treatment")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("red","green","blue"))
p + ggtitle("PCoA 16S data, weighted Unifrac distance")

#PCoA bray distance
JH02_NT_data_phyloseq_rare_bray <- ordinate(JH02_NT_data_phyloseq_rare_2, "PCoA", "bray")
plot_ordination(JH02_NT_data_phyloseq_rare_2, JH02_NT_data_phyloseq_rare_bray , color = "Treatment")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH02_NT_data_phyloseq_rare_2, JH02_NT_data_phyloseq_rare_bray , shape ="B1K", color = "Treatment")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("red","green","blue"))
p + ggtitle("PCoA 16S data, Bray distance")

#adonis calculations

#rhizosphere effect

#WU distance
WU <- phyloseq::distance(JH02_NT_data_phyloseq_rare_2, "wunifrac")
adonis(WU ~ B1K * Microhabitat, data= design_2, permutations = 5000)

#BC distance
BC <- phyloseq::distance(JH02_NT_data_phyloseq_rare_2, "bray")
adonis(BC ~ B1K * Microhabitat, data= design_2, permutations = 5000)

#genotype effect (rhizosphere only)
#Subsetting
JH02_NT_data_phyloseq_rare_rhizo <- subset_samples(JH02_NT_data_phyloseq_rare_2, Microhabitat == "Rhizosphere")
design_rhizosphere <- design_2[colnames(otu_table(JH02_NT_data_phyloseq_prop_rhizo)), ]

#WU distance########check
WU <- phyloseq::distance(JH02_NT_data_phyloseq_rare_rhizo, "wunifrac")
adonis(WU ~ B1K * Treatment, data= design_rhizosphere, permutations = 5000)

#BC distance
BC <- phyloseq::distance(JH02_NT_data_phyloseq_rare_rhizo, "bray")
adonis(BC ~ B1K * Treatment, data= design_rhizosphere, permutations = 5000)

#treatment effect
#WU distance
WU <- phyloseq::distance(JH02_NT_data_phyloseq_rare_rhizo, "unifrac", weighted= TRUE)
adonis(WU ~ B1K * Treatment, data= design_2, permutations = 5000)

#BC distance
BC <- phyloseq::distance(JH02_NT_data_phyloseq_rare_rhizo, "bray")
adonis(BC ~ B1K * Treatment, data= design_2, permutations = 5000)

######try only soils samples

# treatment effect (soil only)
#Subsetting
JH02_NT_data_phyloseq_prop_bulk <- subset_samples(JH02_NT_data_phyloseq_prop, Microhabitat == "Bulk")
design_bulk <- design_2[colnames(otu_table(JH02_NT_data_phyloseq_prop_bulk)), ]

#WU distance
WU <- phyloseq::distance(JH02_NT_data_phyloseq_prop_bulk, "unifrac", weighted= TRUE)
adonis(WU ~ Treatment, data= design_bulk, permutations = 5000)

#BC distance
BC <- phyloseq::distance(JH02_NT_data_phyloseq_prop_bulk, "bray")
adonis(BC ~ Treatment, data= design_bulk, permutations = 5000)



######################################################################################################
#Differentially recruited OTUs
#Data required: design_2;JH02_NT_data_phyloseq_5 ; dat_tax_noPlants
##Hypothesis testing: create a new OTU table to be imported in qiime to test OTUs differentially regulated
#we use DEseq https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y
#############################################################################################################

#extract OTU counts (non rarefied, yet filtered for abundance) from an abundance threshold phyloseq object
dat_count_threshold  <- as.data.frame(otu_table(JH02_NT_data_phyloseq_5))

#extract the cognate taxonomy information
dat_tax_threshold  <- as.data.frame(dat_tax_noPlants[rownames(dat_count_threshold  ), ])
colnames(dat_tax_threshold ) <- c("ConsensusLineage")
colnames(dat_tax_threshold)

#merge count and taxonomy files
dat_info_threshold <- cbind(dat_count_threshold , dat_tax_threshold )
dim(dat_info_threshold)
colnames(dat_info_threshold)

#save the files for the reproducibility of the code and the analysis in QIIME
#write.table(dat_info_threshold, file="JH02_NT_info_threshold.txt")
# save mapping file for DESeq calculations
#write.table(design_2, file="Map_JH02_NT_noCoast_noControls.txt") 

#open the file with excel and add the following name to the first column, #OTU ID
#them convert to biom for QIIME

#Subsetting enriched OTUs from DESeq
#Note when you look at the fold change: positive fold change (>0) is something enriched in the first term of comparison

#alternative method
#D0_rhizo_enriched <- read.delim("JH02_NT_D0_rhizosphere.txt")
#D0_enriched_M0 <- read.delim("JH02_NT_D0_enriched_M0.txt")
#intersection
#D0_rhizo_and_M0_enriched <- as.vector(intersect(D0_enriched_M0[, 1], D0_rhizo_enriched[, 1]))

#proper method


############################################################################################################################
#significance at padj < 0.05
############################################################################################################################

#D0 rhizosphere 669
DESeq_OUTPUT_rhizosphere_D0 <- read.delim("JH02_NT_microhabitat_B0_D0.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_rhizosphere_FDR005_D0 <-DESeq_OUTPUT_rhizosphere_D0[(rownames(DESeq_OUTPUT_rhizosphere_D0)[which(DESeq_OUTPUT_rhizosphere_D0$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_enriched_D0 <- DESeq_OUTPUT_rhizosphere_D0[(rownames(DESeq_OUTPUT_rhizosphere_D0)[which(DESeq_OUTPUT_rhizosphere_D0$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_FDR005_D0 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_D0), rownames(DESeq_OUTPUT_rhizosphere_enriched_D0))
length(DESeq_OUTPUT_rhizosphere_enriched_FDR005_D0)

#N0 rhizosphere 399
DESeq_OUTPUT_rhizosphere_N0 <- read.delim("JH02_NT_microhabitat_B0_N0.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_rhizosphere_FDR005_N0 <-DESeq_OUTPUT_rhizosphere_N0[(rownames(DESeq_OUTPUT_rhizosphere_N0)[which(DESeq_OUTPUT_rhizosphere_N0$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_enriched_N0 <- DESeq_OUTPUT_rhizosphere_N0[(rownames(DESeq_OUTPUT_rhizosphere_N0)[which(DESeq_OUTPUT_rhizosphere_N0$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_FDR005_N0 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_N0), rownames(DESeq_OUTPUT_rhizosphere_enriched_N0))
length(DESeq_OUTPUT_rhizosphere_enriched_FDR005_N0)
#M0 rhizosphere 729
DESeq_OUTPUT_rhizosphere_M0 <- read.delim("JH02_NT_microhabitat_B0_M0.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_rhizosphere_FDR005_M0 <-DESeq_OUTPUT_rhizosphere_M0[(rownames(DESeq_OUTPUT_rhizosphere_M0)[which(DESeq_OUTPUT_rhizosphere_M0$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_enriched_M0 <- DESeq_OUTPUT_rhizosphere_M0[(rownames(DESeq_OUTPUT_rhizosphere_M0)[which(DESeq_OUTPUT_rhizosphere_M0$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_FDR005_M0 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_M0), rownames(DESeq_OUTPUT_rhizosphere_enriched_M0))
length(DESeq_OUTPUT_rhizosphere_enriched_FDR005_M0)
#D25 rhizosphere 512
DESeq_OUTPUT_rhizosphere_D25 <- read.delim("JH02_NT_microhabitat_B25_D25.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_rhizosphere_FDR005_D25 <-DESeq_OUTPUT_rhizosphere_D25[(rownames(DESeq_OUTPUT_rhizosphere_D25)[which(DESeq_OUTPUT_rhizosphere_D25$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_enriched_D25 <- DESeq_OUTPUT_rhizosphere_D25[(rownames(DESeq_OUTPUT_rhizosphere_D25)[which(DESeq_OUTPUT_rhizosphere_D25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_FDR005_D25 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_D25), rownames(DESeq_OUTPUT_rhizosphere_enriched_D25))
length(DESeq_OUTPUT_rhizosphere_enriched_FDR005_D25)
#N25 rhizosphere 400
DESeq_OUTPUT_rhizosphere_N25 <- read.delim("JH02_NT_microhabitat_B25_N25.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_rhizosphere_FDR005_N25 <-DESeq_OUTPUT_rhizosphere_N25[(rownames(DESeq_OUTPUT_rhizosphere_N25)[which(DESeq_OUTPUT_rhizosphere_N25$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_enriched_N25 <- DESeq_OUTPUT_rhizosphere_N25[(rownames(DESeq_OUTPUT_rhizosphere_N25)[which(DESeq_OUTPUT_rhizosphere_N25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_FDR005_N25 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_N25), rownames(DESeq_OUTPUT_rhizosphere_enriched_N0))
#M25 rhizosphere 510
DESeq_OUTPUT_rhizosphere_M25 <- read.delim("JH02_NT_microhabitat_B25_M25.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_rhizosphere_FDR005_M25 <-DESeq_OUTPUT_rhizosphere_M25[(rownames(DESeq_OUTPUT_rhizosphere_M25)[which(DESeq_OUTPUT_rhizosphere_M25$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_enriched_M25 <- DESeq_OUTPUT_rhizosphere_M25[(rownames(DESeq_OUTPUT_rhizosphere_M25)[which(DESeq_OUTPUT_rhizosphere_M25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_FDR005_M25 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_M25), rownames(DESeq_OUTPUT_rhizosphere_enriched_M25))
#D100 rhizosphere 479
DESeq_OUTPUT_rhizosphere_D100 <- read.delim("JH02_NT_microhabitat_B100_D100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_rhizosphere_FDR005_D100 <-DESeq_OUTPUT_rhizosphere_D100[(rownames(DESeq_OUTPUT_rhizosphere_D100)[which(DESeq_OUTPUT_rhizosphere_D100$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_enriched_D100 <- DESeq_OUTPUT_rhizosphere_D100[(rownames(DESeq_OUTPUT_rhizosphere_D100)[which(DESeq_OUTPUT_rhizosphere_D100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_FDR005_D100 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_D100), rownames(DESeq_OUTPUT_rhizosphere_enriched_D100))
#N100 rhizosphere 363
DESeq_OUTPUT_rhizosphere_N100 <- read.delim("JH02_NT_microhabitat_B100_N100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_rhizosphere_FDR005_N100 <-DESeq_OUTPUT_rhizosphere_N100[(rownames(DESeq_OUTPUT_rhizosphere_N100)[which(DESeq_OUTPUT_rhizosphere_N100$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_enriched_N100 <- DESeq_OUTPUT_rhizosphere_N100[(rownames(DESeq_OUTPUT_rhizosphere_N100)[which(DESeq_OUTPUT_rhizosphere_N100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_FDR005_N100 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_N100), rownames(DESeq_OUTPUT_rhizosphere_enriched_N100))
#M100 rhizosphere 407
DESeq_OUTPUT_rhizosphere_M100 <- read.delim("JH02_NT_microhabitat_B100_M100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_rhizosphere_FDR005_M100 <-DESeq_OUTPUT_rhizosphere_M100[(rownames(DESeq_OUTPUT_rhizosphere_M100)[which(DESeq_OUTPUT_rhizosphere_M100$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_enriched_M100 <- DESeq_OUTPUT_rhizosphere_M100[(rownames(DESeq_OUTPUT_rhizosphere_M100)[which(DESeq_OUTPUT_rhizosphere_M100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_FDR005_M100 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_M100), rownames(DESeq_OUTPUT_rhizosphere_enriched_M100))


########################################## 0% #############################################

#Genotype enriched
#D0 vs M0 (54-44)
DESeq_OUTPUT_D0_M0 <- read.delim("JH02_NT_genotype_D0_M0.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_D0_M0 <- DESeq_OUTPUT_D0_M0[(rownames(DESeq_OUTPUT_D0_M0)[which(DESeq_OUTPUT_D0_M0$padj <0.05)]), ]
DESeq_OUTPUT_enriched_D0_M0 <- DESeq_OUTPUT_D0_M0[(rownames(DESeq_OUTPUT_D0_M0)[which(DESeq_OUTPUT_D0_M0$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_M0_D0 <-  DESeq_OUTPUT_D0_M0[(rownames(DESeq_OUTPUT_D0_M0)[which(DESeq_OUTPUT_D0_M0$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D0_M0 <- intersect(rownames(DESeq_OUTPUT_FDR005_D0_M0), rownames(DESeq_OUTPUT_enriched_D0_M0))
DESeq_OUTPUT_enriched_FDR005_M0_D0 <- intersect(rownames(DESeq_OUTPUT_FDR005_D0_M0), rownames(DESeq_OUTPUT_enriched_M0_D0))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_D0_M0 <- intersect(DESeq_OUTPUT_enriched_FDR005_D0_M0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D0_M0_taxonomy <- DESeq_OUTPUT_D0_M0[DESeq_OUTPUT_enriched_rhizo_FDR005_D0_M0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D0_M0_taxonomy, file="DESeq_OUTPUT_enriched_rhizo_FDR005_D0_M0_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_M0_D0 <- intersect(DESeq_OUTPUT_enriched_FDR005_M0_D0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M0_D0_taxonomy <- DESeq_OUTPUT_D0_M0[DESeq_OUTPUT_enriched_rhizo_FDR005_M0_D0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M0_D0_taxonomy, file="DESeq_OUTPUT_enriched_rhizo_FDR005_M0_D0_taxonomy.txt") 

#N0 vs D0 (1-52)
DESeq_OUTPUT_N0_D0 <- read.delim("JH02_NT_genotype_N0_D0.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_N0_D0 <- DESeq_OUTPUT_N0_D0[(rownames(DESeq_OUTPUT_N0_D0)[which(DESeq_OUTPUT_N0_D0$padj <0.05)]), ]
DESeq_OUTPUT_enriched_N0_D0 <- DESeq_OUTPUT_N0_D0[(rownames(DESeq_OUTPUT_N0_D0)[which(DESeq_OUTPUT_N0_D0$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N0_D0 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_D0), rownames(DESeq_OUTPUT_enriched_N0_D0))
DESeq_OUTPUT_enriched_D0_N0 <-  DESeq_OUTPUT_N0_D0[(rownames(DESeq_OUTPUT_N0_D0)[which(DESeq_OUTPUT_N0_D0$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D0_N0 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_D0), rownames(DESeq_OUTPUT_enriched_D0_N0))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_N0_D0 <- intersect(DESeq_OUTPUT_enriched_FDR005_N0_D0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N0_D0_taxonomy <- DESeq_OUTPUT_N0_D0[DESeq_OUTPUT_enriched_rhizo_FDR005_N0_D0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N0_D0_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_N0_D0_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_D0_N0 <- intersect(DESeq_OUTPUT_enriched_FDR005_D0_N0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D0_N0_taxonomy <- DESeq_OUTPUT_N0_D0[DESeq_OUTPUT_enriched_rhizo_FDR005_D0_N0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D0_N0_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_D0_N0_taxonomy.txt")

#N0 vs M0 (9-64)
DESeq_OUTPUT_N0_M0 <- read.delim("JH02_NT_genotype_N0_M0.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_N0_M0<- DESeq_OUTPUT_N0_M0[(rownames(DESeq_OUTPUT_N0_M0)[which(DESeq_OUTPUT_N0_M0$padj <0.05)]), ]
DESeq_OUTPUT_enriched_N0_M0 <- DESeq_OUTPUT_N0_M0[(rownames(DESeq_OUTPUT_N0_M0)[which(DESeq_OUTPUT_N0_M0$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N0_M0 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_M0), rownames(DESeq_OUTPUT_enriched_N0_M0))
DESeq_OUTPUT_enriched_M0_N0 <-  DESeq_OUTPUT_N0_M0[(rownames(DESeq_OUTPUT_N0_M0)[which(DESeq_OUTPUT_N0_M0$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_M0_N0 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_M0), rownames(DESeq_OUTPUT_enriched_M0_N0))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_N0_M0 <- intersect(DESeq_OUTPUT_enriched_FDR005_N0_M0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N0_M0_taxonomy <- DESeq_OUTPUT_N0_M0[DESeq_OUTPUT_enriched_rhizo_FDR005_N0_M0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N0_M0_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_N0_M0_taxonomy.txt" )
DESeq_OUTPUT_enriched_rhizo_FDR005_M0_N0 <- intersect(DESeq_OUTPUT_enriched_FDR005_M0_N0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M0_N0_taxonomy <- DESeq_OUTPUT_N0_M0[DESeq_OUTPUT_enriched_rhizo_FDR005_M0_N0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M0_N0_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_M0_N0_taxonomy.txt" )

#####################25%#########################################################################################
#D25 vs M25 (10-12)
DESeq_OUTPUT_D25_M25 <- read.delim("JH02_NT_genotype_D25_M25.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_D25_M25 <- DESeq_OUTPUT_D25_M25[(rownames(DESeq_OUTPUT_D25_M25)[which(DESeq_OUTPUT_D25_M25$padj <0.05)]), ]
DESeq_OUTPUT_enriched_D25_M25 <- DESeq_OUTPUT_D25_M25[(rownames(DESeq_OUTPUT_D25_M25)[which(DESeq_OUTPUT_D25_M25$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D25_M25 <- intersect(rownames(DESeq_OUTPUT_FDR005_D25_M25), rownames(DESeq_OUTPUT_enriched_D25_M25))
DESeq_OUTPUT_enriched_M25_D25 <-  DESeq_OUTPUT_D25_M25[(rownames(DESeq_OUTPUT_D25_M25)[which(DESeq_OUTPUT_D25_M25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_M25_D25 <- intersect(rownames(DESeq_OUTPUT_FDR005_D25_M25), rownames(DESeq_OUTPUT_enriched_M25_D25))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_D25_M25 <- intersect(DESeq_OUTPUT_enriched_FDR005_D25_M25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D25_M25_taxonomy <- DESeq_OUTPUT_D25_M25[DESeq_OUTPUT_enriched_rhizo_FDR005_D25_M25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D25_M25_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_D25_M25_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_M25_D25 <- intersect(DESeq_OUTPUT_enriched_FDR005_M25_D25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M25_D25_taxonomy <- DESeq_OUTPUT_D25_M25[DESeq_OUTPUT_enriched_rhizo_FDR005_M25_D25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M25_D25_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_M25_D25_taxonomy.txt")

#N25 vs D25 (0-2)
DESeq_OUTPUT_N25_D25 <- read.delim("JH02_NT_genotype_N25_D25.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_N25_D25 <- DESeq_OUTPUT_N25_D25[(rownames(DESeq_OUTPUT_N25_D25)[which(DESeq_OUTPUT_N25_D25$padj <0.05)]), ]
DESeq_OUTPUT_enriched_N25_D25 <- DESeq_OUTPUT_N25_D25[(rownames(DESeq_OUTPUT_N25_D25)[which(DESeq_OUTPUT_N25_D25$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N25_D25 <- intersect(rownames(DESeq_OUTPUT_FDR005_N25_D25), rownames(DESeq_OUTPUT_enriched_N25_D25))
DESeq_OUTPUT_enriched_D25_N25 <-  DESeq_OUTPUT_N25_D25[(rownames(DESeq_OUTPUT_N25_D25)[which(DESeq_OUTPUT_N25_D25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D25_N25 <- intersect(rownames(DESeq_OUTPUT_FDR005_N25_D25), rownames(DESeq_OUTPUT_enriched_D25_N25))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_N25_D25 <- intersect(DESeq_OUTPUT_enriched_FDR005_N25_D25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N25_D25_taxonomy <- DESeq_OUTPUT_N25_D25[DESeq_OUTPUT_enriched_rhizo_FDR005_N25_D25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N25_D25_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_N25_D25_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_D25_N25 <- intersect(DESeq_OUTPUT_enriched_FDR005_D25_N25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D25_N25_taxonomy <- DESeq_OUTPUT_N25_D25[DESeq_OUTPUT_enriched_rhizo_FDR005_D25_N25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D25_N25_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_D25_N25_taxonomy.txt")

#N25 vs M25 (12-16)
DESeq_OUTPUT_N25_M25 <- read.delim("JH02_NT_genotype_N25_M25.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_N25_M25<- DESeq_OUTPUT_N25_M25[(rownames(DESeq_OUTPUT_N25_M25)[which(DESeq_OUTPUT_N25_M25$padj <0.05)]), ]
DESeq_OUTPUT_enriched_N25_M25 <- DESeq_OUTPUT_N25_M25[(rownames(DESeq_OUTPUT_N25_M25)[which(DESeq_OUTPUT_N25_M25$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N25_M25 <- intersect(rownames(DESeq_OUTPUT_FDR005_N25_M25), rownames(DESeq_OUTPUT_enriched_N25_M25))
DESeq_OUTPUT_enriched_M25_N25 <-  DESeq_OUTPUT_N25_M25[(rownames(DESeq_OUTPUT_N25_M25)[which(DESeq_OUTPUT_N25_M25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_M25_N25 <- intersect(rownames(DESeq_OUTPUT_FDR005_N25_M25), rownames(DESeq_OUTPUT_enriched_M25_N25))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_N25_M25 <- intersect(DESeq_OUTPUT_enriched_FDR005_N25_M25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N25_M25_taxonomy <- DESeq_OUTPUT_N25_M25[DESeq_OUTPUT_enriched_rhizo_FDR005_N25_M25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N25_M25_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_N25_M25_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_M25_N25 <- intersect(DESeq_OUTPUT_enriched_FDR005_M25_N25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M25_N25_taxonomy <- DESeq_OUTPUT_N25_M25[DESeq_OUTPUT_enriched_rhizo_FDR005_M25_N25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M25_N25_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_M25_N25_taxonomy.txt")
  
############################################### 100% #####################################################################
#Genotype enriched
#D100 vs M100 (6-0)
DESeq_OUTPUT_D100_M100 <- read.delim("JH02_NT_genotype_D100_M100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_D100_M100 <- DESeq_OUTPUT_D100_M100[(rownames(DESeq_OUTPUT_D100_M100)[which(DESeq_OUTPUT_D100_M100$padj <0.05)]), ]
DESeq_OUTPUT_enriched_D100_M100 <- DESeq_OUTPUT_D100_M100[(rownames(DESeq_OUTPUT_D100_M100)[which(DESeq_OUTPUT_D100_M100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_M100_D100 <-  DESeq_OUTPUT_D100_M100[(rownames(DESeq_OUTPUT_D100_M100)[which(DESeq_OUTPUT_D100_M100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D100_M100 <- intersect(rownames(DESeq_OUTPUT_FDR005_D100_M100), rownames(DESeq_OUTPUT_enriched_D100_M100))
DESeq_OUTPUT_enriched_FDR005_M100_D100 <- intersect(rownames(DESeq_OUTPUT_FDR005_D100_M100), rownames(DESeq_OUTPUT_enriched_M100_D100))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_D100_M100 <- intersect(DESeq_OUTPUT_enriched_FDR005_D100_M100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D100_M100_taxonomy <- DESeq_OUTPUT_D100_M100[DESeq_OUTPUT_enriched_rhizo_FDR005_D100_M100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D100_M100_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_D100_M100_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_M100_D100 <- intersect(DESeq_OUTPUT_enriched_FDR005_M100_D100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M100_D100_taxonomy <- DESeq_OUTPUT_D100_M100[DESeq_OUTPUT_enriched_rhizo_FDR005_M100_D100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M100_D100_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_M100_D100_taxonomy.txt")

#N100 vs D100 (1-5)
DESeq_OUTPUT_N100_D100 <- read.delim("JH02_NT_genotype_N100_D100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_N100_D100 <- DESeq_OUTPUT_N100_D100[(rownames(DESeq_OUTPUT_N100_D100)[which(DESeq_OUTPUT_N100_D100$padj <0.05)]), ]
DESeq_OUTPUT_enriched_N100_D100 <- DESeq_OUTPUT_N100_D100[(rownames(DESeq_OUTPUT_N100_D100)[which(DESeq_OUTPUT_N100_D100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N100_D100 <- intersect(rownames(DESeq_OUTPUT_FDR005_N100_D100), rownames(DESeq_OUTPUT_enriched_N100_D100))
DESeq_OUTPUT_enriched_D100_N100 <-  DESeq_OUTPUT_N100_D100[(rownames(DESeq_OUTPUT_N100_D100)[which(DESeq_OUTPUT_N100_D100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D100_N100 <- intersect(rownames(DESeq_OUTPUT_FDR005_N100_D100), rownames(DESeq_OUTPUT_enriched_D100_N100))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_N100_D100 <- intersect(DESeq_OUTPUT_enriched_FDR005_N100_D100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N100_D100_taxonomy <- DESeq_OUTPUT_N100_D100[DESeq_OUTPUT_enriched_rhizo_FDR005_N100_D100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N100_D100_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_N100_D100_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_D100_N100 <- intersect(DESeq_OUTPUT_enriched_FDR005_D100_N100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D100_N100_taxonomy <- DESeq_OUTPUT_N100_D100[DESeq_OUTPUT_enriched_rhizo_FDR005_D100_N100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D100_N100_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_D100_N100_taxonomy.txt")

#N100 vs M100 (2)
DESeq_OUTPUT_N100_M100 <- read.delim("JH02_NT_genotype_N100_M100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_N100_M100<- DESeq_OUTPUT_N100_M100[(rownames(DESeq_OUTPUT_N100_M100)[which(DESeq_OUTPUT_N100_M100$padj <0.05)]), ]
DESeq_OUTPUT_enriched_N100_M100 <- DESeq_OUTPUT_N100_M100[(rownames(DESeq_OUTPUT_N100_M100)[which(DESeq_OUTPUT_N100_M100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N100_M100 <- intersect(rownames(DESeq_OUTPUT_FDR005_N100_M100), rownames(DESeq_OUTPUT_enriched_N100_M100))
DESeq_OUTPUT_enriched_M100_N100 <-  DESeq_OUTPUT_N100_M100[(rownames(DESeq_OUTPUT_N100_M100)[which(DESeq_OUTPUT_N100_M100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_M100_N100 <- intersect(rownames(DESeq_OUTPUT_FDR005_N100_M100), rownames(DESeq_OUTPUT_enriched_M100_N100))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched in paisr wise comparasion
DESeq_OUTPUT_enriched_rhizo_FDR005_N100_M100 <- intersect(DESeq_OUTPUT_enriched_FDR005_N100_M100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N100_M100_taxonomy <- DESeq_OUTPUT_N100_M100[DESeq_OUTPUT_enriched_rhizo_FDR005_N100_M100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N100_M100_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_N100_M100_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_M100_N100 <- intersect(DESeq_OUTPUT_enriched_FDR005_M100_N100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M100_N100_taxonomy <- DESeq_OUTPUT_N100_M100[DESeq_OUTPUT_enriched_rhizo_FDR005_M100_N100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M100_N100_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_M100_N100_taxonomy.txt")

##########################################################################################################################################################
############################################Treatment################################################################################################

#N0_N25 (2-6)
DESeq_OUTPUT_N0_N25 <- read.delim("JH02_NT_treatment_N0_N25.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_N0_N25 <- DESeq_OUTPUT_N0_N25[(rownames(DESeq_OUTPUT_N0_N25)[which(DESeq_OUTPUT_N0_N25$padj <0.05)]), ]
DESeq_OUTPUT_enriched_N0_N25 <- DESeq_OUTPUT_N0_N25[(rownames(DESeq_OUTPUT_N0_N25)[which(DESeq_OUTPUT_N0_N25$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N0_N25 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_N25), rownames(DESeq_OUTPUT_enriched_N0_N25))
DESeq_OUTPUT_enriched_N25_N0 <-  DESeq_OUTPUT_N0_N25[(rownames(DESeq_OUTPUT_N0_N25)[which(DESeq_OUTPUT_N0_N25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N25_N0 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_N25), rownames(DESeq_OUTPUT_enriched_N25_N0))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_N0_N25 <- intersect(DESeq_OUTPUT_enriched_FDR005_N0_N25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N0_N25_taxonomy <- DESeq_OUTPUT_N0_N25[DESeq_OUTPUT_enriched_rhizo_FDR005_N0_N25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N0_N25_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_N0_N25_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_N25_N0 <- intersect(DESeq_OUTPUT_enriched_FDR005_N25_N0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N25_N0_taxonomy <- DESeq_OUTPUT_N0_N25[DESeq_OUTPUT_enriched_rhizo_FDR005_N25_N0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N25_N0_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_N25_N0_taxonomy.txt")

#N0_N100 (11-21)
DESeq_OUTPUT_N0_N100 <- read.delim("JH02_NT_treatment_N0_N100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_N0_N100 <- DESeq_OUTPUT_N0_N100[(rownames(DESeq_OUTPUT_N0_N100)[which(DESeq_OUTPUT_N0_N100$padj <0.05)]), ]
DESeq_OUTPUT_enriched_N0_N100 <- DESeq_OUTPUT_N0_N100[(rownames(DESeq_OUTPUT_N0_N100)[which(DESeq_OUTPUT_N0_N100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N0_N100 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_N100), rownames(DESeq_OUTPUT_enriched_N0_N100))
DESeq_OUTPUT_enriched_N100_N0 <-  DESeq_OUTPUT_N0_N100[(rownames(DESeq_OUTPUT_N0_N100)[which(DESeq_OUTPUT_N0_N100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N100_N0 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_N100), rownames(DESeq_OUTPUT_enriched_N100_N0))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_N0_N100 <- intersect(DESeq_OUTPUT_enriched_FDR005_N0_N100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N0_N100_taxonomy <- DESeq_OUTPUT_N0_N100[DESeq_OUTPUT_enriched_rhizo_FDR005_N0_N100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N0_N100_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_N0_N100_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_N100_N0 <- intersect(DESeq_OUTPUT_enriched_FDR005_N100_N0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N100_N0_taxonomy <- DESeq_OUTPUT_N0_N100[DESeq_OUTPUT_enriched_rhizo_FDR005_N100_N0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N100_N0_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_N100_N0_taxonomy.txt")

#N25_N100 (6-3)
DESeq_OUTPUT_N25_N100 <- read.delim("JH02_NT_treatment_N25_N100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_N25_N100 <- DESeq_OUTPUT_N25_N100[(rownames(DESeq_OUTPUT_N25_N100)[which(DESeq_OUTPUT_N25_N100$padj <0.05)]), ]
DESeq_OUTPUT_enriched_N25_N100 <- DESeq_OUTPUT_N25_N100[(rownames(DESeq_OUTPUT_N25_N100)[which(DESeq_OUTPUT_N25_N100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N25_N100 <- intersect(rownames(DESeq_OUTPUT_FDR005_N25_N100), rownames(DESeq_OUTPUT_enriched_N25_N100))
DESeq_OUTPUT_enriched_N100_N25 <-  DESeq_OUTPUT_N25_N100[(rownames(DESeq_OUTPUT_N25_N100)[which(DESeq_OUTPUT_N25_N100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_N100_N25 <- intersect(rownames(DESeq_OUTPUT_FDR005_N25_N100), rownames(DESeq_OUTPUT_enriched_N100_N25))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_N25_N100 <- intersect(DESeq_OUTPUT_enriched_FDR005_N25_N100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N25_N100_taxonomy <- DESeq_OUTPUT_N25_N100[DESeq_OUTPUT_enriched_rhizo_FDR005_N25_N100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N25_N100_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_N25_N100_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_N100_N25 <- intersect(DESeq_OUTPUT_enriched_FDR005_N100_N25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_N100_N25_taxonomy <- DESeq_OUTPUT_N25_N100[DESeq_OUTPUT_enriched_rhizo_FDR005_N100_N25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_N100_N25_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_N100_N25_taxonomy.txt")

#D0_D25 (115-33)
DESeq_OUTPUT_D0_D25 <- read.delim("JH02_NT_treatment_D0_D25.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_D0_D25 <- DESeq_OUTPUT_D0_D25[(rownames(DESeq_OUTPUT_D0_D25)[which(DESeq_OUTPUT_D0_D25$padj <0.05)]), ]
DESeq_OUTPUT_enriched_D0_D25 <- DESeq_OUTPUT_D0_D25[(rownames(DESeq_OUTPUT_D0_D25)[which(DESeq_OUTPUT_D0_D25$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D0_D25 <- intersect(rownames(DESeq_OUTPUT_FDR005_D0_D25), rownames(DESeq_OUTPUT_enriched_D0_D25))
DESeq_OUTPUT_enriched_D25_D0 <-  DESeq_OUTPUT_D0_D25[(rownames(DESeq_OUTPUT_D0_D25)[which(DESeq_OUTPUT_D0_D25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D25_D0 <- intersect(rownames(DESeq_OUTPUT_FDR005_D0_D25), rownames(DESeq_OUTPUT_enriched_D25_D0))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_D0_D25 <- intersect(DESeq_OUTPUT_enriched_FDR005_D0_D25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D0_D25_taxonomy <- DESeq_OUTPUT_D0_D25[DESeq_OUTPUT_enriched_rhizo_FDR005_D0_D25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D0_D25_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_D0_D25_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_D25_D0 <- intersect(DESeq_OUTPUT_enriched_FDR005_D25_D0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D25_D0_taxonomy <- DESeq_OUTPUT_D0_D25[DESeq_OUTPUT_enriched_rhizo_FDR005_D25_D0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D25_D0_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_D25_D0_taxonomy.txt")


#D0_D100 (242-37)
DESeq_OUTPUT_D0_D100 <- read.delim("JH02_NT_treatment_D0_D100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_D0_D100 <- DESeq_OUTPUT_D0_D100[(rownames(DESeq_OUTPUT_D0_D100)[which(DESeq_OUTPUT_D0_D100$padj <0.05)]), ]
DESeq_OUTPUT_enriched_D0_D100 <- DESeq_OUTPUT_D0_D100[(rownames(DESeq_OUTPUT_D0_D100)[which(DESeq_OUTPUT_D0_D100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D0_D100 <- intersect(rownames(DESeq_OUTPUT_FDR005_D0_D100), rownames(DESeq_OUTPUT_enriched_D0_D100))
DESeq_OUTPUT_enriched_D100_D0 <-  DESeq_OUTPUT_D0_D100[(rownames(DESeq_OUTPUT_D0_D100)[which(DESeq_OUTPUT_D0_D100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D100_D0  <- intersect(rownames(DESeq_OUTPUT_FDR005_D0_D100 ), rownames(DESeq_OUTPUT_enriched_D100_D0 ))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_D0_D100 <- intersect(DESeq_OUTPUT_enriched_FDR005_D0_D100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D0_D100_taxonomy <- DESeq_OUTPUT_D0_D100[DESeq_OUTPUT_enriched_rhizo_FDR005_D0_D100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D0_D100_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_D0_D100_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_D100_D0 <- intersect(DESeq_OUTPUT_enriched_FDR005_D100_D0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D100_D0_taxonomy <- DESeq_OUTPUT_D0_D100[DESeq_OUTPUT_enriched_rhizo_FDR005_D100_D0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D100_D0_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_D100_D0_taxonomy.txt")

#D25_D100 (27-8)
DESeq_OUTPUT_D25_D100 <- read.delim("JH02_NT_treatment_D25_D100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_D25_D100 <- DESeq_OUTPUT_D25_D100[(rownames(DESeq_OUTPUT_D25_D100)[which(DESeq_OUTPUT_D25_D100$padj <0.05)]), ]
DESeq_OUTPUT_enriched_D25_D100 <- DESeq_OUTPUT_D25_D100[(rownames(DESeq_OUTPUT_D25_D100)[which(DESeq_OUTPUT_D25_D100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D25_D100 <- intersect(rownames(DESeq_OUTPUT_FDR005_D25_D100), rownames(DESeq_OUTPUT_enriched_D25_D100))
DESeq_OUTPUT_enriched_D100_D25 <-  DESeq_OUTPUT_D25_D100[(rownames(DESeq_OUTPUT_D25_D100)[which(DESeq_OUTPUT_D25_D100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_D100_D25 <- intersect(rownames(DESeq_OUTPUT_FDR005_D25_D100), rownames(DESeq_OUTPUT_enriched_D100_D25))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_D25_D100 <- intersect(DESeq_OUTPUT_enriched_FDR005_D25_D100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D25_D100_taxonomy <- DESeq_OUTPUT_D25_D100[DESeq_OUTPUT_enriched_rhizo_FDR005_D25_D100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D25_D100_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_D25_D100_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_D100_D25 <- intersect(DESeq_OUTPUT_enriched_FDR005_D100_D25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D100_D25_taxonomy <- DESeq_OUTPUT_D25_D100[DESeq_OUTPUT_enriched_rhizo_FDR005_D100_D25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_D100_D25_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_D100_D25_taxonomy.txt")

#M0_M25(61-26)
DESeq_OUTPUT_M0_M25 <- read.delim("JH02_NT_treatment_M0_M25.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_M0_M25 <- DESeq_OUTPUT_M0_M25[(rownames(DESeq_OUTPUT_M0_M25)[which(DESeq_OUTPUT_M0_M25$padj <0.05)]), ]
DESeq_OUTPUT_enriched_M0_M25 <- DESeq_OUTPUT_M0_M25[(rownames(DESeq_OUTPUT_M0_M25)[which(DESeq_OUTPUT_M0_M25$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_M0_M25 <- intersect(rownames(DESeq_OUTPUT_FDR005_M0_M25), rownames(DESeq_OUTPUT_enriched_M0_M25))
DESeq_OUTPUT_enriched_M25_M0 <-  DESeq_OUTPUT_M0_M25[(rownames(DESeq_OUTPUT_M0_M25)[which(DESeq_OUTPUT_M0_M25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_M25_M0 <- intersect(rownames(DESeq_OUTPUT_FDR005_M0_M25), rownames(DESeq_OUTPUT_enriched_M25_M0))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_M0_M25 <- intersect(DESeq_OUTPUT_enriched_FDR005_M0_M25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M0_M25_taxonomy <- DESeq_OUTPUT_M0_M25[DESeq_OUTPUT_enriched_rhizo_FDR005_M0_M25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M0_M25_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_M0_M25_taxonomy.txt") 
DESeq_OUTPUT_enriched_rhizo_FDR005_M25_M0 <- intersect(DESeq_OUTPUT_enriched_FDR005_M25_M0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M25_M0_taxonomy <- DESeq_OUTPUT_M0_M25[DESeq_OUTPUT_enriched_rhizo_FDR005_M25_M0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M25_M0_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_M25_M0_taxonomy.txt") 

#M0_M100 (180-32)
DESeq_OUTPUT_M0_M100 <- read.delim("JH02_NT_treatment_M0_M100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_M0_M100 <- DESeq_OUTPUT_M0_M100[(rownames(DESeq_OUTPUT_M0_M100)[which(DESeq_OUTPUT_M0_M100$padj <0.05)]), ]
DESeq_OUTPUT_enriched_M0_M100 <- DESeq_OUTPUT_M0_M100[(rownames(DESeq_OUTPUT_M0_M100)[which(DESeq_OUTPUT_M0_M100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_M0_M100 <- intersect(rownames(DESeq_OUTPUT_FDR005_M0_M100), rownames(DESeq_OUTPUT_enriched_M0_M100))
DESeq_OUTPUT_enriched_M100_M0 <-  DESeq_OUTPUT_M0_M100[(rownames(DESeq_OUTPUT_M0_M100)[which(DESeq_OUTPUT_M0_M100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_M100_M0 <- intersect(rownames(DESeq_OUTPUT_FDR005_M0_M100), rownames(DESeq_OUTPUT_enriched_M100_M0))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_M0_M100 <- intersect(DESeq_OUTPUT_enriched_FDR005_M0_M100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M0_M100_taxonomy <- DESeq_OUTPUT_M0_M100[DESeq_OUTPUT_enriched_rhizo_FDR005_M0_M100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M0_M100_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_M0_M100_taxonomy.txt")
DESeq_OUTPUT_enriched_rhizo_FDR005_M100_M0 <- intersect(DESeq_OUTPUT_enriched_FDR005_M100_M0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M100_M0_taxonomy <- DESeq_OUTPUT_M0_M100[DESeq_OUTPUT_enriched_rhizo_FDR005_M100_M0,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M100_M0_taxonomy, "DESeq_OUTPUT_enriched_rhizo_FDR005_M100_M0_taxonomy.txt")

#M25_M100 (11-3)
DESeq_OUTPUT_M25_M100 <- read.delim("JH02_NT_treatment_M25_M100.txt", row.names = 1, header = TRUE)
DESeq_OUTPUT_FDR005_M25_M100 <- DESeq_OUTPUT_M25_M100[(rownames(DESeq_OUTPUT_M25_M100)[which(DESeq_OUTPUT_M25_M100$padj <0.05)]), ]
DESeq_OUTPUT_enriched_M25_M100 <- DESeq_OUTPUT_M25_M100[(rownames(DESeq_OUTPUT_M25_M100)[which(DESeq_OUTPUT_M25_M100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_FDR005_M25_M100 <- intersect(rownames(DESeq_OUTPUT_FDR005_M25_M100), rownames(DESeq_OUTPUT_enriched_M25_M100))
DESeq_OUTPUT_enriched_M100_M25 <-  DESeq_OUTPUT_M25_M100[(rownames(DESeq_OUTPUT_M25_M100)[which(DESeq_OUTPUT_M25_M100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_enriched_FDR005_M100_M25 <- intersect(rownames(DESeq_OUTPUT_FDR005_M25_M100), rownames(DESeq_OUTPUT_enriched_M100_M25))
#intersection to identify OTUs a) enriched vs rhizosphere and b) differentially enriched Desert vs Morex
DESeq_OUTPUT_enriched_rhizo_FDR005_M25_M100 <- intersect(DESeq_OUTPUT_enriched_FDR005_M25_M100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M25_M100_taxonomy <- DESeq_OUTPUT_M25_M100[DESeq_OUTPUT_enriched_rhizo_FDR005_M25_M100,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M25_M100_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_M25_M100_taxonomy.txt") 
DESeq_OUTPUT_enriched_rhizo_FDR005_M100_M25 <- intersect(DESeq_OUTPUT_enriched_FDR005_M100_M25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M100_M25_taxonomy <- DESeq_OUTPUT_M25_M100[DESeq_OUTPUT_enriched_rhizo_FDR005_M100_M25,  ]
#write.table(DESeq_OUTPUT_enriched_rhizo_FDR005_M100_M25_taxonomy,"DESeq_OUTPUT_enriched_rhizo_FDR005_M100_M25_taxonomy.txt")

#########################################################################################################################################
#TAX4FUN
#PREDICTION OF FUNCTIONAL CAPABILITIES OF THE MICROBIOME
#data required: JH02_NT_silva115_otu_table.txt, 
#silva database 115
##########################################################################################################################################

OTU_TABLE_tax4fun <- importQIIMEData("JH02_NT_silva115_otu_table.txt")
folderReferenceData <- ("C:/Users/rodr/R files/Rodrigo/SILVA115/SILVA115")
Tax4FunOutput <- Tax4Fun(OTU_TABLE_tax4fun, folderReferenceData, fctProfiling = F, refProfile = "UProC", shortReadMode = TRUE, normCopyNo = TRUE)
Tax4FunProfile<-Tax4FunOutput$Tax4FunProfile
Tax4FunProfile <- data.frame(t(Tax4FunOutput$Tax4FunProfile))
write.table(Tax4FunProfile, file="Tax4FunProfile.txt", sep="\t")





#############################################################
#                     ANCOM analysis                        #
#############################################################

# Import OTU table in 'wide format', generated manually in Excel.
OTU_table_transposed <- read.delim("JH02_NT_otu_counts_L2.txt")

# ANCOM calculation
ancom_JH02_NT_phyla <- ANCOM(OTU_table_transposed, multcorr = 2)

# Generate data frame with output of the calculation
ancom_result_dataframe <- as.data.frame(ancom_JH02_NT_phyla$detected)

# Save result of ANCOM analysis as data frame
write.table(ancom_result_dataframe, file="ancom_output.txt", sep="\t")

# End

