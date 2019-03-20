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
# script to reproduce calculations and figures associated to Additional file 5
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
library(vegan)
library(gplots)
library(marray)

rm(list=ls())

#set working directory (from DB's cpu)
setwd("C:/Users/DB42008/Box Sync/Davide_lab_manuscripts/Rodrigo_NT_2018/R_data_no16S/")

#print the session info for reproducibility of the code
sessionInfo()

######################
#Additional file 4
######################

#####################
#Figure 8
#####################

#dry weight data
datafile <- read.delim("Additional_file5_ws1.txt", sep = "\t", header=TRUE )
datafile

#test the normality for both dry weight 
shapiro.test(datafile$Dry.weight..g.)

#Two-way anova 
DW <- aov(datafile$Dry.weight..g.~ Treatment * Genotype, data = datafile)
summary(DW)

#tukey HSD test 
Tukey_HSD_DW <- TukeyHSD(DW)
Tukey_HSD_DW

#data visualistation
#re-order the factors
datafile$Soil <- ordered(datafile$Soil, levels=c("B1K12-sterilised", "B1K12", "Morex-sterilised", "Morex")) 
                                                                 
#Figure 8A
dev.off()
with(datafile, boxplot(Dry.weight..g. ~  Soil , xlab = "Sample", ylab = "dry weight g", ylim = c(0, 0.5),   main = "Transplantation effect"))
with(datafile, stripchart(Dry.weight..g. ~  Soil , xlab = "Sample", ylab = "dry weight g", vertical=T, add=T, method ="jitter", col=c('red'), pch=16))

#######################################################################
##analyses on soil environmental variables
#######################################################################

# upload the data file:mineral nitrogen
datafile_2 <- read.delim("Additional_file5_ws2.txt", sep = "\t", header=TRUE, row.names=1)
dim(datafile_2)


#generate a non-metric multidimensional scaling
datafile_2.mds <- metaMDS(datafile_2, autotransform = FALSE)

#plotting
dev.off()
#save this image with sites name
plot(datafile_2.mds, type= "t", display=c ("sites") )

#fitting environmental parameters
datafile_2.em <- envfit(datafile_2.mds, datafile_2, permutations = 5000)
datafile_2.em

#save the resulting data as Additional_file5_ws3.txt 

#Figure 8A
plot(datafile_2.em, p.max = 0.002)

##end 

