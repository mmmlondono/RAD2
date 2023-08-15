####Using https://fukamilab.github.io/BIO202/06-C-matrix-comparison.html#matrix_comparisons####

# Load packages
library(readr)
library(dplyr)
library(phyloseq)
library(phyloseqGraphTest)
library(vegan)
library(ade4)
library(ggplot2)

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("phyloseq")

#Mantel test
#Compare two distance/similarity matrices that were obtained independently of each other

# Transpose otu tables from both phyloseqs and convert to data frame
data_16s <- as.data.frame(t(otu_table(merged_16s)))
data_CO1 <- as.data.frame(t(otu_table(merged_CO1)))

# Calculate distance using Bray-Curtis
DM16s <- distance(merged_16s,"bray")
DMco1 <- distance(merged_CO1, "bray")

# Mantel test (using three different methods)
mantel(DMco1, DM16s, method = "pearson", permutations = 999)    


#### other way from https://www.molecularecologist.com/2015/04/23/procrustes-analyses-in-r/ ####
#Load packages
library(MCMCpack)
library(maps)
nov<-read.csv("pca_stats/PC_mix.csv",header=TRUE)
X<-as.matrix(cbind(nov$PC1a,nov$PC2a))
Xstar<-as.matrix(cbind(nov$PC1r,nov$PC2r))
p<-procrustes(Xstar,X,translation=TRUE,dilation=TRUE)

# Access the transformed matrix and other results
X.new <- p$X.new
R <- p$R
tt <- p$tt
s <- p$s
