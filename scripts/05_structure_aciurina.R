####READ ME####
#this is the sNMF: population structure for aciurina

####load packages into environment####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")

# Load in additional tools for sNMF
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
# Load main library, LEA
library(LEA)
library("RColorBrewer")
library(viridis)


BlackYellow <- colorRampPalette(c("yellow", "black"))

# Plotting function used to plot sNMF output.
# Input is sNMF object, the k value, and optionally an array of colors (has default of 10).
plot_sNMF <- function(input, k_val, colors = BlackYellow(k_val+1)){
  # picks best run best on cross entropy
  best_run <- which.min(cross.entropy(input, K = k_val))
  # makes q matrix of ancestry coeffs
  q_matrix <- Q(input, K = k_val, run = best_run)
  # plots the output, space makes blank between indivs
  barplot(t(q_matrix), col = colors, border = NA, space = 0.25, xlab = "Individuals", ylab = "Admixture coefficients", horiz=FALSE)
}

vcf2geno("aciurina/aciurina2_.8.filtered.vcf", output = "aciurina_all.geno")
vanua_snmf = snmf("aciurina_all.geno", ploidy=2, 
                  K = 1:10, alpha = 10, project = "new", entropy = T, repetitions = 50)
# plot to decide optimal K, although don't put too much stock in k values
plot(vanua_snmf, cex = 1.2, col = "lightblue", pch = 19)

# here I override the default colors
par(mfrow=c(2,1), mar=c(0,2,0,0), oma=c(1,2,1,0))
plot_sNMF(vanua_snmf, 2, c("yellow", "black"))
plot_sNMF(vanua_snmf, 3, c("black", "yellow", "yellow4"))
plot_sNMF(vanua_snmf, 4)
plot_sNMF(vanua_snmf, 5)

#####don't like ethan's way of plotting, used this instead####
#https://connor-french.github.io/intro-pop-structure-r/


####load packages into environment####
library(adegenet) # for dapc
library(vcfR) # for reading in genetic data
library(tidyverse) # for manipulating and plotting data
library(LEA) # For sNMF
library(rnaturalearth) #for mapping

write.table(vanua_snmf, 
            "vanua_geno.geno",
            col.names = FALSE,
            row.names = FALSE,
            sep = "")

vanua_snmf <- snmf(input.file = "~/Desktop/anole_geno.geno",
                   K = 1:10,
                   entropy = TRUE,
                   repetitions = 5,
                   project = "new",
                   alpha = 100
)

  