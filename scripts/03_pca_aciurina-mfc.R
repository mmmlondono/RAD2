####READ ME####
#this is the pca set up for aciurina mfc data set
####load packages ####
library("adegenet")
library("ade4")
library("vcfR")
library("parallel")
library("viridis")
####adegenet’s glPCA####
#This is run in Rstudio (it being an IDE is especially nice for visualizing)
#Read in VCF
mcf_vcf <- read.vcfR("aciurina/vcf-mfc/mfc.75.filtered.vcf")

#Convert to genlight format
focal_gl <- vcfR2genlight(focal_vcf)
#populations map file
popmap<-data.frame(id=colnames(focal_vcf@gt)[2:length(colnames(focal_vcf@gt))],pop=substr(colnames(focal_vcf@gt)[2:length(colnames(focal_vcf@gt))], 1,8))

#Add population info, setting up the population array
pop(focal_gl) <- c("BOU_FY_C", "BOU_FY_C", "BOU_FY_C", 
                   "GUN_FY_F", "GUN_FY_F", "GUN_FY_F", 
                   "MAL_FY_C", "MAL_FY_C", "MAL_FY_C", 
                   "MAL_FY_M", "MAL_FY_M", "MAL_FY_M", 
                   "MOB_FY_C", "MOB_FY_C", "MOB_FY_C", 
                   "MTB_FY_F", "MTB_FY_F", "MTB_FY_F", 
                   "PAS_FY_X", "PAS_FY_X", "PAS_FY_X", 
                   "RGG_FY_C", "RGG_FY_C", "RGG_FY_C", 
                   "SLV_FY_F", "SLV_FY_F", "SLV_FY_F", 
                   "SLV_FY_M", "SLV_FY_M", "SLV_FY_M", 
                   "SNM_FY_C", "SNM_FY_C", "SNM_FY_C", 
                   "SNM_FY_M", "SNM_FY_M", "SNM_FY_M", 
                   "TES_FY_C", "TES_FY_C", "TES_FY_C", 
                   "TES_FY_M", "TES_FY_M", "TES_FY_M", 
                   "TRE_FY_F", "TRE_FY_F", "TRE_FY_F", 
                   "TRQ_FY_F", "TRQ_FY_F", "TRQ_FY_F")
#Run PCA
focal_pca <- glPca(focal_gl, n.cores=4, nf=4)
#Not really great, but can be helpful for visualizing potentially problematic samples
scatter(focal_pca, cex=.25)
#My preferred option, I generally make plots here and edit in illustrator
s.class(focal_pca$scores[,c(1,2)], pop(focal_gl), 
        col=magma(10, begin=.8, end=0), clab=, cell=2.5)
#col is the color pallet used, magma is a nice default but try viridis
#cell is effectively a confidence interval (2.5 ~= 95%)

#barplot plots the proportion of the total variance explained by given PCs
#Generally, the more a PC sticks out the more it matters
barplot(focal_pca$eig/sum(focal_pca$eig), main="eigenvalues", 
        col=heat.colors(length(focal_pca$eig)))

## Heat map of genotype
glPlot(focal_gl)

pca <- glPca(focal_gl, nf = 30)
# Quick plot
scatter(pca, posi = "none")
# Create named color vector for the legend
#legend_colors <- rainbow(10)

legend_colors <- c("darkblue", "turquoise", "skyblue","seagreen","khaki", "lightgreen", "blue",
                             "hotpink","coral", "gold", "black", "red", "purple",
                             "tan","brown", "magenta")
                             
names(legend_colors) <- c("BOU_FY_C", "GUN_FY_F", "MAL_FY_C", "MAL_FY_M", "MOB_FY_C", 
                          "MTB_FY_F", "PAS_FY_X", "RGG_FY_C", "SLV_FY_F", "SLV_FY_M", 
                          "SNM_FY_C", "SNM_FY_M", "TES_FY_C", "TES_FY_M", "TRE_FY_F", 
                          "TRQ_FY_F")
popmap$pop=as.factor(popmap$pop)
# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("center", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the second and third principal components
plot(x = pca$scores[, 2], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("bottomleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the first and third principal components
plot(x = pca$scores[, 1], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("center", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

