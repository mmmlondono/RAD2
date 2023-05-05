####READ ME####
#this is the pca set up for aciurina smooth data set
####load packages ####
library("adegenet")
library("ade4")
library("vcfR")
library("parallel")
library("viridis")
####adegenet’s glPCA####
#This is run in Rstudio (it being an IDE is especially nice for visualizing)
#Read in VCF
smooth_vcf <- read.vcfR("aciurina/vcf-smooth/smooth.75.filtered.vcf")

#Convert to genlight format
focal_gl <- vcfR2genlight(focal_vcf)
#populations map file
popmap<-data.frame(id=colnames(focal_vcf@gt)[2:length(colnames(focal_vcf@gt))],pop=substr(colnames(focal_vcf@gt)[2:length(colnames(focal_vcf@gt))], 1,8))

#Add population info, setting up the population array
pop(focal_gl) <- c("ABQ_S", "ABQ_S", "ABQ_S", 
                   "BOU_S", "BOU_S", "BOU_S", 
                   "LOG_S", "LOG_S", "LOG_S", 
                   "MOB_G", "MOB_G", "MOB_G", 
                   "MOB_S", "MOB_S", "MOB_S", 
                   "RFY_HS", "RFY_HS", "RFY_HS",
                   "RFY_S", "RFY_S", "RFY_S", 
                   "RGG_S", "RGG_S", "RGG_S", 
                   "SNM_S", "SNM_S", "SNM_S", 
                   "TES_S", "TES_S", "TES_S")
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
legend_colors <- rainbow(10)

legend_colors <- c("darkblue", "turquoise", "skyblue","seagreen","lightgreen",
                        "hotpink","coral", "gold", "black", "red")
                        
names(legend_colors) <- c("ABQ_FY_S", "BOU_FY_S", "LOG_FY_S", "MOB_FY_G",
                          "MOB_FY_S", "RFY_FY_H", "RFY_FY_S", "RGG_FY_S", "SNM_FY_S", "TES_FY_S")
                         
popmap$pop=as.factor(popmap$pop)
# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("bottomleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)
                         
# Plot the second and third principal components
plot(x = pca$scores[, 2], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("bottomleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

