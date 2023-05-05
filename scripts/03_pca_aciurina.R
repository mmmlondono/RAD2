####READ ME####
#this is the pca set up for aciurina
####load packages ####
library("adegenet")
library("ade4")
library("vcfR")
library("parallel")
library("viridis")
####adegenet’s glPCA####
#This is run in Rstudio (it being an IDE is especially nice for visualizing)
#Read in VCF
focal_vcf <- read.vcfR("aciurina/aciurina2_.6.filtered.vcf")
focal_vcf <- read.vcfR("aciurina/aciurina2_.8.filtered.vcf")
focal_vcf <- read.vcfR("aciurina/aciurina4_.8.filtered.vcf")

#Convert to genlight format
focal_gl <- vcfR2genlight(focal_vcf)
#Add population info, setting up the population array
pop(focal_gl) <- c("ABQ-S", "ABQ-S", "ABQ-S", 
                   "BOU-C", "BOU-C", "BOU-C", 
                   "BOU-S", "BOU-S", "BOU-S", 
                   "DRY-N", "DRY-N", "DRY-N", 
                   "GUN-F", "GUN-F", "GUN-F", 
                   "LOG-S", "LOG-S", "LOG-S", 
                   "MAL-C", "MAL-C", "MAL-C", 
                   "MAL-M", "MAL-M", "MAL-M", 
                   "MAL-N", "MAL-N", "MAL-N", 
                   "MOB-C", "MOB-C", "MOB-C", 
                   "MOB-G", "MOB-G", "MOB-G", 
                   "MOB-S", "MOB-S", "MOB-S", 
                   "MTB-F", "MTB-F", "MTB-F", 
                   "PAS-C", "PAS-C", "PAS-C", 
                   "REN-C", "REN-C", "REN-C", 
                   "REN-HC", "REN-HC", "REN-HC", 
                   "REN-HS", "REN-HS", "REN-HS", 
                   "REN-S", "REN-S", "REN-S", 
                   "RGG-C", "RGG-C", "RGG-C", 
                   "RGG-S", "RGG-S", "RGG-S", 
                   "SLV-F", "SLV-F", "SLV-F", 
                   "SLV-M", "SLV-M", "SLV-M", 
                   "SLV-O", "SLV-O", "SLV-O", 
                   "SNM-C", "SNM-C", "SNM-C", 
                   "SNM-M", "SNM-M", "SNM-M", 
                   "SNM-S", "SNM-S", "SNM-S", 
                   "TES-C", "TES-C", "TES-C", 
                   "TES-M", "TES-M", "TES-M", 
                   "TES-N", "TES-N", "TES-N", 
                   "TES-S", "TES-S", "TES-S", 
                   "TRE-F", "TRE-F", "TRE-F", 
                   "TRQ-F", "TRQ-F", "TRQ-F")
#Run PCA
focal_pca <- glPca(focal_gl, n.cores=4, nf=4)
#Not really great, but can be helpful for visualizing potentially problematic samples
scatter(focal_pca, cex=.25)
#My preferred option, I generally make plots here and edit in illustrator
s.class(focal_pca$scores[,c(1,2)], pop(focal_gl), 
        col=magma(9, begin=.8, end=0), clab=1, cell=2.5)
        #col is the color pallet used, magma is a nice default
        #cell is effectively a confidence interval (2.5 ~= 95%)
        
#barplot plots the proportion of the total variance explained by given PCs
#Generally, the more a PC sticks out the more it matters
barplot(focal_pca$eig/sum(focal_pca$eig), main="eigenvalues", 
        col=heat.colors(length(focal_pca$eig)))


## Heat map of genotype
glPlot(focal_gl)

pca <- glPca(focal_gl, nf = 10)
# Quick plot
scatter(pca, posi = "none")
# Create named color vector for the legend
legend_colors <- rainbow(32)

legend_colors <- c("blue", "darkblue","royalblue", "cyan","aquamarine", "turquoise", "navy","lightblue","skyblue", 
                         "green", "darkgreen","limegreen", "seagreen","lightgreen", "khaki","gold","grey","slategray", 
                          "violet", "pink","hotpink","magenta", "purple",
                         "orange", "coral", "orangered" ,
                         "tan","brown", "sienna", "maroon", 
                         "black",
                         "red")
names(legend_colors) <- c("ABQ-S", "BOU-S", "MOB-S","LOG-S", "REN-HS", "REN-S", "RGG-S","SNM-S", "TES-S", 
                          "BOU-C", "MAL-C", "MOB-C", "PAS-C", "REN-C", "REN-HC","RGG-C","SNM-C","TES-C",
                          "GUN-F", "MTB-F",  "SLV-F", "TRE-F", "TRQ-F",
                          "DRY-N", "MAL-N", "TES-N",
                          "MAL-M", "SLV-M","SNM-M","TES-M",  
                          "SLV-O", 
                          "MOB-G")
                         
# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors, cex = 1, pch = 19)
# Add a legend
legend("bottomright", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.5)
                         
# Plot the second and third principal components
plot(x = pca$scores[, 2], y = pca$scores[, 3], col = colors, cex = 1, pch = 19)
# Add a legend
legend("topright", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white")
                         
                         
                         
                         
                         
                         