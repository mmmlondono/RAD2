####READ ME####
#this is the pca set up for rabbit brush 
####load packages ####
library("adegenet")
library("ade4")
library("vcfR")
library("parallel")
library("viridis")
####adegenet’s glPCA####
#This is run in Rstudio (it being an IDE is especially nice for visualizing)
#Read in VCF
focal_vcf <- read.vcfR("rbrush/vcf/aciurinaR_.75.filtered.vcf")
#Convert to genlight format
focal_gl <- vcfR2genlight(focal_vcf)
#Add population info, setting up the population array
pop(focal_gl) <- c("ABQ_S","ABQ_S","ABQ_S", 
                   "MAL_C","MAL_C","MAL_C",
                   "DRY_N","DRY_N",
                   "MAL_M","MAL_M","MAL_M",
                   "RGG_C","RGG_C","RGG_C",
                   "RGG_S","RGG_S","RGG_S",
                   "GUN_F","GUN_F")
#Run PCA
focal_pca <- glPca(focal_gl, n.cores=4, nf=4)
#Not really great, but can be helpful for visualizing potentially problematic samples
scatter(focal_pca, cex=.25)
#My preferred option, I generally make plots here and edit in illustrator
s.class(focal_pca$scores[,c(1,3)], pop(focal_gl), 
        #Col is the color pallete used, magma is a nice default
        #Cell is effectively a confidence interval (2.5 ~= 95%)
        col=magma(7, begin=.8, end=0), clab=.6)
#barplot plots the proportion of the total variance explained by given PCs
#Generally, the more a PC sticks out the more it matters
barplot(focal_pca$eig/sum(focal_pca$eig), main="eigenvalues", 
        col=heat.colors(length(focal_pca$eig)))


pop <- c("ABQ_S","ABQ_S","ABQ_S", 
         "MAL_C","MAL_C","MAL_C",
         "DRY_N","DRY_N",
         "MAL_M","MAL_M","MAL_M",
         "RGG_C","RGG_C","RGG_C",
         "RGG_S","RGG_S","RGG_S",
         "GUN_F","GUN_F")

####other method####

## Heat map of genotype
glPlot(focal_gl)

pca <- glPca(focal_gl, nf = 10)
# Quick plot
scatter(pca, posi = "none")
# Create named color vector for the legend
legend_colors <- c("blue", "green", "darkgreen",
                         "orange", "purple", "red", "black")
names(legend_colors) <- c("ABQ_S","MAL_C","DRY_N", "MAL_M","RGG_C", "RGG_S","GUN_F")
                         
# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = colors, cex = 1, pch = 19)
 # Add a legend
legend("topright", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white")

# Plot the second and third principal components
plot(x = pca$scores[, 2], y = pca$scores[, 3], col = colors, cex = 1, pch = 19)
# Add a legend
legend("topright", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white")
                         
 





