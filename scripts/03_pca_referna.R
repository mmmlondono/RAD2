####READ ME####
#this is the pca set up for erna data set
####load packages ####
library("adegenet")
library("ade4")
library("vcfR")
library("parallel")
library("viridis")
library("paletteer")
####adegenet’s glPCA####
#This is run in Rstudio (it being an IDE is especially nice for visualizing)
#Read in VCF
rbrush_ref <- read.vcfR("rbrush/erna/rbrush_referna.filtered.vcf")

#Convert to genlight format
focal_gl <- vcfR2genlight(rbrush_ref)
#populations map file
popmap<-data.frame(id=colnames(rbrush_ref@gt)[2:length(colnames(rbrush_ref@gt))],pop=substr(colnames(rbrush_ref@gt)[2:length(colnames(rbrush_ref@gt))], 7, 8))
popmap<-read.csv("rbrush/erna/popmap2.csv")
#Add population info, setting up the population array
pop(focal_gl) <- c("arenaria", "arenaria", "arenaria", "arenaria", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "iridis", "iridis", "iridis", "iridis", "iridis", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "graveolens", "graveolens", "graveolens", "graveolens", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "oreophila", "oreophila", "oreophila", "oreophila", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "juncea", "juncea", "juncea", "juncea", "juncea", "juncea", "juncea", "juncea", "juncea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "nitida", "nitida", "nitida", "nitida", "nitida", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "salicifolia", "salicifolia", "graveolens", "graveolens", "graveolens", "turbinata", "turbinata", "turbinata", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "turbinata", "turbinata", "turbinata", "unknown", "unknown", "salicifolia", "salicifolia", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "speciosa", "speciosa", "speciosa", "speciosa", "turbinata", "turbinata", "turbinata", "turbinata", "turbinata", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "speciosa", "speciosa", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "oreophila", "oreophila", "oreophila", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "desert", "desert", "desert", "desert", "desert", "desert", "speciosa", "speciosa", "speciosa", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "oreophila", "oreophila", "oreophila", "oreophila", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "oreophila", "oreophila", "oreophila")

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
#glPlot(focal_gl)

pca <- glPca(focal_gl, nf = 30)
# Quick plot
scatter(pca, posi = "none")
# Create named color vector for the legend

legend_colors <- hcl.colors(15, palette = "spectral")

names(legend_colors) <- c("ammophila", "arenaria", "bigelovii", "desert", "graveolens", "hololeuca", "iridis", "juncea", "latisquamea", "nitida", "oreophila", "salicifolia", "speciosa", "turbinata", "unknown")     
#levels:ammophila arenaria bigelovii desert graveolens hololeuca iridis juncea latisquamea nitida oreophila salicifolia speciosa turbinata unknown

popmap$pop=as.factor(popmap$pop)

# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the second and third principal components
plot(x = pca$scores[, 2], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("topright", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the first and third principal components
plot(x = pca$scores[, 1], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("topright", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

pca_scores <- pca$scores
write.csv(pca_scores, "pca_scores.csv")
