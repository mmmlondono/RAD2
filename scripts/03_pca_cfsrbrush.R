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
focal_vcf <- read.vcfR("rbrush/cfs_rbrush/cfsrbrush.filtered.vcf")
#Convert to genlight format
focal_gl <- vcfR2genlight(focal_vcf)
#Add population info, setting up the population array
pop(focal_gl) <- c("F", "F", "C", "C", "F", "F", "F", "C", "F", "C", "S", "F", "C", "F", "S", "S", "F", "S", "S", "S", "S", "C", "F", "S", "S", "S", "F", "C", "S", "S", "C", "F", "C", "S", "S", "C", "S", "C", "C", "S", "C", "F", "C", "C", "S", "F", "S", "S", "C", "S", "F", "S", "S", "F", "C", "C", "C", "F", "S", "F", "C", "C", "F", "C", "S", "C", "C", "S", "F", "S", "F", "F", "C", "C", "C", "F", "F", "S", "F", "C")

#populations map file
popmap<-data.frame(id=colnames(focal_vcf@gt)[2:length(colnames(focal_vcf@gt))],pop=substr(colnames(focal_vcf@gt)[2:length(colnames(focal_vcf@gt))], 7,8))

#Run PCA
focal_pca <- glPca(focal_gl, n.cores=4, nf=4)
#Not really great, but can be helpful for visualizing potentially problematic samples
scatter(focal_pca, cex=.25)
#My preferred option, I generally make plots here and edit in illustrator
s.class(focal_pca$scores[,c(1,2)], pop(focal_gl), 
        #Col is the color pallete used, magma is a nice default
        #Cell is effectively a confidence interval (2.5 ~= 95%)
        col=magma(10, begin=.8, end=0), clab=1, cell=2.5)
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
#legend_colors <- rainbow(9)

legend_colors <- c("#0079bf", "#d29034", "#8CD17D")
                            #paletteer_d("ggthemes::Tableau_x")

names(legend_colors) <- c("C", "F", "S")     
#levels: C F S

popmap$pop=as.factor(popmap$pop)

# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("bottomleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the second and third principal components
plot(x = pca$scores[, 2], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the first and third principal components
plot(x = pca$scores[, 1], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("bottomleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

pca_scores <- pca$scores
write.csv(pca_scores, "pca_scores.csv")