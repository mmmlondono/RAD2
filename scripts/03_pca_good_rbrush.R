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
rbrush_ref <- read.vcfR("good_rbrush/rbrush_good.filtered.vcf")
rbrush_ref <- read.vcfR("good_rbrush/2goods/rbrush_2good.filtered.vcf") #2
rbrush_ref <- read.vcfR("good_rbrush/3goods/rbrush_3good.filtered.vcf") #3

#Convert to genlight format
focal_gl_rbr <- vcfR2genlight(rbrush_ref)
#populations map file
popmap<-data.frame(id=colnames(rbrush_ref@gt)[2:length(colnames(rbrush_ref@gt))],pop=substr(colnames(rbrush_ref@gt)[2:length(colnames(rbrush_ref@gt))], 7,8))
#Add population info, setting up the population array
pop(focal_gl_rbr) <- c("S", "S", "C", "C", "S", "S", "F", "F", "N", "N", "N", "N", "N", "N", "F", "F", "F", "S", "S", "S", "C", "C", "C", "C", "M", "M", "M", "M", "M", "M", "M", "N", "F", "F", "C", "C", "C", "S", "S", "S", "S", "S", "X", "X", "C", "C", "H", "H", "H", "S", "S", "S", "S", "S", "C", "C", "C", "C", "C", "S", "S", "S", "S", "M", "M", "M", "M", "M", "M", "M", "C", "C", "C", "C", "C", "M", "M", "M", "M", "M", "N", "N", "N", "N", "N", "N", "S", "S", "S", "S", "S", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F")
pop(focal_gl_rbr) <- c("S", "C", "S", "F", "F", "N", "N", "N", "N", "N", "N", "F", "F", "F", "S", "S", "C", "C", "C", "C", "M", "M", "M", "M", "M", "M", "M", "N", "F", "F", "C", "C", "C", "S", "S", "S", "S", "S", "F", "F", "C", "C", "H", "H", "H", "S", "S", "S", "S", "S", "C", "C", "C", "C", "C", "S", "S", "S", "S", "M", "M", "M", "M", "M", "C", "C", "C", "C", "C", "M", "M", "M", "M", "M", "N", "N", "N", "N", "N", "N", "S", "S", "S", "S", "S", "F", "F", "F", "F", "F", "F", "F", "F", "F", "G", "G", "G", "O", "O", "O")
pop(focal_gl_rbr) <- c("C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "H", "H", "H", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "F", "F")

#glPlot(focal_gl_rbr)
pca_rbr <- glPca(focal_gl_rbr, nf = 30)

#barplot plots the proportion of the total variance explained by given PCs
#Generally, the more a PC sticks out the more it matters
barplot(pca_rbr$eig/sum(pca_rbr$eig), main="eigenvalues", 
        col=heat.colors(length(pca_rbr$eig)))

# Quick plot
scatter(pca_rbr, posi = "none")
# Create named color vector for the legend
legend_colors <- spectral(7)
legend_colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                   "#0072B2", "#D55E00", "#CC79A7")
legend_colors <- spectral(9)
legend_colors <- c("#009E73", "#E69F00", "#56B4E9", "#999999",
                   "#0072B2", "#D55E00", "#000000","#CC79A7")
legend_colors <- spectral(6)
legend_colors <- c("#009E73", "#E69F00", "#56B4E9",
                   "#0072B2", "#D55E00", "#CC79A7")

names(legend_colors) <- c("C", "F", "H", "M", "N", "S", "X")     
#levels: C F H M N S X
names(legend_colors) <- c("C", "F", "G", "H", "M", "N", "O", "S")
#levels: C F G H M N O S X
names(legend_colors) <- c("C", "F", "H", "M", "N", "S")  
#levels: C F H M N S

popmap$pop=as.factor(popmap$pop)

# calculate the first two principal components for plotting
prop_var <- pca$eig / sum(pca$eig) * 100
# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors[as.numeric(popmap$pop)],cex = 1, pch = 19,
     xlab=c("0", paste0("PC1 (", round(prop_var[1], 2), "%)")), ylab=c("0", paste0("PC2 (", round(prop_var[2], 2), "%)")))

# Add a legend
legend("topright", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the second and third principal components
plot(x = pca$scores[, 2], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19,       
     xlab=c("0", paste0("PC2 (", round(prop_var[2], 2), "%)")), ylab=c("0", paste0("PC3 (", round(prop_var[3], 2), "%)")))

# Add a legend
legend("topright", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the first and third principal components
plot(x = pca$scores[, 1], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19,
     xlab=c("0", paste0("PC1 (", round(prop_var[1], 2), "%)")), ylab=c("0", paste0("PC3 (", round(prop_var[3], 2), "%)")))

# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

pca_scores <- pca$scores
write.csv(pca_scores, "pca_scores.csv")
