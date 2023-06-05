####READ ME####
#this is the pca set up for aciurina mfc data set
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
rbrush_ref <- read.vcfR("rbrush/reference_based/rbrush_ref90.75.filtered.vcf")
rbrush_ref <- read.vcfR("rbrush_ref60.85.filtered.vcf.gz")
#Convert to genlight format
focal_gl <- vcfR2genlight(rbrush_ref)
#populations map file
popmap<-data.frame(id=colnames(rbrush_ref@gt)[2:length(colnames(rbrush_ref@gt))],pop=substr(colnames(rbrush_ref@gt)[2:length(colnames(rbrush_ref@gt))], 7, 8))
#popmap<-read.csv("rbrush/reference_based/popmap.csv")
#Add population info, setting up the population array
pop(focal_gl) <- c("_S", "_S", "_S", "_C", "_C", "_C", "_S", "_S", "_S", "_F", "_F", "_F", "_N", "_N", "_N", "_N", "_N", "_N", "_N", "_F", "_F", "_F", "_F", "_F", "_S", "_S", "_S", "_S", "_C", "_C", "_C", "_C", "_C", "_C", "_M", "_M", "_M", "_M", "_M", "_M", "_M", "_N", "_N", "_N", "_N", "_F", "_F", "_F", "_F", "_C", "_C", "_C", "_C", "_C", "_C", "_C", "_G", "_G", "_G", "_G", "_G", "_G", "_S", "_S", "_S", "_S", "_S", "_S", "_X", "_X", "_X", "_X", "_C", "_C", "_C", "_H", "_H", "_H", "_H", "_S", "_S", "_S", "_S", "_S", "_S", "_C", "_C", "_C", "_C", "_C", "_C", "_S", "_S", "_S", "_S", "_S", "_S", "_F", "_M", "_M", "_M", "_M", "_M", "_M", "_M", "_O", "_O", "_O", "_O", "_M", "_M", "_M", "_M", "_M", "_M", "_C", "_C", "_C", "_C", "_C", "_C", "_C", "_M", "_M", "_M", "_M", "_M", "_M", "_N", "_N", "_N", "_N", "_N", "_N", "_S", "_S", "_S", "_S", "_S", "_F", "_F", "_F", "_F", "_F", "_F", "_F", "_F", "_F", "_F", "_F", "_F", "_F", "_F")
pop(focal_gl) <- c("S", "S", "S", "C", "C", "S", "S", "F", "F", "N", "N", "N", "N", "N", "N", "F", "F", "F", "S", "S", "S", "C", "C", "C", "C", "M", "M", "M", "M", "M", "M", "M", "N", "N", "N", "F", "F", "C", "C", "C", "G", "G", "S", "S", "S", "S", "S", "X", "X", "C", "C", "H", "H", "H", "S", "S", "S", "S", "S", "C", "C", "C", "C", "C", "S", "S", "S", "S", "M", "M", "M", "M", "M", "M", "M", "O", "M", "M", "M", "M", "C", "C", "C", "C", "C", "M", "M", "M", "M", "M", "N", "N", "N", "N", "N", "N", "S", "S", "S", "S", "S", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F")
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
#legend_colors <- rainbow(10)

legend_colors <- c("turquoise", "black","seagreen","gold", "lightgreen",
                              "purple","coral", "hotpink" , "blue")
                              
legend_colors <- c("#A0CBE8", "#F28E2B","#8CD17D","#F1CE63", "#499894",
                              "#D37295","#BAB0AC", "#D4A6C8" , "#D7B5A6")
                            
legend_colors <- c("#4E79A7", "#F28E2B","#59A14F","#B6992D", "#499894",
                            "#79706E","#D37295", "#B07AA1" , "#D7B5A6")
#paletteer_d("ggthemes::Tableau_x")

legend_colors <- paletteer_d("ggthemes::Miller_Stone")
legend_colors <- paletteer_d("ggthemes::few_Medium")

names(legend_colors) <- c("C", "F", "G", "H", "M", "N", "O", "S", "X")     
#levels: _C _F _G _H _M _N _O _S _X

popmap$pop=as.factor(popmap$pop)

# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("bottomright", legend = names(legend_colors), fill = legend_colors,
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
