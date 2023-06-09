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
focal_vcf <- read.vcfR("rbrush/vcf_all/rbrush.filtered.vcf")
#Convert to genlight format
focal_gl <- vcfR2genlight(focal_vcf)
#Add population info, setting up the population array
pop(focal_gl) <- c("S", "S", "S", "C", "C", "C", "S", "S", "S", "F", "F", "F", "N", "N", "N", "N", "N", "N", "N", "F", "F", "F", "F", "F", "S", "S", "S", "S", "C", "C", "C", "C", "C", "C", "M", "M", "M", "M", "M", "M", "M", "N", "N", "N", "N", "F", "F", "F", "F", "C", "C", "C", "C", "C", "C", "C", "G", "G", "G", "G", "G", "G", "S", "S", "S", "S", "S", "S", "X", "X", "X", "X", "C", "C", "C", "H", "H", "H", "H", "S", "S", "S", "S", "S", "S", "C", "C", "C", "C", "C", "C", "S", "S", "S", "S", "S", "S", "F", "M", "M", "M", "M", "M", "M", "M", "O", "O", "O", "O", "M", "M", "M", "M", "M", "M", "C", "C", "C", "C", "C", "C", "C", "M", "M", "M", "M", "M", "M", "N", "N", "N", "N", "N", "N", "S", "S", "S", "S", "S", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F")
#Run PCA
focal_pca <- glPca(focal_gl, n.cores=4, nf=4)
#Not really great, but can be helpful for visualizing potentially problematic samples
scatter(focal_pca, cex=.25)
#My preferred option, I generally make plots here and edit in illustrator
s.class(focal_pca$scores[,c(1,3)], pop(focal_gl), 
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
                            
legend_colors <- c("#0079bf", "#d29034","#838c91","#b04632", "#89609e",
                            "#cd5a91","#1e663b", "#8CD17D" , "#D7B5A6")
                            #paletteer_d("ggthemes::Tableau_x")

names(legend_colors) <- c("C", "F", "G", "H", "M", "N", "O", "S", "X")     
#levels: C F G H M N O S X

popmap$pop=as.factor(popmap$pop)

# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("topright", legend = names(legend_colors), fill = legend_colors,
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
