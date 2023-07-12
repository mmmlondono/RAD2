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
aciu_vcf <- read.vcfR("aciurina/all/aciurina2_.8.filtered.vcf")
aciu_vcf <- read.vcfR("aciurina/aciurina2/aciurina.filtered.vcf")

#Convert to genlight format
focal_gl <- vcfR2genlight(aciu_vcf)
#populations map file
popmap<-data.frame(id=colnames(aciu_vcf@gt)[2:length(colnames(aciu_vcf@gt))],pop=substr(colnames(aciu_vcf@gt)[2:length(colnames(aciu_vcf@gt))], 8,9))
#Add population info, setting up the population array
pop(focal_gl) <- c("S", "S", "S", "C", "C", "C", "S", "S", "S", "N", "N", "N", "F", "F", "F", "S", "S", "S", "C", "C", "C", "M", "M", "M", "N", "N", "N", "C", "C", "C", "G", "G", "G", "S", "S", "S", "F", "F", "F", "F", "F", "F", "C", "C", "C", "HC", "HC", "HC", "HS", "HS", "HS", "S", "S", "S", "C", "C", "C", "S", "S", "S", "F", "F", "F", "M", "M", "M", "O", "O", "O", "C", "C", "C", "M", "M", "M", "S", "S", "S", "C", "C", "C", "M", "M", "M", "N", "N", "N", "S", "S", "S", "F", "F", "F", "F", "F", "F")
pop(focal_gl) <- c("C_", "C_", "C_", "C_", "C_", "C_", "C_", "C_", "C_", "F_", "F_", "F_", "MC", "MC", "MC", "C_", "C_", "C_", "C_", "C_", "C_", "C_", "C_", "C_", "F_", "F_", "F_", "F_", "F_", "F_", "F_", "F_", "F_", "F_", "F_", "F_", "F_", "F_", "F_", "G_", "G_", "G_", "MC", "MC", "MC", "S_", "S_", "S_", "M_", "M_", "M_", "M_", "M_", "M_", "M_", "M_", "M_", "M_", "M_", "M_", "N_", "N_", "N_", "N_", "N_", "N_", "N_", "N_", "N_", "O_", "O_", "O_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_", "S_")
pop(focal_gl) <- c("C", "C", "C", "C", "C", "C", "C", "C", "C", "X", "X", "X", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "G", "G", "G", "HC", "HC", "HC", "HS", "HS", "HS", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "N", "N", "N", "N", "N", "N", "N", "N", "N", "O", "O", "O", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S")

popfocal<- read.csv("aciurina/aciurina2/popfocal.csv", header=FALSE)
pop(focal_gl) <- popfocal$V1

#run pca#
pca <- glPca(focal_gl, nf = 30)

#barplot plots the proportion of the total variance explained by given PCs
#Generally, the more a PC sticks out the more it matters
barplot(pca$eig/sum(pca$eig), main="eigenvalues", 
        col=heat.colors(length(pca$eig)))

# Quick plot
scatter(pca, posi = "none")
# Create named color vector for the legend

legend_colors <- c("#009E73", "#E69F00", "#56B4E9", "#999999", #"deeppink4",
                   "#F0E442", "#0072B2", "#D55E00", "#000000", "#CC79A7")

legend_colors <- c("#009E73", "#E69F00", "#56B4E9", "#0072B2",
                   "#999999","#D55E00", "#000000", "#CC79A7")

legend_colors <- c("#009E73", "#E69F00", "#56B4E9", "#999999", "deeppink4",
                   "#F0E442", "#0072B2", "#D55E00", "#000000", "#CC79A7")

names(legend_colors) <- c("C", "F", "G", "HC", "HS", "M", "N", "O", "S")     
#levels: C F G HC HS M N O S X
#C F G HC HS M N O S
names(legend_colors) <- c("C", "F", "G", "M", "MC", "N", "O", "S")     
#levels: C_ F_ G_ M_ MC N_ O_ S_
names(legend_colors) <- c("C", "F", "G", "HC", "HS", "M", "N", "O", "S", "X")     
#levels: C F G HC HS M N O S X

popmap$pop=as.factor(popmap$pop)

# Calculate proportion of variance explained
prop_var <- pca$eig / sum(pca$eig) * 100
# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors[popfocal$V1],cex = 1, pch = 19,
     xlab=c("0", paste0("PC1 (", round(prop_var[1], 2), "%)")), ylab=c("0", paste0("PC2 (", round(prop_var[2], 2), "%)")))

# Add a legend
legend("bottomleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the first and third principal components
plot(x = pca$scores[, 1], y = pca$scores[, 3], col = legend_colors[popfocal$V1], cex = 1, pch = 19,
     xlab=c("0", paste0("PC1 (", round(prop_var[1], 2), "%)")), ylab=c("0", paste0("PC3 (", round(prop_var[3], 2), "%)")))

# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the second and third principal components
plot(x = pca$scores[, 2], y = pca$scores[, 3], col = legend_colors[popfocal$V1], cex = 1, pch = 19,       
       xlab=c("0", paste0("PC2 (", round(prop_var[2], 2), "%)")), ylab=c("0", paste0("PC3 (", round(prop_var[3], 2), "%)")))
# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

pca_scores <- pca$scores
write.csv(pca_scores, "pca_scores_aciurina.csv")

                         
                         