
pca_scores <- pca$scores
write.csv(pca_scores, "pca_scores.csv")

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
aciu_vcf <- read.vcfR("wolbachia/wolbachia.filtered.vcf")

#Convert to genlight format
focal_gl <- vcfR2genlight(aciu_vcf)
#populations map file
popmap<-data.frame(id=colnames(aciu_vcf@gt)[2:length(colnames(aciu_vcf@gt))],pop=substr(colnames(aciu_vcf@gt)[2:length(colnames(aciu_vcf@gt))], 1,8))
#Add population info, setting up the population array
pop(focal_gl) <- c("ABQ_FY_S", "ABQ_FY_S", "ABQ_FY_S", "MAL_FY_C", "MAL_FY_C", "MAL_FY_C", "MAL_FY_N", "MAL_FY_N", "MAL_FY_N", "RGG_FY_C", "RGG_FY_C", "RGG_FY_C", "RGG_FY_S", "RGG_FY_S", "RGG_FY_S", "SNM_FY_C", "SNM_FY_C", "SNM_FY_C", "SNM_FY_S", "SNM_FY_S", "SNM_FY_S", "TES_FY_C", "TES_FY_C", "TES_FY_C", "TES_FY_N", "TES_FY_N", "TES_FY_N", "TES_FY_S", "TES_FY_S", "TES_FY_S")

popfocal<- read.csv("wolbachia/popmap.csv", header=TRUE)
pop(focal_gl) <- popfocal$pop

#run pca#
pca <- glPca(focal_gl, nf = 30)

#barplot plots the proportion of the total variance explained by given PCs
#Generally, the more a PC sticks out the more it matters
barplot(pca$eig/sum(pca$eig), main="eigenvalues", 
        col=heat.colors(length(pca$eig)))

# Quick plot
scatter(pca, posi = "none")
# Create named color vector for the legend

legend_colors <- c("#009E73", "#E69F00", "#56B4E9", "#999999", "deeppink4",
                            "#F0E442", "#0072B2", "#D55E00", "#000000", "#CC79A7")
                            
names(legend_colors) <- c("ABQ-S", "MAL-C", "MAL-N", "RGG-C", "RGG-S", "SNM-C", "SNM-S", "TES-C", "TES-N", "TES-S")     
#levels: 

popmap$pop=as.factor(popmap$pop)

# Calculate proportion of variance explained
prop_var <- pca$eig / sum(pca$eig) * 100
# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors[popfocal$pop],cex = 1, pch = 19,
     xlab=c("0", paste0("PC1 (", round(prop_var[1], 2), "%)")), ylab=c("0", paste0("PC2 (", round(prop_var[2], 2), "%)")))

# Add a legend
legend("topright", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.7)

# Plot the first and third principal components
plot(x = pca$scores[, 1], y = pca$scores[, 3], col = legend_colors[popfocal$pop], cex = 1, pch = 19,
     xlab=c("0", paste0("PC1 (", round(prop_var[1], 2), "%)")), ylab=c("0", paste0("PC3 (", round(prop_var[3], 2), "%)")))

# Add a legend
legend("topright", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the second and third principal components
plot(x = pca$scores[, 2], y = pca$scores[, 3], col = legend_colors[popfocal$pop], cex = 1, pch = 19,       
     xlab=c("0", paste0("PC2 (", round(prop_var[2], 2), "%)")), ylab=c("0", paste0("PC3 (", round(prop_var[3], 2), "%)")))
# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

pca_scores <- pca$scores
write.csv(pca_scores, "pca_scores_aciurina.csv")




