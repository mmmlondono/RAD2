
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
aciu_vcf <- read.vcfR("wolbachia/marshmallow.filtered.vcf")

#Convert to genlight format
focal_gl <- vcfR2genlight(aciu_vcf)
#populations map file
popmap<-data.frame(id=colnames(aciu_vcf@gt)[2:length(colnames(aciu_vcf@gt))],pop=substr(colnames(aciu_vcf@gt)[2:length(colnames(aciu_vcf@gt))], 1,8))
#Add population info, setting up the population array
pop(focal_gl) <- c("ABQ_FY_S", "ABQ_FY_S", "ABQ_FY_S", "MAL_FY_C", "MAL_FY_C", "MAL_FY_C", "MAL_FY_N", "MAL_FY_N", "MAL_FY_N", "RGG_FY_C", "RGG_FY_C", "RGG_FY_C", "RGG_FY_S", "RGG_FY_S", "RGG_FY_S", "SNM_FY_C", "SNM_FY_C", "SNM_FY_C", "SNM_FY_S", "SNM_FY_S", "SNM_FY_S", "TES_FY_C", "TES_FY_C", "TES_FY_C", "TES_FY_N", "TES_FY_N", "TES_FY_N", "TES_FY_S", "TES_FY_S", "TES_FY_S")
pop(focal_gl) <- c("ABQ_S", "ABQ_S", "ABQ_S", "MAL_C", "MAL_C", "MAL_C", "MAL_N", "MAL_N", "MAL_N", "RGG_C", "RGG_C", "RGG_C", "RGG_S", "RGG_S", "RGG_S", "SLV_O", "SLV_O", "SLV_O", "SNM_C", "SNM_C", "SNM_C", "SNM_S", "SNM_S", "SNM_S", "TES_C", "TES_C", "TES_C", "TES_N", "TES_N", "TES_N", "TES_S", "TES_S", "TES_S", "MAL_M", "MAL_M", "MAL_M", "SLV_M", "SLV_M", "SLV_M", "SNM_M", "SNM_M", "SNM_M", "TES_M", "TES_M", "TES_M")

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

legend_colors <- c("#009E73", "#65C897", "#A3E0B8","#C2F2D9","#0072B2", "#4C8CCD","#6FA7D1","#99B5E8","#E69F00", "#F4C266")
names(legend_colors) <- c("SNM-S","ABQ-S","TES-S","RGG-S","SNM-C","MAL-C","TES-C","RGG-C","MAL-N","TES-N")     
#levels: ABQ-S MAL-C MAL-N RGG-C RGG-S SNM-C SNM-S TES-C TES-N TES-S

popmap$pop=as.factor(popmap$pop)

# Calculate proportion of variance explained
prop_var <- pca$eig / sum(pca$eig) * 100
# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors[popfocal$pop],cex = 1.5, pch = 19,
     xlab=c("0", paste0("PC1 (", round(prop_var[1], 2), "%)")), ylab=c("0", paste0("PC2 (", round(prop_var[2], 2), "%)")))

# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
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




