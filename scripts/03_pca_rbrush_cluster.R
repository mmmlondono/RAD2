####READ ME####
#this is the pca set up 
####load packages ####
library("adegenet")
library("ade4")
library("vcfR")
library("parallel")
library("viridis")
####adegenet’s glPCA####
#This is run in Rstudio (it being an IDE is especially nice for visualizing)
#Read in VCF
rbrushcluster_vcf <- read.vcfR("rbrush/reference_based/cluster/rbrush_cluster.filtered.vcf")

#Convert to genlight format
focal_gl <- vcfR2genlight(rbrushcluster_vcf)
#populations map file
popmap<-data.frame(id=colnames(rbrushcluster_vcf@gt)[2:length(colnames(rbrushcluster_vcf@gt))],pop=substr(colnames(rbrushcluster_vcf@gt)[2:length(colnames(rbrushcluster_vcf@gt))], 1,3))

#Add population info, setting up the population array
pop(focal_gl) <- c("TES", "TES", "TES", "TES", "TES", "MOB", "TES", "TRE", "TES", "RGG", "RGG", "TES", "BOU", "TRQ", "RGG", "LOG", "SLV", "BOU", "TRE", "TES", "SLV", "MAL", "RGG", "TRQ", "MOB", "MAL", "ABQ", "MAL", "TRE", "TES", "RGG", "LOG", "ABQ", "MAL", "RGG", "TRQ", "BOU", "GUN", "BRG", "PAS", "TRQ", "RGG", "TRE", "GUN", "BOU", "RGG", "MAL", "MOB", "TRE", "SNM", "RGG", "PAS", "MOB", "RGG", "BOU", "MAL", "BRG", "LIM", "MOB", "MOB", "SNM", "SLV", "LOG", "TES", "MOB", "TRE", "MOB", "MOB", "LIM", "GUN", "MOB", "SLV", "MOB", "RGG", "TES", "MEE", "MOB", "TRQ", "PAS", "LOG", "BRG", "MAL", "MEE", "REN", "RGG", "BOU", "REN", "REN", "REN", "REN", "REN", "MOB", "MOB", "TES", "REN", "REN", "REN", "SLV", "MEE", "REN", "TRQ", "REN", "PAS", "MEE", "MOB", "REN", "TRE", "MOB", "MOB", "MOB", "REN")
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
legend_colors <- rainbow(17)

legend_colors <- c("#4E79A7", "#F28E2B","#59A14F","#B6992D", "#499894",
                            "#79706E","#D37295", "#B07AA1" , "#D7B5A6")
                            
legend_colors <- paletteer_d("ggthemes::few_Medium")
                             
names(legend_colors) <- c("C", "F", "G", "H", "M", "N", "O", "S", "X")
#C F G H M N O S X
names(legend_colors) <-c("ABQ", "BOU", "BRG", "GUN", "LIM", "LOG", "MAL", "MEE", "MOB", "PAS", "REN", "RGG", "SLV", "SNM", "TES", "TRE", "TRQ")
popmap$pop=as.factor(popmap$pop)
# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the second and third principal components
plot(x = pca$scores[, 2], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the first and third principal components
plot(x = pca$scores[, 1], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

pca_cluster <- pca$scores
write.csv(pca_cluster, "pca_cluster.csv")
