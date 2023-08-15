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
rbrush_ref <- read.vcfR("rbrush/erna/rbrush_referna2.filtered.vcf")

#Convert to genlight format
focal_gl <- vcfR2genlight(rbrush_ref)
#populations map file
popmap<-read.csv("rbrush/erna/popmap2.csv")
popmap<-read.csv("rbrush/erna/popmap_dataset.csv")
popmap<-data.frame(id=colnames(rbrush_ref@gt)[2:length(colnames(rbrush_ref@gt))],pop=substr(colnames(rbrush_ref@gt)[2:length(colnames(rbrush_ref@gt))], 1,8))
#Add population info, setting up the population array
pop(focal_gl) <- c("arenaria", "arenaria", "arenaria", "arenaria", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "salicifolia", "iridis", "iridis", "iridis", "iridis", "iridis", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "graveolens", "graveolens", "graveolens", "graveolens", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "oreophila", "oreophila", "oreophila", "oreophila", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "juncea", "juncea", "juncea", "juncea", "juncea", "juncea", "juncea", "juncea", "juncea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "nitida", "nitida", "nitida", "nitida", "nitida", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "salicifolia", "salicifolia", "graveolens", "graveolens", "graveolens", "turbinata", "turbinata", "turbinata", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "turbinata", "turbinata", "turbinata", "unknown", "unknown", "salicifolia", "salicifolia", "unknown", "unknown", "unknown", "unknown", "unknown", "unknown", "speciosa", "speciosa", "speciosa", "speciosa", "turbinata", "turbinata", "turbinata", "turbinata", "turbinata", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "speciosa", "speciosa", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "oreophila", "oreophila", "oreophila", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "graveolens", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "oreophila", "desert", "desert", "desert", "desert", "desert", "desert", "speciosa", "speciosa", "speciosa", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "ammophila", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "bigelovii", "oreophila", "oreophila", "oreophila", "oreophila", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "hololeuca", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "latisquamea", "oreophila", "oreophila", "oreophila")
pop(focal_gl) <- c("2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", 
"1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", 
"1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1")

# Create a data.frame with sample names and second identifier
sample_info <- data.frame(sample_names = popmap$id, second_id = popmap$ds)

# Assign the data.frame as an attribute to the genlight object
attr(focal_gl, "sample_info") <- sample_info

# Run PCA
focal_pca <- glPca(focal_gl, n.cores = 4, nf = 4)

# Plot the PCA with differentiation based on population and second identifier
s.class(focal_pca$scores[, c(1, 2)], pop(focal_gl),
        col = magma(15, begin = 0.8, end = 0), clab = , cell = 2.5,
        pch = attr(focal_gl, "sample_info")$second_id)

#barplot plots the proportion of the total variance explained by given PCs
#Generally, the more a PC sticks out the more it matters
barplot(focal_pca$eig/sum(focal_pca$eig), main="eigenvalues", 
        col=heat.colors(length(focal_pca$eig)))
#plot PCA with base
pca <- glPca(focal_gl, nf = 30)
# Quick plot
scatter(pca, posi = "none")
# Create named color vector for the legend
legend_colors <- hcl.colors(15, palette = "spectral")
names(legend_colors) <- c("ammophila", "arenaria", "bigelovii", "desert", "graveolens", "hololeuca", "iridis", "juncea", "latisquamea", "nitida", "oreophila", "salicifolia", "speciosa", "turbinata", "unknown")     
#levels:ammophila arenaria bigelovii desert graveolens hololeuca iridis juncea latisquamea nitida oreophila salicifolia speciosa turbinata unknown

# Create named color vector for the legend for the dataset only
#legend_colors <- c("#009E73", "#0072B2")
#"#E69F00", "#D55E00", "#CC79A7"
#names(legend_colors) <- c("1", "2")     
#levels:1 2

popmap$pop=as.factor(popmap$pop)

# Define a vector of filled symbols
filled_symbols <- c(15, 16)

# Plot the first two principal components with filled symbols
plot(x = focal_pca$scores[, 1], y = focal_pca$scores[, 2],
     col = legend_colors[as.numeric(pop(focal_gl))],
     cex = 1, pch = filled_symbols[attr(focal_gl, "sample_info")$second_id])

# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)
# Add a legend for the filled symbols
legend("topright", legend = levels(attr(focal_gl, "sample_info")$second_id),
       pch = filled_symbols, title = "Second Identifier", bty = "n",
       ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

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

####Nolan's way####
Con.Source.plot <- ggplot() +
  ylab(expression(paste(delta^{2}, "H (\u2030)")))+
  xlab(expression(paste(delta^{18}, "O (\u2030)"))) +
  theme(text = element_text(size=15)) +
  scale_colour_manual(values = col) +
  theme_linedraw(base_size = 22) +
  geom_point(data = Sum_s.data, aes(x = mean.O, y = mean.H, color = Common.Names), size = 10, shape = 18) +
  geom_errorbar(data = Sum_s.data, aes(x = mean.O, y = mean.H,
                                       ymin = mean.H - sd.H,
                                       ymax = mean.H + sd.H, color = Common.Names)) +
  geom_errorbar(data = Sum_s.data, aes(x = mean.O, y = mean.H,
                                       xmin = mean.O - sd.O,
                                       xmax = mean.O + sd.O, color = Common.Names)) +
  geom_point(data = Cotton.Wood.Data.c, aes(x= d18O, y = d2H, color = Species.and.Size,
                                            fill = Species.and.Size, size = 2)) +
  stat_ellipse(data = Cotton.Wood.Data.c, aes(x= d18O, y = d2H,color = Species.and.Size),
               alpha = 0.0  ,
               type = "norm",
               geom = "polygon",
               size = 1)
Con.Source.plot
