#all following from Danielle Land
#except I didn't run fst distances here, I did those in CARC and output the table to a csv
#R Studio Fst Heatmap (Post SNPFiltR, on the same script)

library(adegenet)
library(dplyr)
library(hierfstat)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggplot2)

vcfR<-read.vcfR("wolbachia/wolbachia.filtered.vcf")
popmap <- read.csv("wolbachia/popmap.csv")
#Look at what we have
snps <- vcfR2genind(vcfR)
snps
popmap


x <- popmap$pop
x
x <- as.data.frame(x)
x

strata(snps) <- x
snps
snps@strata
setPop(snps) <- ~x
snps
snps@pop
popNames(snps)


### Identify the sub group that I care about

#snps.sub <- popsub(snps, exclude=c("Somon","Sorcam")) # to exclude specific populations
#snps.sub <- popsub(snps, sublist=c("Sorcae","Sorgra","Soriso", "Sormin",  "Sorspp","Sorung","Soryuk"))
#snps.sub

### Fst
# Compute pairwise Fsts
#fst <- genet.dist(snps, method = "WC84") %>% round(digits = 3)
fst <- read.csv("wolbachia/fst.csv")

popNames(snps)

# Desired order of labels
desired_labels <- c(".", "ABQ-S", "MAL-C", "MAL-N", "RGG-C", "RGG-S", "SNM-C", "SNM-S", "TES-C", "TES-N", "TES-S")

colnames(fst) <- c(".", "ABQ-S","RGG-S","SNM-S","TES-S","MAL-C","RGG-C","SNM-C","TES-C","MAL-N","TES-N")
rownames(fst) <- c("ABQ-S","RGG-S","SNM-S","TES-S","MAL-C","RGG-C","SNM-C","TES-C","MAL-N","TES-N")

nrow(fst)
ncol(fst)

fst <- fst[,2:11]
#lab_order2 = c(",4",)
# Change order of rows and cols
fst.mat = as.matrix(fst)
#fst.mat1 = fst.mat[lab_order, ]
#fst.mat2 = fst.mat1[, lab_order]

# Create a data.frame
ind = which(upper.tri(fst.mat), arr.ind = TRUE)
fst.df = data.frame(Site1 = dimnames(fst.mat)[[2]][ind[,2]],
                    Site2 = dimnames(fst.mat)[[1]][ind[,1]],
                    Fst = fst.mat[ ind ])

# Keep the order of the levels in the data.frame for plotting 
fst.df$Site1 = factor(fst.df$Site1, levels = unique(fst.df$Site1))
fst.df$Site2 = factor(fst.df$Site2, levels = unique(fst.df$Site2))

# Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] = 0

# Print data.frame summary
fst.df %>% str

# Fst italic label
fst.label = expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
fst.df$Fst <- as.numeric(fst.df$Fst)
mid = max(fst.df$Fst) / 2

# Plot heatmap
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 3)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0, 0.25, 0.5, .75, 1.0))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold.italic"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 10))

        
                
        