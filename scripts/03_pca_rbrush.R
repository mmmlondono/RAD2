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
focal_vcf <- read.vcfR("rbrush/vcf/populations.snps.vcf")
#Convert to genlight format
focal_gl <- vcfR2genlight(focal_vcf)
#Add population info, setting up the population array
pop(focal_gl) <- c("ABQ-S", "ABQ-S", "ABQ-S", "ABQ-S", "ABQ-S", "ABQ-S", "ABQ-S", 
                   "BOU-C", "BOU-C", "BOU-C", "BOU-S", "BOU-S", "BOU-S", 
                   "BRG-F", "BRG-F", "BRG-F", 
                   "DRY-N", "DRY-N", "DRY-N", "DRY-N", "DRY-N", "DRY-N", "DRY-N", 
                   "GUN-F", "GUN-F", "GUN-F", "GUN-F", "GUN-F", "GUN-F", "GUN-F", 
                   "LIM-F", "LIM-F", "LIM-F", "LIM-F", 
                   "LOG-S", "LOG-S", "LOG-S", "LOG-S", "LOG-S", "LOG-S", "LOG-S", 
                   "MAL-C", "MAL-C", "MAL-C", "MAL-C", "MAL-C", "MAL-C", "MAL-C",
                   "MAL-M", "MAL-M", "MAL-M", "MAL-M", "MAL-M", "MAL-M", "MAL-M", 
                   "MAL-N", "MAL-N", "MAL-N", "MAL-N", "MAL-N", "MAL-N", "MAL-N", 
                   "MEE-F", "MEE-F", "MEE-F", "MEE-F", "MEE-F", "MEE-F", "MEE-F", 
                   "MOB-C", "MOB-C", "MOB-C", "MOB-C", "MOB-C", "MOB-C", "MOB-C", 
                   "MOB-G", "MOB-G", "MOB-G", "MOB-G", "MOB-G", "MOB-G", 
                   "MOB-S", "MOB-S", "MOB-S", "MOB-S", "MOB-S", "MOB-S", "MOB-S", 
                   "PAS-F", "PAS-F", "PAS-F", "PAS-F", "PAS-F", "PAS-F", "PAS-F", 
                   "REN-C", "REN-C", "REN-C", 
                   "REN-H", "REN-H", "REN-H", "REN-H", "REN-H", "REN-H", "REN-H", 
                   "REN-S", "REN-S", "REN-S", "REN-S", "REN-S", "REN-S", "REN-S", 
                   "RGG-C", "RGG-C", "RGG-C", "RGG-C", "RGG-C", "RGG-C", "RGG-C", 
                   "RGG-S", "RGG-S", "RGG-S", "RGG-S", "RGG-S", "RGG-S", "RGG-S", 
                   "SLV-F", "SLV-F", "SLV-F", "SLV-F", "SLV-F", "SLV-F", "SLV-F", 
                   "SLV-M", "SLV-M", "SLV-M", "SLV-M", "SLV-M", "SLV-M", "SLV-M", 
                   "SLV-O", "SLV-O", "SLV-O", "SLV-O", "SLV-O", "SLV-O", "SLV-O", 
                   "SNM-M", "SNM-M", "SNM-M", "SNM-M", "SNM-M", "SNM-M", "SNM-M", 
                   "TES-C", "TES-C", "TES-C", "TES-C", "TES-C", "TES-C", "TES-C", 
                   "TES-M", "TES-M", "TES-M", "TES-M", "TES-M", "TES-M", "TES-M", 
                   "TES-N", "TES-N", "TES-N", "TES-N", "TES-N", "TES-N", "TES-N", 
                   "TES-S", "TES-S", "TES-S", "TES-S", "TES-S", "TES-S", "TES-S", 
                   "TRE-F", "TRE-F", "TRE-F", "TRE-F", "TRE-F", "TRE-F", "TRE-F", 
                   "TRQ-F", "TRQ-F", "TRQ-F", "TRQ-F", "TRQ-F", "TRQ-F", "TRQ-F")
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


pop <- c("ABQ-S", "ABQ-S", "ABQ-S", "ABQ-S", "ABQ-S", "ABQ-S", "ABQ-S", 
                   "BOU-C", "BOU-C", "BOU-C", "BOU-S", "BOU-S", "BOU-S", 
                   "BRG-F", "BRG-F", "BRG-F", 
                   "DRY-N", "DRY-N", "DRY-N", "DRY-N", "DRY-N", "DRY-N", "DRY-N", 
                   "GUN-F", "GUN-F", "GUN-F", "GUN-F", "GUN-F", "GUN-F", "GUN-F", 
                   "LIM-F", "LIM-F", "LIM-F", "LIM-F", 
                   "LOG-S", "LOG-S", "LOG-S", "LOG-S", "LOG-S", "LOG-S", "LOG-S", 
                   "MAL-C", "MAL-C", "MAL-C", "MAL-C", "MAL-C", "MAL-C", "MAL-C",
                   "MAL-M", "MAL-M", "MAL-M", "MAL-M", "MAL-M", "MAL-M", "MAL-M", 
                   "MAL-N", "MAL-N", "MAL-N", "MAL-N", "MAL-N", "MAL-N", "MAL-N", 
                   "MEE-F", "MEE-F", "MEE-F", "MEE-F", "MEE-F", "MEE-F", "MEE-F", 
                   "MOB-C", "MOB-C", "MOB-C", "MOB-C", "MOB-C", "MOB-C", "MOB-C", 
                   "MOB-G", "MOB-G", "MOB-G", "MOB-G", "MOB-G", "MOB-G", 
                   "MOB-S", "MOB-S", "MOB-S", "MOB-S", "MOB-S", "MOB-S", "MOB-S", 
                   "PAS-F", "PAS-F", "PAS-F", "PAS-F", "PAS-F", "PAS-F", "PAS-F", 
                   "REN-C", "REN-C", "REN-C", 
                   "REN-H", "REN-H", "REN-H", "REN-H", "REN-H", "REN-H", "REN-H", 
                   "REN-S", "REN-S", "REN-S", "REN-S", "REN-S", "REN-S", "REN-S", 
                   "RGG-C", "RGG-C", "RGG-C", "RGG-C", "RGG-C", "RGG-C", "RGG-C", 
                   "RGG-S", "RGG-S", "RGG-S", "RGG-S", "RGG-S", "RGG-S", "RGG-S", 
                   "SLV-F", "SLV-F", "SLV-F", "SLV-F", "SLV-F", "SLV-F", "SLV-F", 
                   "SLV-M", "SLV-M", "SLV-M", "SLV-M", "SLV-M", "SLV-M", "SLV-M", 
                   "SLV-O", "SLV-O", "SLV-O", "SLV-O", "SLV-O", "SLV-O", "SLV-O", 
                   "SNM-M", "SNM-M", "SNM-M", "SNM-M", "SNM-M", "SNM-M", "SNM-M", 
                   "TES-C", "TES-C", "TES-C", "TES-C", "TES-C", "TES-C", "TES-C", 
                   "TES-M", "TES-M", "TES-M", "TES-M", "TES-M", "TES-M", "TES-M", 
                   "TES-N", "TES-N", "TES-N", "TES-N", "TES-N", "TES-N", "TES-N", 
                   "TES-S", "TES-S", "TES-S", "TES-S", "TES-S", "TES-S", "TES-S", 
                   "TRE-F", "TRE-F", "TRE-F", "TRE-F", "TRE-F", "TRE-F", "TRE-F", 
                   "TRQ-F", "TRQ-F", "TRQ-F", "TRQ-F", "TRQ-F", "TRQ-F", "TRQ-F")
