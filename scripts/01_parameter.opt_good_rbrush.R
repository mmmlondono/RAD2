####READ ME####
# The purpose of this script is to optimize the different de novo parameters for rbrsh RAD data a small subset of the samples
#Osuna-Mascaro:
#We allowed a maximum missing data of 20 % (--max-missing 0.8), 
#a minimum minor allele frequency of 0.03 (--maf 0.03) and 
#specified a thin value of 5 (--thin 5), which allows that no two sites are within the specified distance from one another. 
#Also, we only included sites with quality scores above 10 (--minQ 10).

##Mean depth
#Generates a file containing the mean depth per site averaged across all individuals.
#the number of reads that have mapped to this position (coverage)
#the more reads that cover a site, the higher confidence our basecall is. 
#Should we set a maximum and a minimum depth?
mean_depth

##Mean depth per individual
#Generates a file containing the mean depth per individual.
depth_indiv

##Variant missingness
#Generates a file reporting the missingness on a per-site basis.
#This is a measure of how many individuals lack a genotype at a call site.
#Typically missingness of 75-95% is used.
variant_missingness

##Minor alelle frecuencies
#Outputs the allele frequency for each site (distribution of allele frequencies)
#This will help inform our minor-allele frequency (MAF) thresholds
#Usually then, it is best practice to produce one dataset with a good MAF threshold and keep another without any MAF filtering at all.
#Ideally you want an idea of the distribution of your allelic frequencies but 0.05 to 0.10 is a reasonable cut-off.
minor_allele

##Proportion of missing data per individual
#Generates a file reporting the missingness on a per-individual basis = proportion of missing data per individual
missingness_indiv 


####libraries####
library(tidyverse)
library(readr)
library(ggplot2)
library(ggpubr)
####????####
##Mean depth
var_depth <- read_delim("good_rbrush/vcftools/good_sitedepth.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip =1)
a <- ggplot(var_depth, aes(mean_depth)) +
  geom_density(fill = "dodgerblue1",colour = "black", alpha = 0.3)+
  xlim(0, 200)
a + theme_light()
summary(var_depth$mean_depth)
a + theme_light() + xlim(0, 100)
a
##Variant missingness
var_miss <- read_delim("good_rbrush/vcftools/good_miss_site.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss","fmiss"), skip = 1)
b <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
b + theme_light()
summary(var_miss$fmiss)

##Minor alelle frecuencies
var_freq <- read_delim("good_rbrush/vcftools/good_alleles.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
c <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
c + theme_light()
summary(var_freq$maf)

##Mean depth per individual
ind_depth <- read_delim("good_rbrush/vcftools/good_depth.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
d <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
d + theme_light()

##Proportion of missing data per individual
ind_miss  <- read_delim("good_rbrush/vcftools/good_miss_indiv.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
e <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
e + theme_light()

####2good with R 0.80####
var_depth1 <- read_delim("good_rbrush/vcftools/good2_sitedepth.ldepth.mean", delim = "\t",
                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip =1)
a1 <- ggplot(var_depth1, aes(mean_depth)) +
  xlim(0, 100)+ geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a1 + theme_light()
summary(var_depth1$mean_depth)
var_miss1 <- read_delim("good_rbrush/vcftools/good2_miss_site.lmiss", delim = "\t",
                        col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss","fmiss"), skip = 1)
b1 <- ggplot(var_miss1, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
b1 + theme_light()
summary(var_miss1$fmiss)
var_freq1 <- read_delim("good_rbrush/vcftools/good2_alleles.frq", delim = "\t",
                        col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq1$maf <- var_freq1 %>% select(a1, a2) %>% apply(1, function(z) min(z))
c1 <- ggplot(var_freq1, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
c1 + theme_light()
summary(var_freq1$maf)
ind_depth1 <- read_delim("good_rbrush/vcftools/good2_depth.idepth", delim = "\t",
                         col_names = c("ind", "nsites", "depth"), skip = 1)
d1 <- ggplot(ind_depth1, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
d1 + theme_light()
ind_miss1  <- read_delim("good_rbrush/vcftools/good2_miss_indiv.imiss", delim = "\t",
                         col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
e1 <- ggplot(ind_miss1, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
e1 + theme_light()

####3good with R 0.80####
var_depth2 <- read_delim("good_rbrush/vcftools/good3_sitedepth.ldepth.mean", delim = "\t",
                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip =1)
a2 <- ggplot(var_depth2, aes(mean_depth)) +
  xlim(0, 100)+ geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a2 + theme_light()
summary(var_depth2$mean_depth)
var_miss2 <- read_delim("good_rbrush/vcftools/good3_miss_site.lmiss", delim = "\t",
                        col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss","fmiss"), skip = 1)
b2 <- ggplot(var_miss2, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
b2 + theme_light()
summary(var_miss1$fmiss)
var_freq2 <- read_delim("good_rbrush/vcftools/good3_alleles.frq", delim = "\t",
                        col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq2$maf <- var_freq2 %>% select(a1, a2) %>% apply(1, function(z) min(z))
c2 <- ggplot(var_freq1, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
c2 + theme_light()
summary(var_freq1$maf)
ind_depth2 <- read_delim("good_rbrush/vcftools/good3_depth.idepth", delim = "\t",
                         col_names = c("ind", "nsites", "depth"), skip = 1)
d2 <- ggplot(ind_depth2, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
d2 + theme_light()
ind_miss2  <- read_delim("good_rbrush/vcftools/good3_miss_indiv.imiss", delim = "\t",
                         col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
e2 <- ggplot(ind_miss2, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
e2 + theme_light()

####plot them all together####
mean_depth <- ggarrange(a, a1, a2,
                        labels = c("no .80 1", ".80 2", ".80 3"),
                        ncol = 3, nrow = 1)
variant_missingness <- ggarrange(b, b1, b2,
                                 labels = c("no .80", ".80 2", ".80 3"),
                                 ncol = 3, nrow = 1)
minor_allele <- ggarrange(c, c1, c2,
                          labels = c("no .80", ".80 2", ".80 3"),
                          ncol = 3, nrow = 1)
depth_indiv <- ggarrange(d, d1, d3,
                         labels = c("no .80", ".80 2", ".80 3"),
                         ncol = 3, nrow = 1)
missingness_indiv <- ggarrange(e, e1, e4,
                               labels = c("no .80", ".80 2", ".80 3"),
                               ncol = 3, nrow = 1)
