####READ ME####
# The purpose of this script is to optimize the different de novo parameters for rabbit brush RAD data 
#Osuna-Mascaro:
#We allowed a maximum missing data of 20 % (--max-missing 0.8), 
#a minimum minor allele frequency of 0.03 (--maf 0.03) and 
#specified a thin value of 5 (--thin 5), which allows that no two sites are within the specified distance from one another. 
#Also, we only included sites with quality scores above 10 (--minQ 10).

####libraries####
library(tidyverse)

####m5_M3_n3####
##Mean depth
#Generates a file containing the mean depth per site averaged across all individuals.
#the number of reads that have mapped to this position (coverage)
#the more reads that cover a site, the higher confidence our basecall is. 
#Should we set a maximum and a minimum depth?
var_depth1 <- read_delim("rbrush/vcf/m7_M3_n3_sitedepth.ldepth.mean", delim = "\t",
                         col_names = c("chr", "pos", "mean_depth", "var_depth"), skip =1)
a <- ggplot(var_depth1, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_depth1$mean_depth)
a + theme_light() + xlim(0, 100)

##Variant missingness
#Generates a file reporting the missingness on a per-site basis.
#This is a measure of how many individuals lack a genotype at a call site.
#Typically missingness of 75-95% is used.
var_miss1 <- read_delim("rbrush/vcf/m7_M3_n3_miss_site.lmiss", delim = "\t",
                        col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss","fmiss"), skip = 1)
a <- ggplot(var_miss1, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0, 0.3)
summary(var_miss1$fmiss)

##Minor alelle frecuencies
#Outputs the allele frequency for each site in a file with the suffix ".frq".
#distribution of allele frequencies. This will help inform our minor-allele frequency (MAF) thresholds
#Usually then, it is best practice to produce one dataset with a good MAF threshold and keep another without any MAF filtering at all.
#Ideally you want an idea of the distribution of your allelic frequencies but 0.05 to 0.10 is a reasonable cut-off.
var_freq1 <- read_delim("rbrush/vcf/m7_M3_n3_alleles.frq", delim = "\t",
                        col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq1$maf <- var_freq1 %>% select(a1, a2) %>% apply(1, function(z) min(z))
a <- ggplot(var_freq1, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_freq1$maf)

##Mean depth per individual
#Generates a file containing the mean depth per individual.
#distribution of mean depth among individuals.
ind_depth1 <- read_delim("rbrush/vcf/m7_M3_n3_depth.idepth", delim = "\t",
                         col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth1, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
#most individuals were sequenced at around 35X - 45X.

##Proportion of missing data per individual
#Generates a file reporting the missingness on a per-individual basis = proportion of missing data per individual
ind_miss1  <- read_delim("rbrush/vcf/m7_M3_n3_miss_indiv.imiss", delim = "\t",
                         col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss1, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
#most individuals have 0 - 20% missing data but some have over 30% up to ~65%. have to look into those
