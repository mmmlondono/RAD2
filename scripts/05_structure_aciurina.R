####READ ME####
# The purpose of this script is to make a structure-like plot called sNMF.
#sNMF is a fast and efficient method for estimating individual ancestry coefficients based on sparse non-negative matrix factorization algorithms. 

####load packages into environment####
library(adegenet) # for dapc
library(vcfR) # for reading in genetic data
library(tidyverse) # for manipulating and plotting data
library(LEA) # For sNMF
library(rnaturalearth) #for mapping
library(RColorBrewer)


#Read in VCF
aciurina_vcf <- read.vcfR("aciurina/all/aciurina2_.8.filtered.vcf")
#convert to genlight object
aciurina_genlight <- vcfR2genlight(aciurina_vcf)

####create snmf file####
#takes for ever#
vcf2geno("aciurina/aciurina2_.8.filtered.vcf", output = "aciurina/all/snmf/aciurina_all.geno")
aciurina_snmf = snmf("aciurina_all.geno", ploidy=2, 
                     K = 1:10, alpha = 10, project = "new", entropy = T, repetitions = 50)
#so save and load#
saveRDS(aciurina_snmf, file = "aciurina_snmf")

best_run <- which.min(cross.entropy(aciurina_snmf, K = 10))
#select K
q_mat <- Q(aciurina_snmf, K = 10, run = best_run)
colnames(q_mat) <- paste0("P", 1:10)
pops <- read.csv("popmap_order.csv")
q_df <- q_mat %>%
  as_tibble() %>%
  # add the pops data for plotting
  mutate(individual = pops$id,
         region = pops$pop,
         order = pops$order)
q_df
q_df_long <- q_df %>%
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q")
q_df_prates <- q_df_long %>%
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>%
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))
q_palette <- c("darkblue", "turquoise", "skyblue","seagreen","khaki", "lightgreen", "brown",
                         "hotpink","coral", "gold")
                         q_palette <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2")
                         q_palette <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd")   
                         
                         q_df_prates %>%
                           ggplot() +
                           geom_col(aes(x = individual, y = q, fill = pop)) +
                           scale_fill_manual(values = q_palette) +
                           #scale_fill_viridis_d() 
                           labs(fill = "Region") +
                           theme_minimal() +
                           # some formatting details to make it pretty
                           theme(panel.spacing.x = unit(0, "lines"),
                                 axis.line = element_blank(),
                                 axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                 strip.background = element_rect(fill = "transparent", color = "black"),
                                 panel.background = element_blank(),
                                 axis.title = element_blank(),
                                 panel.grid = element_blank()
                           )