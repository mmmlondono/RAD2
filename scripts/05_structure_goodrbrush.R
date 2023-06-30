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
#devtools::install_github("katiejolly/nationalparkcolors") #ntl parks colors
library(nationalparkcolors)

#Read in VCF
goodrbrush_vcf <- read.vcfR("good_rbrush/rbrush_good.filtered.vcf")
popmap<-data.frame(id=colnames(goodrbrush_vcf@gt)[2:length(colnames(goodrbrush_vcf@gt))],pop=substr(colnames(goodrbrush_vcf@gt)[2:length(colnames(goodrbrush_vcf@gt))], 7,8))
#convert to genlight object
goodrbrush_genlight <- vcfR2genlight(goodrbrush_vcf)

####create snmf file####
#takes for ever#
vcf2geno("good_rbrush/rbrush_good.filtered.vcf", output = "good_rbrush/rbrush_good.geno")
goodrbrush_snmf = snmf("good_rbrush/rbrush_good.geno", ploidy=2, 
                     K = 1:10, alpha = 10, project = "new", entropy = T, repetitions = 50)
#so save and load#
saveRDS(goodrbrush_snmf, file = "goodrbrush_snmf")

best_run <- which.min(cross.entropy(goodrbrush_snmf, K = 8))
#select K 7####
q_mat <- Q(goodrbrush_snmf, K = 8, run = best_run)
colnames(q_mat) <- paste0("P", 1:8)
pops <- read.csv("good_rbrush/popmap_ord.csv")
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

q_palette <- hcl.colors(9, palette = "spectral")
q_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                        "#F0E442", "#0072B2", "#CC79A7")
#q_palette  <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#"gold", "#bebada", "#8dd3c7", "#bc80bd", "darkblue", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#8dd3c7","#fb8072"
                        
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

#------------------------------------------------------------------------------#
####Ethan's way####
                        
q_palette <- c("#8dd3c7", "gold", "#bebada", "#fb8072", "#80b1d3",
                        "#fdb462", "#b3de69", "#fccde5", "#bc80bd", "darkblue") 
                                                
# Plotting function used to plot sNMF output.
# Input is sNMF object, the k value, and optionally an array of colors (has default of 10).
plot_sNMF <- function(input, k_val, colors = q_palette(k_val+1)){
  # picks best run best on cross entropy
  best_run2 <- which.min(cross.entropy(input, K = k_val))
  # makes q matrix of ancestry coeffs
  q_matrix <- Q(input, K = k_val, run = best_run2)
  # plots the output, space makes blank between indivs
  barplot(t(q_matrix), col = colors, border = NA, space = 0.25, xlab = "Individuals", ylab = "Admixture coefficients", horiz=FALSE)
  }
                        
# plot to decide optimal K, although don't put too much stock in k values
plot(aciurina_snmf, cex = 1.2, col = "lightblue", pch = 19)
                        
# here I override the default colors
par(mfrow=c(6,1), mar=c(0,2,0,0), oma=c(1,2,1,0))
plot_sNMF(aciurina_snmf, 2, c("gold", "#bebada"))
plot_sNMF(aciurina_snmf, 3, c("gold", "#bebada", "#b3de69"))
plot_sNMF(aciurina_snmf, 4, c("gold", "#bebada", "#b3de69", "#80b1d3"))
plot_sNMF(aciurina_snmf, 5, c("gold", "#bebada", "#b3de69", "#80b1d3", "#fb8072"))
plot_sNMF(aciurina_snmf, 6, c("gold", "#bebada", "#b3de69", "#80b1d3", "#fb8072","#fdb462"))
plot_sNMF(aciurina_snmf, 7, c("gold", "#bebada", "#b3de69", "#80b1d3", "#fb8072","#fdb462", "#fccde5"))                         
                        