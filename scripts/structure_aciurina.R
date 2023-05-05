####load packages into environment####
library(adegenet) # for dapc
library(vcfR) # for reading in genetic data
library(tidyverse) # for manipulating and plotting data
library(LEA) # For sNMF
library(rnaturalearth) #for mapping

#Read in VCF
aciurina_vcf <- read.vcfR("aciurina/all/aciurina2_.8.filtered.vcf")
#convert to genlight object
aciurina_genlight <- vcfR2genlight(aciurina_vcf)

# Input is sNMF object, the k value, and optionally an array of colors (has default of 10).
plot_sNMF <- function(input, k_val, colors = BlackYellow(k_val+1)){
  # picks best run best on cross entropy
  best_run <- which.min(cross.entropy(input, K = k_val))
  # makes q matrix of ancestry coeffs
  q_mat <- Q(input, K = 10, run = best)
  # plots the output, space makes blank between indivs
  barplot(t(q_matrix), col = colors, border = NA, space = 0.25, xlab = "Individuals", ylab = "Admixture coefficients", horiz=FALSE)
}

best_run <- which.min(cross.entropy(vanua_snmf, K = 10))
#select K
q_mat <- Q(vanua_snmf, K = 10, run = best_run)
colnames(q_mat) <- paste0("P", 1:3)
pops <- read.csv("pop_list_str.csv")
q_df <- q_mat %>%
  as_tibble() %>%
  # add the pops data for plotting
  mutate(individual = pops$ID,
         region = pops$pop,
         order = pops$plot_order)
q_df
q_df_long <- q_df %>%
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q")
q_df_prates <- q_df_long %>%
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>%
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))
q_palette <- c("#FDE725", "#35B779", "#440154")
q_df_prates %>%
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette) +
  #scale_fill_viridis_d() +
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


