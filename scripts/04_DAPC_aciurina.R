####READ ME####
#this is the DAPC for quick population structure inference for aciurina
#from this site: https://connor-french.github.io/intro-pop-structure-r/ 

####load packages into environment####
library(adegenet) # for dapc
library(vcfR) # for reading in genetic data
library(tidyverse) # for manipulating and plotting data
library(LEA) # For sNMF
library(rnaturalearth) #for mapping

####read in####
#aciruina vcf file
focal_vcf <- read.vcfR("aciurina/all/aciurina2_.8.filtered.vcf")
aciurina_genlight <- vcfR2genlight(focal_vcf)
#make popmap
popmap<-data.frame(id=colnames(focal_vcf@gt)[2:length(colnames(focal_vcf@gt))],pop=substr(colnames(focal_vcf@gt)[2:length(colnames(focal_vcf@gt))], 1,8))

####DAPC####
#decide on number of PC's and clusters to retain
num_clust <- find.clusters(aciurina_genlight)
#run the DAPC
aciurina_dapc <- dapc(aciurina_genlight, num_clust$grp)
#plot the results
scatter(aciurina_dapc, posi.da="bottomright")

#tidy the data for plotting. 
dapc_data_df <-
  # as_tibble() converts the ind.coord matrix to a special data frame called a tibble. 
  as_tibble(aciurina_dapc$ind.coord, rownames = "individual") %>%
  # mutate changes or adds columns to a data frame. Here, we're adding the population and group assignment columns to the data frame
  mutate(population = popmap$pop,
         group = aciurina_dapc$grp)

dapc_data_df

#plot the data
#Can color the points according to your pre-defined populations and the dapc groups to see if it conforms to your hypothesis.
dapc_plot <-
  ggplot(dapc_data_df, aes(
    x = LD1,
    y = LD2,
    fill = population
  )) +
  geom_point(shape = 21, size = 3) +
  #reverse the color direction to better reflect Prates et al. 2018
  scale_fill_viridis_d(direction = -1) + 
  theme_bw(base_size = 16)

dapc_plot
