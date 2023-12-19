# the source is the path to the plotting_funcs.R file
# (see TreemixManual for more info: 
## http://gensoft.pasteur.fr/docs/treemix/1.12/treemix_manual_10_1_2012.pdf)
# the input for "plot_tree" is the path to the folder your treemix output
## with the stem of the output file names (e.g. sympos for me)
# includes the 
#source("D:/Documents/Projects/PachyHybrid/treemix/plotting_funcs.R")
#setwd("D:/Documents/Projects/PachyHybrid/treemix")

source("wolbachia/treemix/treemix-1.13/src/plotting_funcs.R")

####Unsupervised####
par(mfrow=c(2,3))

plot_tree("wolbachia/treemix/out/wol_jk0")
plot_tree("wolbachia/treemix/out/wol_jk1")
plot_tree("wolbachia/treemix/out/wol_jk2")
plot_tree("wolbachia/treemix/out/wol_jk3")
plot_tree("wolbachia/treemix/out/wol_jk4")
plot_tree("wolbachia/treemix/out/wol_jk5")

plot_resid("wolbachia/treemix/out/wol_jk3", "wolbachia/treemix/out/poporder1.csv")


####SNM####
par(mfrow=c(1,3))

plot_tree("wolbachia/treemix/out_snm/wol_jk0")
plot_tree("wolbachia/treemix/out_snm/wol_jk1")
plot_tree("wolbachia/treemix/out_snm/wol_jk2")

####TES####
par(mfrow=c(1,3))

plot_tree("wolbachia/treemix/out_tes/wol_jk0")
plot_tree("wolbachia/treemix/out_tes/wol_jk1")
plot_tree("wolbachia/treemix/out_tes/wol_jk2")

####SLV####
par(mfrow=c(1,3))

plot_tree("wolbachia/treemix/out_slv/wol_jk0")
plot_tree("wolbachia/treemix/out_slv/wol_jk1")
plot_tree("wolbachia/treemix/out_slv/wol_jk2")
plot_tree("wolbachia/treemix/out_slv/wol_jk3")
plot_tree("wolbachia/treemix/out_slv/wol_jk4")

par(mfrow=c(1,3))

plot_tree("wolbachia/treemix/out_slv2/wol_jk0")
plot_tree("wolbachia/treemix/out_slv2/wol_jk1")
plot_tree("wolbachia/treemix/out_slv2/wol_jk2")
plot_tree("wolbachia/treemix/out_slv2/wol_jk3")
plot_tree("wolbachia/treemix/out_slv2/wol_jk4")
