# the source is the path to the plotting_funcs.R file
# (see TreemixManual for more info: 
## http://gensoft.pasteur.fr/docs/treemix/1.12/treemix_manual_10_1_2012.pdf)
# the input for "plot_tree" is the path to the folder your treemix output
## with the stem of the output file names (e.g. sympos for me)
# includes the 
#source("D:/Documents/Projects/PachyHybrid/treemix/plotting_funcs.R")
#setwd("D:/Documents/Projects/PachyHybrid/treemix")

source("wolbachia/treemix/treemix-1.13/src/plotting_funcs.R")
plot_tree("outstem")

par(mfrow=c(2,2), mar=c(0,2,0,0), oma=c(1,2,1,0))

plot_tree("wolbachia/treemix/out/wol_jk0")
plot_tree("wolbachia/treemix/out/wol_jk1")
plot_tree("wolbachia/treemix/out/wol_jk2")
plot_tree("wolbachia/treemix/out/wol_jk3")
plot_tree("wolbachia/treemix/out/wol_jk4")

source("wolbachia/treemix/treemix-1.13/src/plotting_funcs.R")
plot_resid("wolbachia/treemix/out/wol_jk3", "wolbachia/treemix/out/poporder")



