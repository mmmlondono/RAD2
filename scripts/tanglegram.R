####READ ME####
#Co-phylogeny plot: easy tanglegram in R
####load packages ####
#to install ggtree
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
#libraries
library(ggplot2)
library(ggtree)
library(phangorn)
library(dplyr)
library(ape)
library(phytools)

####first methods, tanglegram and paco ####
# Meta file
meta <- read.table('tanglegram/meta.csv', sep=',', header = T)
meta1 <- read.table('tanglegram/meta.csv', sep=',', header = T)

# Load tree 1
tree1 <- read.nexus('tanglegram/aciunex.tre')
tree1 <- midpoint(tree1)

# visualize labels of internal nodes:
plot(tree1, use.edge.length=FALSE)
nodelabels()
# rotate clades around nodes
tre1.new <- rotate(tree1, 33)
tre1.new <- rotate(tre1.new, 36)

# compare the results:
par(mfrow=c(1,2)) # split graphical device
plot(tree1) # plot old tre
plot(tre1.new) # plot new tree

# Load tree 2
tree2 <- read.nexus('tanglegram/rbrnex.tre')
tree2 <- midpoint(tree2)

# visualize labels of internal nodes:
plot(tree2, use.edge.length=FALSE)
nodelabels()
# rotate clades around node
tre2.new <- rotate(tree2, 30)
tre2.new <- rotate(tre2.new, 38)
# compare the results:
par(mfrow=c(1,2)) # split graphical device
plot(tree2) # plot old tre
plot(tre2.new) # plot new tree


#Combine the meta feature dataset with both phylogenetic trees and visualize how they look.
t1 <-ggtree(tre1.new)  %<+%  meta + geom_tiplab()
t2 <- ggtree(tre2.new)  %<+%  meta + geom_tiplab()

tre1.new <- rotate(t1, 33)
t1
t2
#Draw both trees in a single figure. We also want to flip tree 2, for which we need to change the x-coordinates in that tree.
d1 <- t1$data
d2 <- t2$data

d1$tree <-'t1'
d2$tree <-'t2'

d2$x <- max(d2$x) - d2$x + max(d1$x) +  max(d1$x)*0.3
pp <- t1 + geom_tree(data=d2)
pp 

# join d1 and d2 for dataset so that we can use the coordinates of the tips for making connections between both of the trees.
dd <- bind_rows(d1, d2) %>% 
  filter(isTip == TRUE)
dd1 <- as.data.frame(dd) 

#conditionally join the tips of both trees for the feature we are interested in. Connected tips will represent the same isolates.
green_tree <- dd1[which(dd1$label == 'Green'), c('label', 'x', 'y', 'tree')]
pp + geom_line(aes(x, y, group=label), data=green_tree, color='#009E73')

pp + geom_line(aes(x, y, group=label), data=dd1, color='#56B4E9')

#### using phytools ####
ass <- read.csv('tanglegram/ass.csv')

# Load tree 1
tree1 <- read.nexus('tanglegram/aciurina.tre')
tree1 <- midpoint(tree1)

# visualize labels of internal nodes:
plot(tree1, use.edge.length=FALSE)
nodelabels()
# rotate clades around nodes
#tre1.new <- rotate(tree1, 33)
#tre1.new <- rotate(tre1.new, 36)

# Load tree 2
tree2 <- read.nexus('tanglegram/rbrnex.tre')
tree2 <- midpoint(tree2)

# visualize labels of internal nodes:
plot(tree2, use.edge.length=FALSE)
nodelabels()

##cophylo function##
vcophylo<-cophylo(tree1,tree2, assoc=ass)
n_colors <- 36
random_colors <- RColorBrewer::brewer.pal(8, name = "Dark2")
cols <- setNames(colorRampPalette(random_colors)(n_colors), tree1$tip.label)

# Create the cophylo object with bat.tree, betaCoV.tree, and bat_virus.data
bat.cophylo <- cophylo(tree1, tree2, assoc = ass)

par(lend=3) 
plot(vcophylo,link.type="curved",
     fsize=c(0.7,0.6),link.lwd=2, link.lty="solid",pts=FALSE, link.col=make.transparent(cols[ass[,1]],0.5))
pies<-diag(1,Ntip(tree1)) 
colnames(pies)<-rownames(pies)<-names(cols) 
tiplabels.cophylo(pie=pies,piecol=cols[ vcophylo$trees[[1]]$tip.label], which="left",cex=0.2)


#### original phytools code ####
data(bat.tree)
data(betaCoV.tree)
data(bat_virus.data) 
head(bat_virus.data)

bat.cophylo<-cophylo(bat.tree,betaCoV.tree, assoc=bat_virus.data)
cols<-setNames(RColorBrewer::brewer.pal( n=7,name="Dark2"),bat.tree$tip.label)
par(lend=3) 
plot(bat.cophylo,link.type="curved",
                 fsize=c(0.7,0.6),link.lwd=2, link.lty="solid",pts=FALSE, link.col=make.transparent(cols[ bat_virus.data[,1]],0.5))
pies<-diag(1,Ntip(bat.tree)) 
colnames(pies)<-rownames(pies)<-names(cols) 
tiplabels.cophylo(pie=pies,piecol=cols[ bat.cophylo$trees[[1]]$tip.label], which="left",cex=0.2)
