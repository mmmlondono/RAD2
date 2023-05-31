####READ ME####
#this script is to install DartR (since it is so problematic)

####libraries####
library(devtools)
install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SNPRelate")

install.packages("terra")
install.packages("rlang")
library(rlang)

#if trouble with SNPRelate run this 
install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "qvalue")) 

install.packages("dartR")
library("dartR")

library(vcfR)

####read in vcf####
vcf <- read.vcfR(file="/Users/martinson_lab/Documents/RAD2/ipyrad/vcf/iptest.vcf")
vcf_gl <- gl.read.vcf("/Users/martinson_lab/Documents/RAD2/ipyrad/vcf/iptest.vcf")

gl <- vcfR2genlight(vcf)

test <- gl.compliance.check(gl)

glmono <- gl.filter.monomorphs(gl, verbose = NULL)
