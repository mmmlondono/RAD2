####High-Dimensional Data Analysis####

####packages####
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)

####Distance exercises #1####
data(tissuesGeneExpression)
head(e)
head(tissue)
#1. How many biological replicates are there for hippocampus?
length(grep("hippocampus", tissue)) 
#2. What is the distance between samples 3 and 45?
x <- e[,3]
y <- e[,45]
sqrt(sum((x-y)^2))
#or
sqrt(crossprod(x-y))
#3. What is the distance between gene 210486_at and 200805_at?
grep("210486_at", colnames(tissue))
x <- e["210486_at",]
y <- e["200805_at",]
sqrt(crossprod(x-y))
#4. If I run the command (don't run it!):
#d = as.matrix(dist(e))
#How many cells (number of rows times number of columns) would this matrix have?
#Compute the distance between all pairs of samples:
d = dist(t(e))
#How many distances are stored in d? (Hint: What is the length of d)?
#Why is the answer above not
ncol(e)^2

####Projections Exercises####
#packages
library(Biobase)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
library(tissuesGeneExpression)

#data
data(GSE5859Subset)
y = geneExpression[,1:2]
data(tissuesGeneExpression)

s = svd(e)
signflips = sample(c(-1,1),ncol(e),replace=TRUE)
signflips

newu= sweep(s$u,2,signflips,FUN="*")
newv= sweep(s$v,2,signflips,FUN="*" )
all.equal( s$u %*% diag(s$d) %*% t(s$v), newu %*% diag(s$d) %*% t(newv))

s = svd(e)
m = rowMeans(e)

cor(s[["u"]], m) 

newmeans = rnorm(nrow(e)) ##random values we will add to create new means
newe = e+newmeans ##we change the means
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(newe[,3]-newe[,45]))

y = e - rowMeans(e)
s = svd(y)

resid = y - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))

x=matrix(rep(c(1,2),each=5),5,2)
x
x*c(1:5)

sweep(x,1,1:5,"*")

diag(s$d)%*%t(s$v)

z = s$d * t(s$v)
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(y[,3]-y[,45]))
sqrt(crossprod(z[,3]-z[,45]))

ed <- sqrt(crossprod(e[,3]-e[,45]))
zd <- sqrt(crossprod(e[,3]-e[,45]))

distances = sqrt(apply(e[,-3]-e[,3],2,crossprod))
