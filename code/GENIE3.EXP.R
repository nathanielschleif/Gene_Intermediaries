library(netbenchmark)
library(GENIE3)
library(doParallel)
library(doRNG)
source("https://bioconductor.org/biocLite.R")
biocLite("GENIE3")


setwd('/home/wcic/Documents/GENIE3-Maize')
data <- read.csv(file='../data/walley.data.csv',header=TRUE,sep=',',row.names=1)
data <- read.csv(file='../data/arabidopsis.meristem.expression.csv',header=TRUE,sep=',',row.names=1)
data_no0 <- data[rowSums(data)!=0,]
head(as.matrix(data))
weights <- GENIE3(exprMatrix=as.matrix(data_no0),nTrees=10000,verbose=TRUE)
weight_test = GENIE3(exprMatr)

# Test - Expected form
exprMatr <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
rownames(exprMatr) <- paste("Gene", 1:20, sep="")
colnames(exprMatr) <- paste("Sample", 1:5, sep="")
head(exprMatr)
max(weights)
hist(weights)
flat_w = as.vector(weights)
summary(flat_w)
sum(flat_w>0.007)
top_10k = flat_w[rank(-flat_w)<=10000]
top_100k = flat_w[rank(-flat_w)<=100000]

# Null Weights:
null_data = t(apply(data_no0, 1, sample))
colnames(null_data) = colnames(data_no0)
dim(null_data)
null_weights = GENIE3(null_data,verbose=TRUE)

hist(null_weights)
max(null_weights)
flat_null_w = as.vector(null_weights)
summary(flat_null_w)
sum(flat_null_w>0.007)
null_top_10k = flat_null_w[rank(-flat_null_w)<=10000]
null_top_100k = flat_null_w[rank(-flat_null_w)<=100000]

# Complete Random Weights:
comp_null_data = matrix(sample(as.vector(as.matrix(data_no0))),nrow=nrow(data_no0),ncol=ncol(data_no0))
colnames(comp_null_data) = colnames(data_no0)
rownames(comp_null_data) = rownames(data_no0)
dim(comp_null_data)
comp_null_weights = GENIE3(comp_null_data,verbose=TRUE)

hist(comp_null_weights)
max(comp_null_weights)
flat_comp_null = as.vector(comp_null_weights)
comp_null_10k = flat_comp_null[rank(-flat_comp_null)<=10000]

data_rand = data.frame(data_no0)
data_rand = slice(data_rand,sample(1:n()))
data_rand = apply(as.matrix(data_rand),1,sample)
data_rand = sapply(data_rand,sample)
dim(data_rand)
data_rand[25000]

# Diversity Dataset

test_set = read.table('../data/test.txt',header=TRUE)
#temp = list.files(pattern="*")
#myfiles = lapply(temp, read.delim)
#test2 = myfiles[[1]]
row.names(test_set) = test_set$Gene
names(test_set) = c('bur1','bur2','can1','can2','col1','col2','ct1','ct2','edi1','edi2','hi1','hi2','kn1','kn2','ler1','ler2','mt1','mt2','no1','no2','oy1','oy2','po1','po2','rsch1','rsch2','sf1','sf2','tsu1','tsu2','wil1','wil2','ws1','ws2','wu1','wu2','zu1','zu2')
test_set = test_set[,even_indexes]
div_data_no0 <- test_set[rowSums(test_set)!=0,]
div_data_no0 <- div_data_no0[rowSums(div_data_no0)>=0,]
even_indexes<-seq(2,76,2)

test2 = div_data_no0[grepl("*TE*",test_set$Gene),]

## CLEAN IMPORT OF DATA ##
even_indexes<-seq(2,76,2)
div_data = read.table('../data/comb.GEO.Diversity.txt',header=TRUE)
row.names(div_data) = div_data$Gene
div_data = div_data[!grepl(".*TE.*",div_data$Gene),]
div_data = div_data[!grepl(".*new.*",div_data$Gene),]
div_data = div_data[,even_indexes]
div_data = div_data[rowSums(div_data)>0,]
names(div_data) = c('bur1','bur2','can1','can2','col1','col2','ct1','ct2','edi1','edi2','hi1','hi2','kn1','kn2','ler1','ler2','mt1','mt2','no1','no2','oy1','oy2','po1','po2','rsch1','rsch2','sf1','sf2','tsu1','tsu2','wil1','wil2','ws1','ws2','wu1','wu2','zu1','zu2')
mat_div_data = as.matrix(div_data)

phos_regs = as.character(read.table(file = '../data/genelist.protein.phos.txt')$V1)

div_weights = GENIE3(mat_div_data,regulators=phos_regs,verbose=TRUE)
hist(div_weights)
max(div_weights)
flat_div = as.vector(div_weights)
div_10k = flat_div[rank(-flat_div)<=10000]
div_100k = flat_div[rank(-flat_div)<=100000]
min(div_10k)
length(div_10k[div_10k>0.009])

# Now need to randomize the data
div_null_data = t(apply(mat_div_data, 1, sample))
colnames(div_null_data) = colnames(mat_div_data)
dim(div_null_data)
null_weights = GENIE3(div_null_data,verbose=TRUE)
hist(null_weights)
max(null_weights)
flat_null = as.vector(null_weights)
comp_null_10k = flat_null[rank(-flat_null)<=10000]
min(comp_null_10k)
length(comp_null_10k[comp_null_10k>0.009])
null_100k = flat_null[rank(-flat_null)<=100000]


names(test_set) = c('bur1','bur2','can1','can2','col1','col2','ct1','ct2','edi1','edi2','hi1','hi2','kn1','kn2','ler1','ler2','mt1','mt2','no1','no2','oy1','oy2','po1','po2','rsch1','rsch2','sf1','sf2','tsu1','tsu2','wil1','wil2','ws1','ws2','wu1','wu2','zu1','zu2')
test_set = test_set[,even_indexes]
div_data_no0 <- test_set[rowSums(test_set)!=0,]
div_data_no0 <- div_data_no0[rowSums(div_data_no0)>=0,]
test2 = div_data_no0[grepl("*TE*",test_set$Gene),]

div_weights = GENIE3(as.matrix(div_data_no0),verbose=TRUE)


