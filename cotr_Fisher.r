#!/usr/bin/env Rscript 
#
# Fisher probability from a gene co-transition table 
#
library(data.table)# Fast fread
library(corpora)   # Fast fisher.pval
library(argparse) 

ap <- ArgumentParser()
ap$add_argument("-p", "--p.cutoff", type="double", default=1e-3, help="p cutoff")
ap$add_argument("-pa", "--padj.cutoff", type="double", default=1e-3, help="p.adj cutoff")
ap$add_argument("csv", nargs=1, help="*.transitions tsv file")
 
args <- ap$parse_args()

csv<-args$csv;
if (csv == '-') csv='file:///dev/stdin' 

p.cutoff=args$p.cutoff
padj.cutoff=args$padj.cutoff

##===================================
##    Read *.transitions file       =
##===================================

column_classes <- c(rep("factor",2),rep("integer",6))
D <- fread(csv, header=T, colClasses = column_classes, showProgress=T)

##===================================
##     Calculate scores             =
##===================================

D$k_score <- D$k/(D$t1+D$t2-abs(D$k)) # Morandin formula

##===================================
##   Calculate p.values and p.adj   =
##===================================
write(paste0("calculating probabilities (output table to STDOUT)"), stderr())
norg <- D$orgs
k <- abs(D$k)
M<-as.matrix(cbind(k, D$t2, D$t1-k, norg-D$t2)) #k1 n1 k2 n2 for fisher.pval
p=c() # vector of p.values 
for(i in 1:nrow(M)) 
  p[i] <- tryCatch(fisher.pval(M[i,1],M[i,2],M[i,3],M[i,4], alternative="greater"), 
                   error=function(error_message){return(1)})

p.adj <- p.adjust(p, method="holm")

##=========================================
##  Write results if p and pa < p.cutoff  =
##=========================================

D <- cbind(D,p,p.adj)
D <- D[p.adj <= padj.cutoff,]
D <- D[p <= p.cutoff,]
D <- D[order(D$p),]

write.table(format(D,digits=5), file = "" , quote=F, sep="\t", row.names=F)

#all done:
write(paste0(nrow(D) ," gene pairs with p <= ",p.cutoff, " and p.adj <= ",padj.cutoff), stderr())




