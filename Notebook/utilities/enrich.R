#!/usr/bin/env Rscript

library(clusterProfiler)
library(dplyr)
library(data.table)
library("argparse")
library(rlang)

parser <- ArgumentParser()
parser$add_argument("-i", "--input_file", help="input file")
parser$add_argument("-cut", "--pval_cutoff", help="P-value threshold")
args <- parser$parse_args()

file = args$input_file
cutoff = args$pval_cutoff

BP = read.csv(file, sep='\t')
BP = BP[order(-BP$score),]

BP_nodup = distinct(BP, name, .keep_all = TRUE)
geneList = setNames(as.numeric(BP_nodup$score), as.character(BP_nodup$name))

geneData = as.data.frame(BP[,c('feature','name')])
geneDataNames = as.data.frame(BP[,c('feature','description')])

BPgsea = GSEA(geneList, TERM2GENE=geneData, TERM2NAME=geneDataNames, 
              pvalueCutoff=1, seed=2022, maxGSSize=100000, minGSSize=1, eps = 0)

print(paste(basename(file), '.enriched', sep=''))
write.csv(as.data.frame(BPgsea), paste(basename(file), '.enriched', sep=''))
