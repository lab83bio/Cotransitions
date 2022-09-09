#!/usr/bin/env Rscript 
#
# Read tables from the OrthoDB server
#
library(data.table)
library(tibble)
library(taxize)
library(treeio)
library(ape)
library(phangorn)
library(argparse)

# create parser object
parser <- ArgumentParser()

parser$add_argument("-l", "--level", type="character", default="Eukaryota",
                    help="taxonomic level, [default %(default)s]")
parser$add_argument("-s", "--species taxid", type="character",
                    help="only outgrups present in the species taxid, [default %(default)s]")
parser$add_argument("-m", "--min_percentage", type="double", default=1,
                    help="include orthogroups present in  at least x percent of taxa, [default %(default)s]")
parser$add_argument("-t", "--tree", type="character", default="raxml", choices=c('ncbi','raxml', 'random'),
                    help="order taxa by ncbi tree (partially resolved), ncbi-constrained raxml tree (fully resolved), or randomly  [default %(default)s]")
parser$add_argument("-o", "--out_root (none=no rooting)", type="character", default="Viridiplantae",
                    help="Group for tree rooting, [default %(default)s]")
parser$add_argument("-d", "--ladderize", type="character", default="RL", choices=c('LL','RL', 'NL'),
                    help="ladderize the tree for taxa order, [default %(default)s]")
args <- parser$parse_args()


oDB <- "https://v101.orthodb.org/download/odb10v1"
write(paste0("reading data from: ",oDB),stderr())

#NCBI taxonomy nodes (levels) where orthologous groups (OGs) are calculated
lev <- fread(paste0(oDB,'_levels.tab.gz'),header = F, 
            col.names=c('Taxid','Group','num.genes','num.orthogroups','num.orgs'))
#OrthoDB orthologous groups
OGs <- fread(paste0(oDB,'_OGs.tab.gz'),header = F,
            col.names=c('Orthogroup','Taxid','Description'))
#OGs to genes correspondence
Gcs <- fread(paste0(oDB,'_OG2genes.tab.gz'), header=F, showProgress =T,
            col.names=c('Orthogroup','GeneID'))

#Select level
sOGs <- OGs[Taxid %in% lev[Group %in% args$level,Taxid],Orthogroup]
Gs <- Gcs[Orthogroup %in% sOGs]
Gs$GenomeID=gsub(':.*','',Gs$GeneID)

if(!is.null(args$species)) {
  Gs$SpeciesID=gsub('_.*','',Gs$GenomeID)
  sOGs <- sOGs[sOGs %in% Gs[SpeciesID %in% args$species,Orthogroup]]
  Gs <- Gs[Orthogroup %in% sOGs]
}
  
#number of genes and genomes per Orthogroup
Ngenes <- Gs[,.(Total=.N),by="Orthogroup"]
Ngenomes <- Gs[,.N,by=c("Orthogroup","GenomeID")][,.N,by="Orthogroup"]
Gs_count <- merge(Ngenes,Ngenomes,by="Orthogroup")

#number of Orthogroups per genome
Ngroups<- Gs[,.N,by=c("Orthogroup","GenomeID")][,.N,by="GenomeID"]
write(paste0("Ortogroups: ",dim(Gs_count)[1]," Genomes: ",dim(Ngroups)[1], " Genes: ",dim(Gs)[1]),stderr())
Ngroups$Taxid <-  gsub( '_.*','',Ngroups$GenomeID )
Ngroups <- Ngroups[order(Taxid,-N)] #order by Num orthogroups
Dup <- Ngroups[duplicated(Ngroups$Taxid)]$GenomeID # duplicated species with less orthogroups
Ngroups <- Ngroups[!(GenomeID %in% Dup)]
write(paste0("Removed duplicated Taxid: ",Dup),stderr())


#############################################
#       Get taxonomic classification        #
#############################################
cl <- suppressMessages(taxize::classification(Ngroups$Taxid, db='ncbi'))
names(cl) <- Ngroups$GenomeID
classified <- names(cl[lengths(cl)>1])
unclassified <- names(cl[lengths(cl)<2])
write(paste0(length(classified),"/",length(Ngroups$Taxid)," organisms classified"),stderr())
if (length(unclassified)>0) {
  warning(paste0(length(unclassified)," unclassified taxid will be dropped: ",paste(unclassified,collapse=',')))
}
#################################################
# Calculate tree from taxonomic classification  #
#################################################
tr_cl <- class2tree(cl[classified])
tr <- tr_cl$phylo #tree

Ngroups$species <- tr$tip.label[match(Ngroups$GenomeID, tr_cl$names)]
tr$tip.label <- tr_cl$names #GenomeID as tip label
cl_d<-do.call("rbind", cl)
rank <- unlist(lapply(c(Ngroups$species,tr$node.label), function(x) cl_d[which(cl_d$name==x),"rank"][1]))

#################################################
#   order GenomeIDs according to tree order     #
#################################################
ordered_tips <- tr$edge[tr$edge[,2]<=length(tr$tip.label), 2]   #tips in tree order
Ngroups <- Ngroups[match(tr$tip.label[ordered_tips],GenomeID)] 
Gs <- Gs[Gs[.(Ngroups$GenomeID), on=.(GenomeID),which=T]] #order and drop duplicated species

#Discard orthogroups with <MIN_PERCENT  genomes
write(paste0("Ortogroups: ",dim(Gs_count)[1]," Genomes: ",dim(Ngroups)[1], " Genes: ",dim(Gs)[1]),stderr())
Keep <- Gs_count[Gs_count$N/nrow(Ngroups)*100 > args$min_percent,Orthogroup]
Gs<-Gs[Orthogroup %in% Keep]
write(paste0("Ortogroups in ", args$min_percent, "% of Genomes: ",length(Keep)," Genomes: ",dim(Ngroups)[1], " Genes: ",dim(Gs)[1]),stderr())

#Group GenomeID by Orthogroup in wide format (column order by tree)
Gst <- dcast(Gs, Orthogroup~factor(GenomeID,levels = Ngroups$GenomeID), value.var = "GeneID", function(x) list(x))

#Convert to numbers
Gst_num_m <- cbind(Gst[,1],Gst[,lapply(.SD,lengths),.SDcols=2:NCOL(Gst)]) #Number of genes
Gst_num <- cbind(Gst[,1],Gst_num_m[,lapply(.SD,function (x) ifelse(x>1,1,x)),.SDcols=2:NCOL(Gst)]) # 0,1 


base_name=paste(args$level, args$tree, args$ladderize, sep='.')

#################################################
#   taxonomy-constrained raxml tree             #
#################################################

phylip_tree=paste0(base_name,".tree")
write.tree(tr,file=phylip_tree) #tree in phylip fortmat
tG <- data.table::transpose(Gst_num[,2:ncol(Gst_num)],keep.names="S")
tG$S <- unlist (lapply(tG$S, function (x) paste(x,paste(rep(' ',10-nchar(x)),collapse=""))))
phylip_data=paste0(base_name,".data")
sink(phylip_data) # presence/absence data in phylip format
cat("\t",paste(c(dim(tG[,2:ncol(tG)]))),'\n') #nOrgs nOrthogroups
fwrite(tG, file ="", sep='#',quote=F, col.names=F)
sink()
system(paste("sed -i s/#//g",phylip_data)) #could not use empty sep

if (args$tree=='raxml'){
  
  RAxMLfile <- paste0("RAxML_bestTree.",base_name)
  if ( !file.exists(RAxMLfile) ){
    RAxML_opt <- "-m BINCAT -p 33 -T 8";
    RAXml <- paste("raxmlHPC-PTHREADS","-g",phylip_tree,"-s",phylip_data,"-n",base_name, RAxML_opt," > /dev/null")
  }
  
  tr <- read.tree(RAxMLfile)

  #assign internal nodes
  LCA <- function(taxids, cl){ #Lowest common ancestor of a taxid vector
    LCA_row <- length(Reduce(intersect, lapply(cl[taxids], "[[", 1)))
    LCA <- cl[taxids][[1]][LCA_row,]
    return (LCA)
  } 
  node_ids <- sort(unique(tr$edge[,1][tr$edge[,1]>length(tr$tip.label)])) #internal nodes
  tr$node.label <- unlist(lapply(phangorn::Descendants(tr, node_ids, "tips"), function(x) LCA(tr$tip.label[x],cl)$name))
  rank <- unlist(lapply(c(Ngroups$species,tr$node.label), function(x) cl_d[which(cl_d$name==x),"rank"][1]))
  
}

if(args$out_root != "none" && which(tr$node.label==args$out_root)>0){
  tr <- ape::root(tr,node=length(tr$tip.label)+which(tr$node.label==args$out_root)) #euk. root is uncertain see to https://doi.org/10.1016/j.tree.2019.08.008
}

if (args$ladderize=="RL") {
  tr <- ladderize(tr,right=T)
}
if (args$ladderize=="LL") {
  tr <- ladderize(tr,right=F)
}




#################################################
#   order GenomeIDs according to tree order     #
#################################################
ordered_tips <- tr$edge[tr$edge[,2]<=length(tr$tip.label), 2]   #tips in tree order
if (args$tree=='random'){
  ordered_tips <- ordered_tips[rev(sample(length(ordered_tips)))] #random order
}
Ngroups <- Ngroups[match(tr$tip.label[ordered_tips],Ngroups$GenomeID)] 
tr <- rotateConstr(tr,Ngroups$GenomeID) # order according the tree



setcolorder(Gst,c("Orthogroup",Ngroups$GenomeID))
setcolorder(Gst_num,c("Orthogroup",Ngroups$GenomeID))
setcolorder(Gst_num_m,c("Orthogroup",Ngroups$GenomeID))

#Add orthogroup description
Gst$Description <- OGs[match(Gst$Orthogroup,Orthogroup), Description]
setcolorder(Gst,c("Orthogroup","Description"))


#Save nexus_tree
all_nodes <- length(tr$tip.label) + tr$Nnode
taxa=unlist(c(Ngroups[match(tr$tip.label, Ngroups$GenomeID)]$species,tr$node.label))
d <- tibble(node=1:all_nodes,Taxa=taxa,Rank=rank) #assign also internal nodes
d$Taxa <- gsub("\\[|\\]", "_", d$Taxa) #square brackets not allowed in nexus names
tr_d <- treedata(phylo=tr,data=d) #tree with taxonomic information
treefile=paste0(base_name,".nexus")
write.beast(tr_d, tree.name="Species_tree", file=treefile) #Beast Nexus format (Figtree compatible)
write(paste0("Tree written in: ",treefile),stderr())


#write tab-separated table flattening the list
csv_file=paste0(base_name,".csv")
fwrite(Gst, file =csv_file, sep="\t", sep2 = c("",",",""),quote=F)
write(paste0("csv_file written in: ",csv_file),stderr())
csv_file_num_m=paste0(base_name,".csv.num_m")
fwrite(Gst_num_m, file =csv_file_num_m, sep="\t",quote=F, col.names=F)
write(paste0("dataset with gene numbers written in: ",csv_file_num_m),stderr())
csv_file_num=paste0(base_name,".csv.num")
fwrite(Gst_num, file =csv_file_num, sep="\t",quote=F, col.names=F)
write(paste0("dataset with presence/absence (1/0) written in: ",csv_file_num),stderr())

#all done





