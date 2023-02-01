#!/usr/bin/env Rscript 
#
# Order species according to a phylogenetic tree
#
library(data.table)


                    help="ladderize the tree for taxa order, [default %(default)s]")
ap$add_argument("RData", nargs=1,help=".Rdata file from procedure_Orthodb_read_tables.r")
args <- ap$parse_args()

print(args,stderr())

load(args$RData)

base_name=gsub('.RData','',args$RData)

trp <- tr #phylip tree

#################################################
#   taxonomy-constrained raxml tree             #
#################################################

if (args$tree=='raxml'){
  RAxMLfile <- paste0("RAxML_bestTree.",base_name)
  if (file.exists(RAxMLfile) ){
    tr <- read.tree(RAxMLfile)
  }
  else {
	stop(paste("File ",RAxMLfile, "not found"))
  }
}

#assign names to internal nodes based on phylip tree
tr_md5sum <- makeNodeLabel(tr, method = "md5sum")
trp_md5sum <- makeNodeLabel(trp, method = "md5sum")
tr$node.label <- trp$node.label[match(tr_md5sum$node.label,trp_md5sum$node.label)]
rank <- unlist(lapply(c(Ngroups$species,tr$node.label), function(x) cl_d[which(cl_d$name==x),"rank"][1]))
  


if(args$out_root != "none" && args$out_root %in% tr$node.label){
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

setcolorder(Gst,c("Orthogroup", "Description",Ngroups$GenomeID))

#Convert to numbers
Gst_num_m <- cbind(Gst[,1],Gst[,lapply(.SD,lengths),.SDcols=3:NCOL(Gst)]) #Number of genes
Gst_num <- cbind(Gst[,1],Gst_num_m[,lapply(.SD,function (x) ifelse(x>1,1,x)),.SDcols=2:NCOL(Gst_num_m)]) # binary (0,1) 

base_name=paste(paste(base_name, args$tree, args$ladderize, sep='.'))

#Save nexus_tree
all_nodes <- length(tr$tip.label) + tr$Nnode
taxa=unlist(c(Ngroups[match(tr$tip.label, Ngroups$GenomeID)]$species,tr$node.label))
d <- tibble(node=1:all_nodes,Taxa=taxa,Rank=rank) #assign also internal nodes
d$Taxa <- gsub("[^[:alnum:][:blank:]-]", "_", d$Taxa) #some special chars (e.g. brackets) not allowed in nexus names
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
