#!/usr/bin/env Rscript 
#
# Read tables from the OrthoDB server
# save .RData object
#
library(data.table)
library(taxize)
library(ape) # write.tree
library(argparse)

# create parser object
ap <- ArgumentParser()
ap$add_argument("-l", "--level", type="character", default="Eukaryota",
                    help="taxonomic level, [default %(default)s]")
ap$add_argument("-s", "--species taxid", type="character",
                    help="only outgrups present in the species taxid, [default %(default)s]")
ap$add_argument("-m", "--min_percentage", type="double", default=1,
                    help="include orthogroups present in  at least x percent of taxa, [default %(default)s]")
args <- ap$parse_args()

print(args,stderr())

oDB <- "https://v101.orthodb.org/download/odb10v1"
write(paste0("reading data from: ",oDB),stderr())

#NCBI taxonomy nodes (levels) where orthologous groups (OGs) are calculated
lev <- fread(paste0(oDB,'_levels.tab.gz'),header = F, 
            col.names=c('Taxid','Group','num.genes','num.orthogroups','num.orgs'), tmpdir='.')
#OrthoDB orthologous groups
OGs <- fread(paste0(oDB,'_OGs.tab.gz'),header = F,
            col.names=c('Orthogroup','Taxid','Description'), tmpdir='.')
#OGs to genes correspondence
Gcs <- fread(paste0(oDB,'_OG2genes.tab.gz'), header=F, showProgress =T,
            col.names=c('Orthogroup','GeneID'), tmpdir='.')

#Select level
sOGs <- OGs[Taxid %in% lev[Group %in% args$level,Taxid],Orthogroup]
Gs <- Gcs[Orthogroup %in% sOGs]
Gs$GenomeID=gsub(':.*','',Gs$GeneID)

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
#   order GenomeIDs according to ncbi tree order     #
#################################################
ordered_tips <- tr$edge[tr$edge[,2]<=length(tr$tip.label), 2]   #tips in tree order
Ngroups <- Ngroups[match(tr$tip.label[ordered_tips],GenomeID)] 
Gs <- Gs[Gs[.(Ngroups$GenomeID), on=.(GenomeID),which=T]] #order and drop duplicated species

#Discard orthogroups with <MIN_PERCENT genomes
write(paste0("Ortogroups: ",dim(Gs_count)[1]," Genomes: ",dim(Ngroups)[1], " Genes: ",dim(Gs)[1]),stderr())
Keep <- Gs_count[Gs_count$N/nrow(Ngroups)*100 > args$min_percent,Orthogroup]
Gs<-Gs[Orthogroup %in% Keep]
write(paste0("Ortogroups in ", args$min_percent, "% of Genomes: ",length(Keep)," Genomes: ",dim(Ngroups)[1], " Genes: ",dim(Gs)[1]),stderr())

#Group GenomeID by Orthogroup in wide format (column order by ncbi tree)
Gst <- dcast(Gs, Orthogroup~factor(GenomeID,levels = Ngroups$GenomeID), value.var = "GeneID", list) 

#Add orthogroup description
Gst$Description <- OGs[match(Gst$Orthogroup,Orthogroup), Description]
setcolorder(Gst,c("Orthogroup","Description"))

#bincat -  binary absence/presence data for tree resolution
bincat <- dcast(Gs[,1,by=c("GenomeID","Orthogroup")],GenomeID ~ Orthogroup, value.var="V1",fill=0)


######################################################
#   save tree and 0/1 data in phylip format          #
######################################################
base_name=paste(args$level, "phylip", sep='.')

phylip_tree=paste0(base_name,".tree")
write.tree(tr,file=phylip_tree) #tree in phylip fortmat
bincat$GenomeID <- unlist (lapply(bincat$GenomeID, function (x) paste(x,paste(rep(' ',10-nchar(x)),collapse="")))) #names in phylip format
phylip_data=paste0(base_name,".data")
sink(phylip_data) # presence/absence data in phylip format
cat("\t",paste(c(dim(bincat[,2:ncol(bincat)]))),'\n') #nOrgs nOrthogroups
fwrite(bincat, file ="", sep='#',quote=F, col.names=F)
sink()
system(paste("sed -i s/#//g",phylip_data)) #could not use empty sep

#Delete missing orthogroups if -s
if(length(args$species)>0){
  scol <- grep(paste0(args$species,"_"),colnames(Gst))
  Gst <- Gst[which(lengths(Gst[,..scol][[1]])>0),]
  write(paste0("Ortogroups in ",taxid2name(args$species),": ",dim(Gst)[1]), stderr())
}

outfile = paste0(args$level,".RData")

save(Gst,
     tr,
     Ngroups,
     cl_d,
     file = outfile,
     compress =T)

write(paste0("Data written in: ",outfile),stderr())

#all done





