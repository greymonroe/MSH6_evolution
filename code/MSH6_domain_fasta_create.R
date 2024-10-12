library(Biostrings)
library(seqinr)
library(ape)
library(phylogram)
library(dendextend)
library(msa)

sequences <- readAAStringSet(filepath = "~/Documents/all.faa")

interpro<-read_interproscan("~/Documents/all.faa.tsv")

#orthos<-fread("~/Documents/potential_orthologs_mismatch.txt")
orthos<-fread("~/Documents/predicted_msh6_blast.txt")
orthos<-rbind(data.table(Organism=predicted_msh6_blast_annotated$Organism, Protein=predicted_msh6_blast_annotated$protein),
              data.table(Organism=orthos$Organism, Protein=orthos$Protein))

interpro$Organism<-orthos$Organism[match(interpro$Protein, orthos$Protein)]
interpro_sub<-interpro[Protein %in% orthos$Protein & Analysis!="SUPERFAMILY"]
interpro_sub<-interpro[Analysis %in% c("CDD", "Pfam","SMART", "ProSiteProfiles")
                       & Protein %in% interpro_sub$Protein
                       & (grepl("MutS", Signature_Description)
                         )]

interpro_sub$Name<-paste(interpro_sub$Organism, interpro_sub$Protein, sep=" | ")
best_model<-interpro_sub[grepl("MutS", Signature_Description) & grepl("domain", Signature_Description),
                         .(I=sum(Signature_Description=="MutS domain I"),
                            II=sum(Signature_Description=="MutS domain II"),
                            III=sum(Signature_Description=="MutS domain III"),
                            IV=sum(Signature_Description=="MutS family domain IV"),
                            V=sum(Signature_Description=="MutS domain V"),
                            Length=max(Stop_Location),
                            start=min(Start_Location),
                            stop=max(Stop_Location)), by=.(Organism, Protein)]

# Count the number of zeros in columns I, II, III, IV, and V
best_model[, zero_count := rowSums(.SD == 0), .SDcols = c("I", "II", "III", "IV", "V")]
best_model[, greater_than_one_count := rowSums(.SD > 1), .SDcols = c("I", "II", "III", "IV", "V")]

# Order by zero_count and Length
setorder(best_model, Organism, zero_count, Length)

# Subset to get the top row for each Organism
#result <- best_model[zero_count==0 &greater_than_one_count==0 , .SD[1], by = Organism]
result <- best_model[zero_count==0 &greater_than_one_count==0 ,]
result<-result[!duplicated(result[,-2])]

tree<-read.tree("files/phyloT_generated_tree_Aug12_23.txt")
tree$tip.label<-gsub("_"," ",tree$tip.label)
plants<-organisms_in_node(tree, msh6_good, c("Viridiplantae"))
aglae<-organisms_in_node(tree, msh6_good, c("Chlorophyta"))
meta<-organisms_in_node(tree, msh6_good, "Metazoa")
fungi<-organisms_in_node(tree, msh6_good_genera, "Fungi")
other<-organisms_in_node(tree, msh6_good, c("Fungi","Metazoa","Viridiplantae"), exclude = T)
small_tree<-rbindlist(list(fungi, plants, meta, other), fill=T)

subresult<-result[Organism %in% c(plants$Organism, other$Organism)]

MutS_domains<-interpro_sub[Organism %in% subresult$Organism & grepl("MutS", Signature_Description) & grepl("domain", Signature_Description) &!grepl("smr",Signature_Description, ignore.case = T )]
table(MutS_domains$Signature_Description)

MutS_domains$lab<-paste0(MutS_domains$Protein,"-", gsub(" ", "_", MutS_domains$Organism))


msh6_domain_sequences<-apply(subresult, 1, function(r){
  protein=r["Protein"]
  protein_MutS_domains<-MutS_domains[Protein==protein][order(Signature_Description)]

  #### this is for extracting individual domains
  # pos<-sort(unique(unlist(apply(protein_MutS_domains, 1, function(x){
  #   start=as.numeric(x["Start_Location"])
  #   stop=as.numeric(x["Stop_Location"])
  #   start:stop
  # }))))
  #

  ###extracting whole MutS...spanning across all 5 domains
  start=min(protein_MutS_domains$Start_Location)
  stop=max(protein_MutS_domains$Stop_Location)
  pos<-start:stop

  organism<-r["Organism"]
  protein=r["Protein"]
  index<-grep(protein, names(sequences))
  seq<-sequences[[index]]
  seq[pos]
})

names(msh6_domain_sequences)<-paste0(gsub(" ", "_", subresult$Protein),"-", gsub(" ", "_", subresult$Organism))
lapply(msh6_domain_sequences, length)
length(msh6_domain_sequences)
seqinr::write.fasta(msh6_domain_sequences, names=names(msh6_domain_sequences), file.out = "files/msh6_domain_sequences_plants.faa")


# alignment and tree (archived) -------------------------------------------


dom <- readAAStringSet(filepath = "files/msh6_domain_sequences_plants.faa")
alignment <- msa(dom, method = "ClustalW")
alignment_AA <- as(alignment, "AAStringSet")
writeXStringSet(alignment_AA, filepath = "~/Desktop/alignment.fasta")

alignment_matrix <- as.character(alignment)
conservation <- calculateConservation(alignment_matrix)
conservation$POS<-1:nrow(conservation)

alignment_trimmed<-lapply(alignment_AA, function(x) x[conservation[consensus!="-"]$POS])

n<-names(alignment_trimmed)[1]
alignment_matrix <- rbindlist(lapply(names(alignment_trimmed), function(n){
  AA<-unlist(strsplit(as.character(alignment_trimmed[[n]]), split=""))
  if(sum(AA=="-")>50) return(NULL)
  data.table(Organism=n, POS=1:length(AA),AA)
}))


alignment_trimmed <- as(alignment_trimmed, "AAStringSet")

writeXStringSet(alignment_trimmed, filepath = "~/Desktop/alignment_trimmed.fasta")

#read in the data from a format with the aligned sequences (eg msf)
msf.res <- read.alignment(file = "~/Desktop/alignment_trimmed.fasta", format = "fasta")
#dist.alignment and not dna.dist commonly used in tutorials
msf.res.dist.alignment = dist.alignment(msf.res, matrix = "identity", gap=F)
#perform neighbor joining upon the distance set
msf.res.dist.alignment.nj = nj(msf.res.dist.alignment)

x <- as_tibble(msf.res.dist.alignment.nj)
x$plant<-x$label %in% annot_plants$lab
x$other<-x$label %in% annot_other$lab
x$SAR<-x$label %in% annot_sar$lab

x$tudor<-x$label %in% predicted_msh6_blast_annotated[Tudor_domain==T | Tudor_interpro==T]$lab


pdf("~/Desktop/msh6_treeplant.pdf", width=7, height=7)

p0<-ggtree(as.treedata(x), size = 0.1, layout = "circular", color="gray80")+
  theme(plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))+
  geom_tree(size = 0.1, layout = "circular", color="gray80")+
  geom_tippoint(aes(col=tudor), size=0.15)+
  #geom_tippoint(data=as.treedata(x[x$tudor==T,]), shape=21, col="black", size=0.15, fill=NA)+
  scale_color_manual(values = c("gray90","green3"), name="Tudor")+
  geom_tiplab(size=0.5, aes(col=tudor))+
  theme(legend.position = "none")
p0

dev.off()
