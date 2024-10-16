library(data.table)
library(ggplot2)
library(ggtree)
library(ape)
library(parallel)
library(dplyr)
library(tidytree)
library(bio3d)


source("code/organisms_in_node.R")
source("code/read_interproscan.R")

msh6<-fread("data/msh6_organism_table.txt")
potential_orthologs_mismatch<-fread("data/potential_orthologs_mismatch.txt")
predicted_msh6_blast<-fread("data/predicted_msh6_blast.txt")
tblastn_msh6_blast_dist<-fread("data/tblastn_msh6_blast_dist.txt")
predicted_msh6_blast_annotated<-fread("data/predicted_msh6_blast_annotated.txt")
MSH6_human_rev_blastp_all<-fread("data/MSH6_human_rev_blastp_all.txt", fill=T)
interpro<-read_interproscan("data/interpro_all.faa.tsv")
ncbi<-fread("data/ncbi_dataset.tsv")
not_in_PhyloT<-fread("data/not_in_PhyloT.txt", sep="\t", header = F)
CDS_stats2<-fread("data/CDS_stats_all.csv")
CDS_stats<-fread("data/genome_stats.csv")
AA_freq<-fread("data/AA_freq_all.csv")
write(ncbi$`Organism Name`[-1], "data/species.txt")


#fwrite(data.table(msh6[!Organism %in% not_in_PhyloT$V1]$Organism), "data/species_for_tree.csv")
tree<-read.tree("data/phyloT_generated_tree_Aug12_23.txt")
tree$tip.label<-gsub("_"," ",tree$tip.label)

tblastn_msh6_blast_dist$Organism<-ncbi$`Organism Name`[match(tblastn_msh6_blast_dist$GC, ncbi$`Assembly Accession`)]
tblastn_msh6_blast_dist<-tblastn_msh6_blast_dist[,.(Tudortblastn=min(Tudortblastn), PWWPtblastn=min(PWWPtblastn)), by=Organism]
tblastn_msh6_blast_dist<-merge(tblastn_msh6_blast_dist, msh6, by="Organism", all=F)

tblastn_msh6_blast_dist$tipcolor <- ifelse(tblastn_msh6_blast_dist$Tudor_call == "Present", "Tudor",
                                           ifelse(tblastn_msh6_blast_dist$PWWP_call == "Present", "PWWP", "None"))

msh6_good<-tblastn_msh6_blast_dist[!(Tudor_call=="Absent" & Tudortblastn<5000) & !(PWWP_call=="Absent" & PWWPtblastn<5000)]

msh6_good$Genus<-sapply(strsplit(msh6_good$Organism, " "), function(x) x[1])
msh6_good_genera<-msh6_good[!duplicated(Genus)]

plants<-organisms_in_node(tree, msh6_good, c("Viridiplantae"))
landplants<-organisms_in_node(tree, msh6_good, c("Streptophyta"))
prop.table(table(landplants$tipcolor))
meta<-organisms_in_node(tree, msh6_good, "Metazoa")
deuts<-organisms_in_node(tree, msh6_good, "Deuterostomia")
prop.table(table(deuts$tipcolor))
prots<-organisms_in_node(tree, msh6_good, "Protostomia")
chordates<-organisms_in_node(tree, msh6_good, "Chordata")

prop.table(table(prots$tipcolor))

fungi<-organisms_in_node(tree, msh6_good_genera, "Fungi")
other<-organisms_in_node(tree, msh6_good, c("Fungi","Metazoa","Viridiplantae"), exclude = T)
SAR<-organisms_in_node(tree, msh6_good, c("Sar"), exclude = F)

small_tree<-rbindlist(list(fungi, plants, meta, other), fill=T)

#summary stats
uniqueN(ncbi$`Organism Name`)
uniqueN(msh6$Organism)
uniqueN(msh6_good$Organism)

predicted_msh6_blast_annotated<-fread("data/predicted_msh6_blast_annotated.txt")
table(domain=predicted_msh6_blast_annotated$Tudor_domain, interpro=predicted_msh6_blast_annotated$Tudor_interpro)
table(domain=predicted_msh6_blast_annotated$PWWP_domain, interpro=predicted_msh6_blast_annotated$PWWP_interpro)
table(domain=predicted_msh6_blast_annotated$Tudor_interpro, interpro=predicted_msh6_blast_annotated$PWWP_interpro)
table(domain=predicted_msh6_blast_annotated$PWWP_domain, interpro=predicted_msh6_blast_annotated$Tudor_domain)

### Tables
# fwrite(predicted_msh6_blast_annotated, "tables/S2_msh6_orthos_domain_evidence.csv")
# fwrite(msh6_good, "tables/S1_msh6_domains_annotated.csv")
# fwrite(interpro, "tables/S3_interproscan_results.csv")
potential<-fread("data/potential_orthologs_mismatch.txt")
orthos<-fread("data/predicted_msh6_blast_annotated.txt")


