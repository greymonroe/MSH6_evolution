### Analysis of orthoDB (accessed May 10, 2024)
# this code examines matrix of orthogroups in each species
# to test for correlations between orthogroup PAV and MSH6 histone reader
####

# Functions and libraries -------------------------------------------------
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
i_chitest<-function(i){

  message(i)
  Ortho_ID<-colnames(wide)[i]
  cand<-wide[[i]]
  cand<-as.numeric(cand>0)
  names(cand)<-wide$Organism
  traits_time$p<-cand[match(traits_time$Organism, names(cand))]
  traits_time<-traits_time[!is.na(p),]
  sum<-sum(traits_time$p)
  if(sum<20|uniqueN(traits_time$p)<2) return(NULL)
  trait<-as.numeric(traits_time$p>0)
  names(trait)<-traits_time$Organism
  reader<-traits_time$tip!="None"
  names(reader)<-traits_time$Organism
  chitest<- fisher.test(as.matrix(table(reader, trait)))
  as.matrix(table(reader, trait))
  data.table(i, sum, Ortho_ID, p=chitest$p.value, or=chitest$estimate)
}
i_fitPagel<-function(i){

  Ortho_ID<-colnames(wide)[i]
  message(paste(i, Ortho_ID))
  cand<-wide[[i]]
  cand<-as.numeric(cand>0)
  names(cand)<-wide$Organism
  traits_time$p<-cand[match(traits_time$Organism, names(cand))]
  traits_time<-traits_time[!is.na(p),]
  sum<-sum(traits_time$p)
  if(sum<20|uniqueN(traits_time$p)<2) return(NULL)
  trait<-as.numeric(traits_time$p>0)
  names(trait)<-traits_time$Organism
  reader<-traits_time$tip!="None"
  names(reader)<-traits_time$Organism
  chitest<- fisher.test(as.matrix(table(reader, trait)))
  as.matrix(table(reader, trait))
  tree_time_sub <- drop.tip(tree_time, tree_time$tip.label[!tree_time$tip.label %in% (traits_time$Organism)])

  fit.ape<-fitPagel(tree_time_sub,reader,trait, model = "ARD")
  return(fit.ape)
}
i_tree<-function(i){
  message(i)
  Ortho_ID<-colnames(wide)[i]
  cand<-wide[[i]]
  cand<-as.numeric(cand>0)
  names(cand)<-wide$Organism
  traits_time$p<-cand[match(traits_time$Organism, names(cand))]
  traits_time<-traits_time[!is.na(p),]
  sum<-sum(traits_time$p)
  if(sum<20|uniqueN(traits_time$p)<2) return(NULL)
  trait<-as.numeric(traits_time$p>0)
  names(trait)<-traits_time$Organism
  reader<-traits_time$tip!="None"
  names(reader)<-traits_time$Organism
  chitest<- fisher.test(as.matrix(table(reader, trait)))
  as.matrix(table(reader, trait))
  tree_time_sub <- drop.tip(tree_time, tree_time$tip.label[!tree_time$tip.label %in% (traits_time$Organism)])

  x <- as_tibble(tree_time_sub)
  x$pwwp<-sapply(x$label, function(n){
    if(n %in% tree_time_sub$node.label){
      clade<-extract.clade(tree_time_sub, n)
      sub_traits_time<-traits_time[Organism %in% clade$tip.label]
      pwwp<-sum(sub_traits_time$tip!="None")
    } else {
      sub_traits_time<-traits_time[Organism %in% n]
      pwwp<-sum(sub_traits_time$tip!="None")
    }
    return(pwwp)
  })

  x$ortho<-sapply(x$label, function(n){
    if(n %in% tree_time_sub$node.label){
      clade<-extract.clade(tree_time_sub, n)
      sub_traits_time<-traits_time[Organism %in% clade$tip.label]
      ortho<-sum(sub_traits_time$p==1)
    } else {
      sub_traits_time<-traits_time[Organism %in% n]
      ortho<-sum(sub_traits_time$p==1)
    }
    return(ortho)
  })

  x$N<-sapply(x$label, function(n){
    if(n %in% tree_time_sub$node.label){
      clade<-extract.clade(tree_time_sub, n)
      sub_traits_time<-traits_time[Organism %in% clade$tip.label]
      N<-nrow(sub_traits_time)
    } else {
      sub_traits_time<-traits_time[Organism %in% n]
      N<-nrow(sub_traits_time)
    }
    return(N)
  })


  tree_small<-drop.tip(tree_time_sub, tree_time_sub$tip.label[!tree_time_sub$tip.label %in% traits_time$Organism])
  y <- as_tibble(tree_small)

  d<-x %>% dplyr::select(label, pwwp, ortho, N)
  z <- left_join(y, unique(d), by = 'label')

  return(z)
}
plot_i_tree<-function(z){
  hjustment=-1
  p1<-ggtree(as.treedata(z), size = 0.1, layout = "circular", color="gray80")+
    theme(plot.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"))+
    geom_tiplab(size=0.2, alpha=1, hjust = hjustment)

  p2<-ggtree(as.treedata(z), size = 0.4, layout = "circular", aes(col=pwwp/N, alpha=pwwp/N))+
    scale_color_continuous(high='orange2', low='gray70', guide="none")+
    scale_alpha_continuous(range = c(0,1), guide="none")+
    theme(plot.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"))+
    geom_tiplab(size=0.2, alpha=1, hjust = hjustment)

  p3<-ggtree(as.treedata(z), size = 0.4, layout = "circular", aes(col=ortho/N, alpha=ortho/N))+
    scale_color_continuous(high='purple', low='gray70', guide="none")+
    scale_alpha_continuous(range = c(0,1), guide="none")+
    theme(plot.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"))+
    geom_tiplab(size=0.2, alpha=1, hjust = hjustment)

  p4<-ggtree(as.treedata(z), size = 0.1, layout = "circular", color="gray80", alpha=0)+
    theme(plot.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"))+
    geom_tiplab(size=0.2, alpha=1, hjust = hjustment)


  p5<-ggtree(as.treedata(z), size = 0.1, layout = "circular", color="gray80")+
    theme(plot.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"))+
    geom_tippoint(aes(col=factor(ortho)), size=0.5)+
    scale_color_manual(values=c("white","purple"), guide="none")+
    geom_tiplab(size=0.2, alpha=1, hjust = hjustment)


  list(p1, p2, p3, p4, p5)
}
Multiple_threshold_test <- function(allScore) {
  return(allScore == 1)
}
library(topGO)
library(polymorphology2)


# Load data ---------------------------------------------------------------

# tree from timetree
tree_time<-read.tree("files/time_tree_may4.nwk") #from http://www.timetree.org/

tree_time$tip.label<-gsub("_"," ",tree_time$tip.label)
tree_time <- multi2di(tree_time)
tree_time$edge.length <- tree_time$edge.length + 1e-5

tree<-read.tree("files/phyloT_generated_tree_Aug12_23.txt")
tree$tip.label<-gsub("_"," ",tree$tip.label)
ncbi<-fread("files/ncbi_dataset.tsv")

ncbi_sub<-ncbi[`Organism Name` %in% tree$tip.label]

species<-fread("~/Desktop/odb11v0_species.tab.gz")
setnames(species, old = c("V1", "V2", "V3", "V4", "V5", "V6", "V7"),
         new = c("NCBI_tax_id", "OrthoDB_org_id", "Scientific_name",
                 "Genome_assembly_id", "Total_clustered_genes",
                 "Total_OGs_participated", "Mapping_type"))

species<-species[Genome_assembly_id %in% ncbi_sub$`Assembly Accession`]
species$Organism<-ncbi$`Organism Name`[match(species$Genome_assembly_id, ncbi$`Assembly Accession`)]

wide<-fread("files/orthoDB/OrthoDB_species_wide.txt")
wide$Organism<-species$Organism[match(wide$OrthoDB_individual_organism_id, species$OrthoDB_org_id)]

OG_info<-fread("files/orthoDB/odb11v0_OGs.tab.gz")
colnames(OG_info)<-c("Ortho_ID","Tax","Name")
OG_info<-OG_info[Ortho_ID %in% colnames(wide)]

OG_Xref<-fread("files/orthoDB/odb11v0_OG_xrefs.tab.gz")
colnames(OG_Xref)<-c("Ortho_ID","exDB","exDB_ID","N")
OG_Xref<-OG_Xref[Ortho_ID %in% colnames(wide)]

# run chitests ---------------------------------------------------------

traits<-fread("tables/S4_species_traits.csv")
traits_time<-traits
#Protostomes only
traits_time<-data.table(ncbi_sub_summary[Organism %in% Protostomia$Organism])
#All species
#traits_time<-ncbi_sub_summary

#### below takes ~2.5hrs to run
# t<-Sys.time()
# chitests<-rbindlist(lapply(2:(ncol(wide)-1), function(i){
#   i_chitest(i)
# }))
# Sys.time()-t
# fwrite(chitests, "files/orthoDB/chitests_all_species.csv")
#fwrite(chitests, "files/orthoDB/chitests_protostomes.csv")

# run GO enrichment tests -------------------------------------------------

## Analysis scope. The Protostome co-evolution analysis reported in the manuscript
## uses chitests_protostomes.csv. Swap the comments (and only then) for the
## all-species run -- previously BOTH lines were live, so the protostome table was
## silently overwritten by the all-species table.
chitests<-fread("files/orthoDB/chitests_protostomes.csv")
#chitests<-fread("files/orthoDB/chitests_all_species.csv")

chitests$Name<-OG_info$Name[match(chitests$Ortho_ID, OG_info$Ortho_ID)]
chitests$p.adjust<-p.adjust(chitests$p)

chitests_Xref<-merge(chitests, OG_Xref, by="Ortho_ID")
chitests_Xref_Xref_bio<-chitests_Xref[exDB=="biological_process"]


## gene2GO: ONE entry per ortholog, holding ALL of that ortholog's GO BP terms.
## Do not build this with as.list() on the long-format table: that yields duplicate
## list names, and annFUN.gene2GO subsets via gene2GO[intersect(names, feasible)],
## which silently returns only the FIRST GO term per ortholog. That collapses the
## annotation (883 -> 732 genes on GO:0006281) and inflates the enrichment p-value
## by many orders of magnitude.
geneID2GO<-split(chitests_Xref_Xref_bio$exDB_ID, chitests_Xref_Xref_bio$Ortho_ID)

## geneList likewise needs one entry per ortholog, not one per (ortholog, GO term) row.
chitests_bio_unique<-unique(chitests_Xref_Xref_bio[,.(Ortho_ID, or, p)])
geneList<-as.numeric(chitests_bio_unique$or>1 & chitests_bio_unique$p<0.05)
names(geneList)<-chitests_bio_unique$Ortho_ID


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel = Multiple_threshold_test,
              annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allResOGs <- GenTable(GOdata, classicFisher = resultClassic, weight01 = resultFisher,
                      ranksOf = "classicFisher", topNodes = 20)
print(allResOGs)
#fwrite(allResOGs, file="files/orthoDB/GO_enrich_protostomes.csv")
#fwrite(allResOGs, file="files/orthoDB/GO_enrich_all.csv")

allGO = genesInTerm(GOdata)

repair_candidates<-chitests_Xref_Xref_bio[Ortho_ID %in% c(allGO$`GO:0006281`) & or>1 & p<0.05]
other_candidates<-chitests_Xref_Xref_bio[Ortho_ID %in% c(allGO$`GO:0006334`, allGO$`GO:0016567`) & or>1 & p<0.05]
unique(repair_candidates$Name)
unique(other_candidates$Name)

# run Pagel phylogenetic co-evolution ----------------------------------------

Is<-unique(c(repair_candidates$i, other_candidates$i))

repair_Pagel<-lapply(Is, function(i){
  t<-Sys.time()
  if(!file.exists(paste0("files/orthoDB/Pagel_RData/protosome_",i,".RData"))){
  Pagel<-i_fitPagel(i)
  save(Pagel, file=paste0("files/orthoDB/Pagel_RData/protosome_",i,".RData"))
  message(Sys.time()-t)
  return(Pagel)
  } else{
    load(paste0("files/orthoDB/Pagel_RData/protosome_",i,".RData"))
    return(Pagel)}
})
names(repair_Pagel)<-Is

repair_Pagel_P<-rbindlist(lapply(Is, function(i){
  load(paste0("files/orthoDB/Pagel_RData/protosome_",i,".RData"))
  P<-Pagel$P[1]
  data.table(i, P)
}))
repair_Pagel_P<-merge(repair_Pagel_P, chitests_Xref_Xref_bio, by="i")
## BH, matching the "Pagel P (BH-adjusted)" column header in Table S5.
## p.adjust() defaults to "holm" -- leaving this bare mislabels the reported values.
repair_Pagel_P$P.adjust<-p.adjust(repair_Pagel_P$P, method="BH")
repair_Pagel_P$repair<-repair_Pagel_P$Name %in% repair_candidates$Name
repair_Pagel_P_sig<-repair_Pagel_P[P<0.05]

# Plot tree ---------------------------------------------------------------

pdf("Figures/time_treecircular_ortho.pdf", width=3, height=3)
i=373208
z<-i_tree(i)
plot_i_tree(z)
#chitests[i==i2]$Name
# chitests[i==i2]$p
# load(paste0("files/orthoDB/Pagel_RData/protosome_",i,".RData"))
# Pagel
dev.off()


# SET domains tests -------------------------------------------------------
chitests$Name<-OG_info$Name[match(chitests$Ortho_ID, OG_info$Ortho_ID)]
OG_Xref$Name<-OG_info$Name[match(OG_Xref$Ortho_ID, OG_info$Ortho_ID)]
OG_Xref_reader<-OG_Xref[exDB_ID=="IPR001025"]
unique(OG_Xref_reader$Name)
OG_Xref_reader[Name=="RSC complex subunit"]
chitests[Name=="RSC complex subunit"]
OG_Xref_repair<-OG_Xref[exDB_ID=="GO:0007062"]
OG_Xref_reader_repair<-chitests[Ortho_ID %in% OG_Xref_repair$Ortho_ID &Ortho_ID %in% OG_Xref_reader$Ortho_ID]
unique(OG_Xref_reader_repair$Name)
chitests_Xref_PWWP_repair

chitests_sig_SET2<-chitests_sig[grepl("SET domain containing 2", Name, ignore.case = T)]
chitests_Xref_interpro<-chitests_Xref[exDB=="interpro_domains"]
#IPR001214 = SET domain
ggplot(chitests_Xref_interpro, aes(x=log(or), group=exDB_ID=="IPR001214", fill=exDB_ID=="IPR001214"))+
  geom_density(alpha=0.5)

chitests_Xref_SET<-chitests_Xref_interpro[exDB_ID=="IPR001214"]
View(chitests_Xref_SET)
chitests_sig<-chitests
SETtest<-chitests_sig[,.(enrich=sum(or<1), N=.N), by=.(SET=i  %in% chitests_Xref_SET$i)]
chisq.test(SETtest[,2:3])

chitests_sig<-chitests

SETtest<-chitests[p.adjust<0.05,.(SET=sum(i  %in% chitests_Xref_SET$i), N=.N), by=.(or=as.numeric(cut2(or, g=100)))]
ggplot(SETtest, aes(x=or, y=SET/N))+
  geom_point()

t.test(-log10(chitests_Xref_SET$p)~chitests_Xref_SET$exDB_ID=="IPR001214")
SETtest<-chitests_Xref_SET[,.(enrich=sum(or<1), N=.N), by=.(exDB_ID)]
SETtest$rat<-SETtest$enrich/SETtest$N


# Integrate chisqtests ----------------------------------------------------

chitests_prot<-fread("files/orthoDB/chitests_protostomes.csv")
chitests_all<-fread("files/orthoDB/chitests_all_species.csv")
chitests<-merge(chitests_prot, chitests_all, by="i")

allResOGs_prot<-fread("files/orthoDB/GO_enrich_protostomes.csv")
allResOGs_prot$oe<-log(allResOGs_prot$Significant/allResOGs_prot$Expected)
allResOGs_all<-fread("files/orthoDB/GO_enrich_all.csv")
allResOGs_all$oe<-log(allResOGs_all$Significant/allResOGs_all$Expected)

allResOGs_merge<-merge(allResOGs_prot, allResOGs_all, by="GO.ID")
ggplot(allResOGs_merge, aes(x=(oe.x), y=oe.y))+
  geom_point()
View(allResOGs_merge)

# screen orthogroups for reader/repairs -----------------------------------


## PWWP: IPR000313
## Tudor: IPR002999
## Zinc finger, RING/FYVE/PHD-type: IPR013083
## Zinc Fingers: IPR013083
OG_Xref[exDB=="interpro_domains"]
OG_Xref_reader<-OG_Xref[Ortho_ID %in% OG_Xref[exDB_ID=="IPR011011"]$Ortho_ID]
OG_Xref_reader$Name<-OG_info$Name[match(OG_Xref_reader$Ortho_ID, OG_info$Ortho_ID)]


OG_Xref_reader_repair<-OG_Xref[Ortho_ID %in% OG_Xref[exDB_ID=="IPR002999"]$Ortho_ID & exDB_ID=="GO:0006281"]
OG_Xref_reader_repair$Name<-OG_info$Name[match(OG_Xref_reader_repair$Ortho_ID, OG_info$Ortho_ID)]
OG_Xref_reader_repair$Ortho_ID
View(wide[,.(`1023467at32523`, Organism)])

OG_Xref_repair<-OG_Xref[Ortho_ID %in% OG_Xref[exDB_ID=="GO:0006281"]$Ortho_ID]
OG_Xref_repair$Name<-OG_info$Name[match(OG_Xref_repair$Ortho_ID, OG_info$Ortho_ID)]




# Archive -----------------------------------------------------------------


### Data prep (filtering to only look at species in phylogeny)
# awk 'NR == FNR { ids[$1] = 1; next } $2 in ids' filtered_odb11v0_genes.OrgID.txt odb11v0_genes.tab > filtered_odb11v0_genes.tab

# IDfile<-fread("~/Desktop/odb11v0_genes.OrgID.txt", header = F)
# IDfile$row<-1:nrow(IDfile)
#
# cat filtered_odb11v0_genes.tab | cut -f2 > allfiltered_odb11v0_genes.OrgID.txt
# cat filtered_odb11v0_genes.tab | cut -f1 > filtered_odb11v0_genes.geneID.txt
#
# gzcat odb11v0_OG2genes.tab.gz | awk 'NR == FNR { ids[$1] = 1; next } $2 in ids' filtered_odb11v0_genes.geneID.txt - | gzip > filtered_odb11v0_OG2genes.tab.gz
#
# #rm("merged")
# OGs<-fread("~/Desktop/filtered_odb11v0_OG2genes.tab.gz", header=F)
# colnames(OGs)<-c("OrthoDB_orthogroup_id","OrthoDB_unique_gene_id")
# t<-Sys.time()
# genes<-fread("~/Desktop/odb11v0_genes.tab", select = c(1, 2))
# Sys.time()-t
# colnames(genes)<-c("OrthoDB_unique_gene_id","OrthoDB_individual_organism_id")
# #genes$orthogroup<-OGs$`0at1028384`[match(genes$OrthoDB_unique_gene_id, OGs$`1173701_1:00197a`)]
# merged<-merge(genes, OGs, by="OrthoDB_unique_gene_id")
# #rm("genes")
# #rm("OGs")
# fwrite(merged, "~/Desktop/OrthoDB_msh6.txt")
# Ortho_table<-data.table(table(merged$OrthoDB_orthogroup_id))
# hist(log10(Ortho_table$N))
# fwrite(Ortho_table, "~/Desktop/OrthoDB_table.txt")
# Ortho_table2<-Ortho_table[N>50] #only consider orthogroups with at least 50 proteins
# merged<-merged[OrthoDB_orthogroup_id %in% Ortho_table2$V1]
# wide<-dcast(merged, formula = OrthoDB_individual_organism_id~OrthoDB_orthogroup_id)
# fwrite(wide, "~/Desktop/OrthoDB_species_wide.txt")


