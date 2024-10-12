library(cowplot)
library(phytools)
library(phylolm)
library(polymorphology2)
source("code/organisms_in_node.R")
source("code/read_interproscan.R")

tree<-read.tree("files/phyloT_generated_tree_Aug12_23.txt")
tree$tip.label<-gsub("_"," ",tree$tip.label)

subtree<-extract.clade(tree, "Eukaryota")
ncbi<-fread("files/ncbi_dataset.tsv")
sub_ncbi<-ncbi[`Organism Name` %in% subtree$tip.label]

msh6_good<-fread("tables/S1_msh6_domains_annotated.csv")
CDS_stats<-fread("files/genome_stats.csv")


ncbi_sub<-ncbi[`Organism Name` %in% msh6_good$Organism]
ncbi_sub$tip<-msh6_good$tipcolor[match(ncbi_sub$`Organism Name`, msh6_good$Organism)]
ncbi_sub<-merge(ncbi_sub, CDS_stats, by.x = "Assembly Accession", by.y="GC")
ncbi_sub_summary<-ncbi_sub[,.(Genome_size=mean(`Assembly Stats Total Sequence Length`, na.rm=T),
                              N50=mean(`Assembly Stats Contig N50`, na.rm=T),
                              BUSCO=mean(`Annotation BUSCO Complete`, na.rm=T),
                              Protein_number=mean(`Annotation Count Gene Protein-coding`, na.rm=T),
                              CDS_genome_total=mean(total_CDS_length_genome, na.rm=T),
                              mean_CDS_length=mean(mean_CDS_exon_length, na.rm=T),
                              mean_intron_length=mean(mean_intron_length, na.rm=T)/mean(mean_intron_number, na.rm=T),
                              mean_total_intron_length=mean(mean_intron_length, na.rm=T),
                              prop_intron=mean(prop_intron, na.rm=T),
                              mean_exon_number=mean(mean_exon_number, na.rm=T),
                              gene_lengths=mean(mean_total_length, na.rm=T),
                              mean_total_CDS_length=mean(mean_total_CDS_length, na.rm=T)),
                           by=.(Organism=`Organism Name`, tip )]



ncbi_sub_summary$exonvsintron<-(ncbi_sub_summary$gene_lengths-ncbi_sub_summary$mean_total_intron_length)/ncbi_sub_summary$mean_total_intron_length
ncbi_sub_summary$CDS_percent<-ncbi_sub_summary$CDS_genome_total/ncbi_sub_summary$Genome_size
ncbi_sub_summary$pct_intron_genic<-ncbi_sub_summary$mean_total_intron_length/ncbi_sub_summary$gene_lengths
ncbi_sub_summary$pct_CDS_genic<-ncbi_sub_summary$mean_total_CDS_length/ncbi_sub_summary$gene_lengths

ncbi_sub_summary$genic_percent<-(ncbi_sub_summary$gene_lengths*ncbi_sub_summary$Protein_number)/ncbi_sub_summary$Genome_size
ncbi_sub_summary$nongenic<-ncbi_sub_summary$Genome_size-ncbi_sub_summary$CDS_genome_total
ncbi_sub_summary$pctnongenic<-ncbi_sub_summary$nongenic/ncbi_sub_summary$Genome_size
ncbi_sub_summary$codingnoncoding<-ncbi_sub_summary$CDS_genome_total/ncbi_sub_summary$nongenic

age<-fread("~/Documents/dataset/anage_data.txt")
age$Organism<-paste(age$Genus, age$Species)
ncbi_sub_summary$longevity<-age$`Maximum longevity (yrs)`[match(ncbi_sub_summary$Organism, age$Organism)]

mutation<-fread("~/Documents/mutation_rate_literature_updating3.csv")
mutation$Organism<-mutation$Species
snps<-mutation[TYPE=="snp", .(u=mean(u_mean)), by=Organism]
indels<-mutation[TYPE=="indel", .(u=mean(u_mean)), by=Organism]
ncbi_sub_summary$snp_u<-snps$u[match(ncbi_sub_summary$Organism, snps$Organism)]
ncbi_sub_summary$indel_u<-indels$u[match(ncbi_sub_summary$Organism, indels$Organism)]

pop<-fread("~/Documents/combined_data.tsv")
pop$Organism<-pop$species
ncbi_sub_summary$body_mass<-10^pop$pred_log10_body_mass[match(ncbi_sub_summary$Organism, pop$Organism)]
ncbi_sub_summary$Ne<-10^pop$pred_log10_N[match(ncbi_sub_summary$Organism, pop$Organism)]

#confirming concordance of assembly size with C-value
Marino_table1<-fread("local_data/Marino_TableS1.tsv")
cor.test(log10(Marino_table1$Cvalue), log10(Marino_table1$Assembly_size))

# Fetch the organisms in each node
Hexapoda <- organisms_in_node(tree, msh6_good, c("Hexapoda"))
Chlorophyta <- organisms_in_node(tree, msh6_good, c("Chlorophyta"))
Arthropoda <- organisms_in_node(tree, msh6_good, c("Arthropoda"))
Protostomia <- organisms_in_node(tree, msh6_good, c("Protostomia"))
Sar <- organisms_in_node(tree, msh6_good, c("Sar"))
Chordata <- organisms_in_node(tree, msh6_good, c("Chordata"))
Mammalia <- organisms_in_node(tree, msh6_good, c("Mammalia"))
Fungi <- organisms_in_node(tree, msh6_good, "Fungi")
Cnidaria <- organisms_in_node(tree, msh6_good, "Cnidaria")
Spiralia <- organisms_in_node(tree, msh6_good, "Spiralia")
Viridiplantae <- organisms_in_node(tree, msh6_good, c("Viridiplantae"))
Metazoa <- organisms_in_node(tree, msh6_good, "Metazoa")
Deuterostomia <- organisms_in_node(tree, msh6_good, "Deuterostomia")

# Add columns to ncbi_sub_summary
ncbi_sub_summary$Hexapoda <- ncbi_sub_summary$Organism %in% Hexapoda$Organism
ncbi_sub_summary$Chlorophyta <- ncbi_sub_summary$Organism %in% Chlorophyta$Organism
ncbi_sub_summary$Arthropoda <- ncbi_sub_summary$Organism %in% Arthropoda$Organism
ncbi_sub_summary$Protostomia <- ncbi_sub_summary$Organism %in% Protostomia$Organism
ncbi_sub_summary$Sar <- ncbi_sub_summary$Organism %in% Sar$Organism
ncbi_sub_summary$Chordata <- ncbi_sub_summary$Organism %in% Chordata$Organism
ncbi_sub_summary$Mammalia <- ncbi_sub_summary$Organism %in% Mammalia$Organism
ncbi_sub_summary$Fungi <- ncbi_sub_summary$Organism %in% Fungi$Organism
ncbi_sub_summary$Cnidaria <- ncbi_sub_summary$Organism %in% Cnidaria$Organism
ncbi_sub_summary$Spiralia <- ncbi_sub_summary$Organism %in% Spiralia$Organism
ncbi_sub_summary$Viridiplantae <- ncbi_sub_summary$Organism %in% Viridiplantae$Organism
ncbi_sub_summary$Metazoa <- ncbi_sub_summary$Organism %in% Metazoa$Organism
ncbi_sub_summary$Deuterostomia <- ncbi_sub_summary$Organism %in% Deuterostomia$Organism

## species from OrthoDB
species<-fread("files/orthoDB/odb11v0_species.tab.gz")
setnames(species, old = c("V1", "V2", "V3", "V4", "V5", "V6", "V7"),
         new = c("NCBI_tax_id", "OrthoDB_org_id", "Scientific_name",
                 "Genome_assembly_id", "Total_clustered_genes",
                 "Total_OGs_participated", "Mapping_type"))

species<-species[Genome_assembly_id %in% ncbi_sub$`Assembly Accession`]
species$Organism<-ncbi$`Organism Name`[match(species$Genome_assembly_id, ncbi$`Assembly Accession`)]

ncbi_sub_summary$OrthoDB_avail <- ncbi_sub_summary$Organism %in% species$Organism

# tree from timetree
tree_time<-read.tree("files/time_tree_may4.nwk") #from http://www.timetree.org/

tree_time$tip.label<-gsub("_"," ",tree_time$tip.label)
tree_time <- multi2di(tree_time)
tree_time$edge.length <- tree_time$edge.length + 1e-5

ncbi_sub_summary$TimeTree_avail <- ncbi_sub_summary$Organism %in% tree_time$tip.label
ncbi_sub_summary$MSH6_reader<-ncbi_sub_summary$tip
fwrite(ncbi_sub_summary,"tables/S4_species_traits.csv")

