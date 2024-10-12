
source("github/code/libraries_functions.R")


interpro<-read_interproscan("results/predicted_msh6/all.faa.tsv")
GCs<-fread("done.txt", header = F)
ncbi<-fread("github/files/ncbi_dataset.tsv", quote="")
predicted_msh6<-fread("results/predicted_msh6_blast.txt")

interproscan_tudor<-interpro[Analysis!="SUPERFAMILY" & (grepl("tudor", InterPro_Description, ignore.case = T) | grepl("tudor", Signature_Description, ignore.case = T))]
interproscan_pwwp<-interpro[Analysis!="SUPERFAMILY" & (grepl("pwwp", InterPro_Description, ignore.case = T) | grepl("pwwp", Signature_Description, ignore.case = T))]

MutS<-interpro[,.(MutS=sum(grepl("Muts", InterPro_Description, ignore.case = T))>0), by=Protein]

predicted_msh6$Organism<-ncbi$`Organism Name`[match(predicted_msh6$GC, ncbi$`Assembly Accession`)]

predicted_msh6$MutS=predicted_msh6$protein %in% MutS$Protein

num_cores <- 32

# Parallelize the apply function for Tudor
predicted_msh6$Tudor <- mclapply(1:nrow(predicted_msh6), function(i) {
  r <- predicted_msh6[i,]
  annotate_domain(r$GC, r$protein, Q = "arabidopsis_msh6_tudor", proteinQ = "arabidopsis_msh6", interproscan_domain = interproscan_tudor)
}, mc.cores = num_cores)

# Parallelize the apply function for PWWP

predicted_msh6$PWWP <- mclapply(1:nrow(predicted_msh6), function(i) {
  r <- predicted_msh6[i,]
  annotate_domain(r$GC, r$protein, Q = "human_msh6_pwwp", proteinQ = "human_msh6", interproscan_domain = interproscan_pwwp)
}, mc.cores = num_cores)

predicted_msh6$tudor_tblastn<-as.numeric(gsub(".+dist_tblastn dist: (.+);neighbor.+","\\1", predicted_msh6$Tudor))

predicted_msh6$tudor_neighbor<-as.numeric(gsub(".+neighbor protein dist: (.+)","\\1", predicted_msh6$Tudor))

predicted_msh6$PWWP_tblastn<-as.numeric(gsub(".+dist_tblastn dist: (.+);neighbor.+","\\1", predicted_msh6$PWWP))

predicted_msh6$PWWP_neighbor<-as.numeric(gsub(".+neighbor protein dist: (.+)","\\1", predicted_msh6$PWWP))

predicted_msh6$Tudor_domain<-grepl("domain_blast", predicted_msh6$Tudor)
predicted_msh6$Tudor_interpro<-grepl("interpro", predicted_msh6$Tudor)

predicted_msh6$PWWP_domain<-grepl("domain_blast", predicted_msh6$PWWP)
predicted_msh6$PWWP_interpro<-grepl("interpro", predicted_msh6$PWWP)



# archive non-parallel version
# predicted_msh6$Tudor<-apply(predicted_msh6, 1, function(r){
#   annotate_domain(r["GC"], r["protein"], "arabidopsis_msh6_tudor","arabidopsis_msh6", interproscan_tudor)
# })
#
# predicted_msh6$PWWP<-apply(predicted_msh6, 1, function(r){
#   annotate_domain(r["GC"], r["protein"], Q = "human_msh6_pwwp",proteinQ = "human_msh6", interproscan_domain = interproscan_pwwp)
# })

# Parallelize the sapply function for PWWPtblastn
GCs$PWWPtblastn <- unlist(mclapply(GCs$V1, function(GC) {
  dist <- min(check_tblastn(GC, "human_msh6_pwwp", proteinQ = "human_msh6")$dist_to_domain)
  return(dist)
}, mc.cores = num_cores))

# Parallelize the sapply function for Tudortblastn
GCs$Tudortblastn <- unlist(mclapply(GCs$V1, function(GC) {
  dist <- min(check_tblastn(GC, "arabidopsis_msh6_tudor", proteinQ = "arabidopsis_msh6")$dist_to_domain)
  return(dist)
}, mc.cores = num_cores))
colnames(GCs)[1]<-"GC"
# archive non-parallel version

# GCs$PWWPtblastn<-sapply(GCs$V1, function(GC){
#   dist<-min(check_tblastn(GC, "human_msh6_pwwp", proteinQ = "human_msh6")$dist_to_domain)
# })
# GCs$Tudortblastn<-sapply(GCs$V1, function(GC){
#   dist<-min(check_tblastn(GC, "arabidopsis_msh6_tudor", proteinQ = "arabidopsis_msh6")$dist_to_domain)
# })




fwrite(predicted_msh6, "results/predicted_msh6_blast_annotated.txt", sep=" ")
msh6<-make_MSH6_organism_table("results/predicted_msh6_blast_annotated.txt", 10000)
fwrite(GCs, "results/tblastn_msh6_blast_dist.txt", sep=" ")
fwrite(msh6, "results/msh6_organism_table.txt", sep=" ")



