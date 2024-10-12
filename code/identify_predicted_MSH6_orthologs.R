source("github/code/libraries_functions.R")

ncbi<-fread("github/files/ncbi_dataset.tsv", quote="")
system("cat batches_progress/* > done.txt")
complete<-fread("done.txt", header = F)

arabidopsis_names<-fread("data/references/GCA_000001735.2/protein_names.txt", sep="\t", header = F)
arabidopsis_names[, ID := sub(" .*", "", V1)]
arabidopsis_names[, Description := sub("^[^ ]* ", "", V1)]

human_names<-fread("data/references/GCF_000001405.40/protein_names.txt", sep="\t", header = F)
human_names[, ID := sub(" .*", "", V1)]
human_names[, Description := sub("^[^ ]* ", "", V1)]

total <- length(complete$V1)

# Set the number of cores
num_cores <- detectCores()

# Use mclapply to parallelize the computation
msh6_blast_results <- mclapply(1:total, function(i) {
  GC <- complete$V1[i]
  predict_MSH6_orthologs(GC, arabidopsis_names, human_names)
}, mc.cores = num_cores)


predicted_msh6<-rbindlist(lapply(msh6_blast_results, function(x) x$msh6_candidates))
predicted_msh6 <- predicted_msh6[!is.na(protein)]

human_rev_blastp<-rbindlist(lapply(msh6_blast_results, function(x) x$human_rev_blastp))
arabidopsis_rev_blastp<-rbindlist(lapply(msh6_blast_results, function(x) x$arabidopsis_rev_blastp))

potential_orthologs<-rbindlist(lapply(unique(c(human_rev_blastp$sseqid, arabidopsis_rev_blastp$sseqid)), function(ssid){
  h<-human_rev_blastp[sseqid==ssid][1]
  a<-arabidopsis_rev_blastp[sseqid==ssid][1]
  data.table(GC=a$GC, Protein=ssid, human_name=h$name, arabidopsis_name=a$name)
}))

potential_orthologs$human_name<-gsub("(.+) isoform.+","\\1", potential_orthologs$human_name)

potential_orthologs_mismatch<-potential_orthologs[grepl("muts", arabidopsis_name, ignore.case = T) &  grepl("MSH", human_name, ignore.case = T)][!(grepl("MUTS homolog 6", arabidopsis_name, ignore.case = T) &  grepl("MSH6", human_name, ignore.case = T))]

potential_orthologs_mismatch$organism<-ncbi$`Organism Name`[match(potential_orthologs_mismatch$GC, ncbi$`Assembly Accession`)]


fwrite(predicted_msh6, "results/predicted_msh6_blast.txt", sep=" ")
fwrite(human_rev_blastp, "results/MSH6_human_rev_blastp_all.txt", sep=" ")
fwrite(arabidopsis_rev_blastp, "results/MSH6_arabidopsis_rev_blastp_all.txt", sep=" ")
fwrite(potential_orthologs_mismatch, "results/potential_orthologs_mismatch.txt", sep=" ")


