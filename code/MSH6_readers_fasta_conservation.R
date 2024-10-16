library(seqinr)
library(data.table)
library(msa)
library(Biostrings)
source("code/read_blast.R")


# functions ---------------------------------------------------------------

shannon_entropy <- function(char_vector) {
  # Calculate the frequency of each character
  freq <- table(char_vector)

  # Calculate the probability of each character
  prob <- freq / sum(freq)

  # Compute the entropy
  entropy <- -sum(prob * log2(prob))

  return(entropy)
}


calculateConservation <- function(alignment) {
  alignedSeqs <- as.character(alignment)
  n <- length(alignedSeqs)
  l <- nchar(alignedSeqs[1])

  conservation <- numeric(l)
  aa <- character(l)
  entropy <- numeric(l)

  # Create a list to store amino acid frequencies for each position
  aa_frequencies <- vector("list", l)

  for (i in 1:l) {
    position <- substr(alignedSeqs, i, i)
    freq_table <- table(position)
    mostCommon <- names(sort(freq_table, decreasing = TRUE))[1]
    aa[i] <- mostCommon
    conservation[i] <- sum(position == mostCommon) / n
    entropy[i] <- shannon_entropy(position)

    freq_table2<-table(position)
    # Store the frequencies for each amino acid at this position as a data.table
    aa_frequencies[[i]] <- as.data.table(as.list(freq_table2 / n))
  }

  # Convert the list of frequencies to a data.table
  freq_dt <- rbindlist(aa_frequencies, fill = TRUE, use.names = TRUE)

  # Combine the results
  result <- data.table(consensus = aa, conservation = conservation, entropy = entropy)
  result <- cbind(result, freq_dt)

  return(result)
}

align_conserv <- function(fasta_file, domain, method, ref_sequence_name) {

  # Define the output filepath
  output_filepath <- paste0("data/", method, "_msh6_", domain, ".fasta")

    sequences <- readAAStringSet(filepath = fasta_file)
    sequences<-unique(sequences)
    alignment <- msa(sequences, method = method)
    alignment_matrix <- as.character(alignment)
    alignment_AAStringSet <- AAStringSet(alignment_matrix)
    writeXStringSet(alignment_AAStringSet, filepath = output_filepath)

  conservation <- calculateConservation(alignment_matrix)

  # Extract the sequence by its name
  extracted_sequence <- unlist(strsplit(alignment_matrix[ref_sequence_name], split=""))

  # Create a data table with amino acid and conservation values
  conservation_dt <- cbind(data.table(refAA = extracted_sequence), conservation)
  conservation_dt <- conservation_dt
  conservation_dt$POS <- 1:nrow(conservation_dt)

  # Write the data table to a CSV file
  output_csv <- paste0("data/", method, "_msh6_", domain, "_conservation.csv")
  fwrite(conservation_dt, output_csv)

  return(conservation_dt)
}


# run code below

predicted_msh6_blast_annotated<-fread("data/predicted_msh6_blast_annotated.txt")

# Initialize lists to store the domains
tudors <- list()
pwwps <- list()

ncbi<-fread("data/ncbi_dataset.tsv")
# Initialize the progress bar
#pb <- txtProgressBar(min = 0, max = nrow(predicted_msh6_blast_annotated), style = 3)

system("cat batches_progress/* > done.txt")
complete<-fread("done.txt", header = F)
GCs <- complete$V1

predicted_msh6_blast_annotated<-predicted_msh6_blast_annotated[GC %in% GCs]
for (i in 1:nrow(predicted_msh6_blast_annotated)) {

  row <- predicted_msh6_blast_annotated[i, ]
  GC <- row$GC
  protein <- row$protein

  # Extract organism name based on Assembly Accession
  organism_name <- ncbi$`Organism Name`[ncbi$`Assembly Accession` == GC]
message(paste("i", i))
  # Check for Tudor domain
  if (row$Tudor_domain == TRUE) {
    fasta <- read.fasta(paste0("results/blastp/", GC, "/arabidopsis_msh6_tudor/proteins/", protein, ".fasta"))
    blast <- read_blast(paste0("results/blastp/", GC, "/arabidopsis_msh6_tudor/blastp_results.txt"))
    protein_blast <- blast[blast$sseqid == protein, ][1]
    min<-(protein_blast$sstart-5)
    if(min<0) min<-1
    max<-(protein_blast$send+5)
    tudor <- fasta[[1]][min:max]
    names(tudor) <- paste(organism_name, GC, protein)
    tudors[[paste(organism_name, GC, protein)]] <- tudor
  }

  # Check for PWWP domain
  if (row$PWWP_domain == TRUE) {
    fasta <- read.fasta(paste0("results/blastp/", GC, "/human_msh6_pwwp/proteins/", protein, ".fasta"))
    blast <- read_blast(paste0("results/blastp/", GC, "/human_msh6_pwwp/blastp_results.txt"))
    protein_blast <- blast[blast$sseqid == protein, ][1]
    min<-(protein_blast$sstart-5)
    if(min<0) min<-1
    max<-(protein_blast$send+5)
    pwwp <- fasta[[1]][min:max]
    names(pwwp) <- paste(organism_name, GC, protein)
    pwwps[[paste(organism_name, GC, protein)]] <- pwwp
  }

  # Update the progress bar
  #setTxtProgressBar(pb, i)
}

# Close the progress bar
#close(pb)


# Convert lists to fasta objects
# write.fasta(tudors, file = "data/msh6_tudors.fasta", names = names(tudors))
# write.fasta(pwwps, file = "data/msh6_pwwps.fasta", names = names(pwwps))

tudors<-read.fasta("data/msh6_tudors.fasta")
tudors_lengths<-unlist(lapply(tudors, length))

pwwps<-read.fasta("data/msh6_pwwps.fasta")
pwwps_lengths<-unlist(lapply(pwwps, length))


methods <- c("Muscle")

# Loop through each method
for (method in methods) {

  # Apply align_conserv for msh6_pwwps.fasta
  align_conserv("results/msh6_pwwps.fasta",
                "pwwps",
                method = method,
                "Homo sapiens GCF_000001405.40 NP_001394291.1")

  # Apply align_conserv for msh6_tudors.fasta
  align_conserv("results/msh6_tudors.fasta",
                "tudors",
                method = method,
                "Arabidopsis thaliana GCF_000001735.4 NP_192116.1")
}

