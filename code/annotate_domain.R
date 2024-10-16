#' Annotate Domains in Proteins
#'
#' This function checks whether a given protein was detected by blastp for a specific domain
#' and if the location of the protein in the genome is neighboring a region that shows
#' the domain in another protein from blastp or sequence from tblastn.
#'
#' @param GC Character. The accession from NCBI.
#' @param protein Character. The protein ID for the given accession.
#' @param Q Character. The query for the domain used by blastp and tblastn.
#' @param proteinQ Character. The blastp query originally used to identify the protein.
#' @param interproscan_domain data.table. A data table containing proteins annotated as containing the specified domain.
#'
#' @return
#' If the protein is detected by blastp for the domain, it returns "domain_blast".
#' If the protein is found in the interproscan_domain data table, it returns "interpro".
#' If the protein is neighboring another protein that shows the domain in blastp, it returns a string indicating "neighbor protein dist:" followed by the distance.
#' If the protein is neighboring a sequence that shows the domain in tblastn, it returns a string indicating "dist_tblastn dist:" followed by the distance.
#' If there are no tblastn or blastp hits on the same chromosome, it returns "no tblastn or blastp hits on same chromosome".
#' Otherwise, it returns "error".
#'
#' @examples
#' # Assuming you have the necessary data and files set up:
#' result <- annotate_domain("GC12345", "proteinID123", "domain_query", "protein_query", interproscan_data)
#'
annotate_domain <- function(GC, protein, Q, proteinQ, interproscan_domain) {
  # Check if protein is also found by blast of the domain Q
  evidence = ""
  domain_blast <-
    file.exists(
      paste0(
        "results/blastp/",
        GC,
        "/",
        Q,
        "/proteins/",
        protein,
        ".fasta"
      )
    )

  if (domain_blast)
    evidence=paste(evidence, "domain_blast", sep=";")

  # Check interproscan
  if (protein %in% interproscan_domain$Protein) evidence=paste(evidence, "interpro", sep=";")

  protein_location<-get_protein_location(paste0("results/blastp/", GC, "/", proteinQ, "/gff/",protein,".gff"))

  # Check tblastn for all domains
  tblastn <- read_blast(file = paste0("results/tblastn/", GC, "/", Q, "/tblastn_results.txt"))
  tblastn<-tblastn[sseqid==unique(protein_location$CHROM)]
  if(nrow(tblastn)==0) {
    dist_tblastn<-Inf
  }else dist_tblastn = min(abs(protein_location$START - tblastn$sstart))

  # Check gffs for all domains
  domain_gff_files <- list.files(paste0("results/blastp/", GC, "/", Q, "/gff"),
                                 full.names = T)
  if(length(domain_gff_files)==0) {
    dist<-Inf
  } else {
  all_domain_gff <- rbindlist(lapply(domain_gff_files, get_protein_location))

  other_domain_gff <-
    all_domain_gff[CHROM == unique(protein_location$CHROM) &
                     !grepl(protein, INFO)]
  # if (nrow(other_domain_gff) == 0)
  #   return("no domain gff on same chromosme")

  if(nrow(other_domain_gff)==0) {
    dist<-Inf
  } else dist = min(abs(protein_location$START - other_domain_gff$START))
  }

  #if(is.finite(dist)) dist = min(abs(protein_location$START - other_domain_gff$START))


  evidence=paste(evidence, paste("dist_tblastn dist:", dist_tblastn), sep=";")
  evidence=paste(evidence, paste("neighbor protein dist:", dist), sep=";")
  evidence
  return(evidence)
}
#MSH6 evolution project

