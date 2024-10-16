#' Check TBLASTN Results for Protein and Domain Queries
#'
#' This function processes TBLASTN results for protein and domain queries. It filters the protein query results based on length and start position criteria, and then calculates the minimum distance to the closest domain query result for each protein query result.
#'
#' @param GC A character string specifying the genome context or identifier.
#' @param domainQ A character string specifying the domain query.
#' @param proteinQ A character string specifying the protein query.
#'
#' @return A data.table containing the protein query TBLASTN results with an additional column indicating the distance to the closest domain query result. If no domain query result is found for a protein query result, the distance is set to 'Inf'.
#'
#' @examples
#' # Assuming you have set up the necessary files and directory structure:
#' # check_tblastn("someGC", "someDomainQ", "someProteinQ")
#'
check_tblastn <- function(GC, domainQ, proteinQ) {

  # Read TBLASTN results
  # Check if proteinQ_tblastn file exists
  proteinQ_file_path <- paste0("results/tblastn/", GC, "/", proteinQ, "/tblastn_results.txt")
  if (!file.exists(proteinQ_file_path)) {
    warning(paste("File not found:", proteinQ_file_path))
    return(NULL)
  } else {
    proteinQ_tblastn <- read_blast(file = proteinQ_file_path)
    proteinQ_tblastn$ID <- 1:nrow(proteinQ_tblastn)
  }

  # Check if domainQ_tblastn file exists
  domainQ_file_path <- paste0("results/tblastn/", GC, "/", domainQ, "/tblastn_results.txt")
  if (!file.exists(domainQ_file_path)) {
    warning(paste("File not found:", domainQ_file_path))
    return(NULL)
  } else {
    domainQ_tblastn <- read_blast(file = domainQ_file_path)
  }
  # Calculate minimum distance to domainQ_tblastn for each proteinQ_tblastn result
  dists <- sapply(proteinQ_tblastn$ID, function(i) {
    sseq<-proteinQ_tblastn$sseqid[i]
    matching_domains <- domainQ_tblastn[sseqid == sseq]

    if (nrow(matching_domains) == 0) {
      return(Inf)
    } else {
      return(min(abs(proteinQ_tblastn$sstart[proteinQ_tblastn$ID == i] - matching_domains$sstart)))
    }
  })

  # Add distance column to proteinQ_tblastn
  proteinQ_tblastn$dist_to_domain <- dists

  # Add GC column as the first column
  proteinQ_tblastn <- data.table(GC = GC, proteinQ_tblastn)
  proteinQ_tblastn<-proteinQ_tblastn[qstart>200]

  return(proteinQ_tblastn)
}

#MSH6 evolution project
