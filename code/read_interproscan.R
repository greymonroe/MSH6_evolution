#' Read InterProScan Results
#'
#' This function reads the results from an InterProScan analysis.
#' InterProScan provides comprehensive protein domain and functional site predictions.
#' The function reads the file and assigns appropriate column names to the data.
#'
#' @param file Character. The path to the InterProScan results file to be read.
#'
#' @return
#' A data.table containing the results from the InterProScan analysis with appropriate column names.
#'
#' @examples
#' # Assuming you have an InterProScan results file named "interpro_results.txt":
#' interpro_data <- read_interproscan("interpro_results.txt")
#'
#'
read_interproscan<-function(file){
interpro<-fread(file, header = F, quote="", sep="\t", fill=T)
interproscan_columns <- c(
  "Protein",
  "Sequence_MD5_digest",
  "Sequence_Length",
  "Analysis",
  "Signature_Accession",
  "Signature_Description",
  "Start_Location",
  "Stop_Location",
  "Score",
  "Status",
  "Date",
  "InterPro_Accession",
  "InterPro_Description",
  "GO_Annotations"
)
colnames(interpro)<-interproscan_columns
return(interpro)
}

#MSH6 evolution project

