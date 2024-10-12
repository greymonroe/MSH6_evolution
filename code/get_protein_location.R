#' Extract Protein Location from GFF File
#'
#' This function reads a GFF file and extracts the location information of the protein.
#' It returns the chromosome (CHROM), start (START), and stop (STOP) positions,
#' as well as the minimum and maximum of these positions. The INFO column contains the file name.
#'
#' @param file Character. The path to the GFF file to be read.
#'
#' @return
#' A data.table containing the columns CHROM, START, STOP, and INFO.
#' If the GFF file is empty, the function returns NULL.
#'
#' @examples
#' # Assuming you have a GFF file named "example.gff":
#' protein_location <- get_protein_location("example.gff")
#'
get_protein_location<-function(file){
  gff<-fread(file, sep="\t")
  if(nrow(gff)==0) return(NULL)
  colnames(gff) <- c("CHROM","SOURCE","TYPE","START","STOP","SCORE","DIRECTION","PHASE","INFO")
  location=data.table(CHROM=unique(gff$CHROM), START=min(gff$START), STOP=max(gff$STOP), INFO=file)
  return(location)
}
