
#' Create a Table for MSH6 Organisms
#'
#' This function processes a file containing predicted MSH6 annotations and
#' returns a table summarizing the presence or absence of MSH6 and its domains
#' in different organisms. The function also considers the distance between
#' MSH6 orthologs and possible misannotated domains in neighboring regions.
#'
#' @param file Character. The path to the file containing predicted MSH6 annotations.
#'             Default is "predicted_msh6_blast_annotated.txt".
#' @param dist_threshold Numeric. The distance threshold to consider for
#'                        distinguishing between present, absent, and ambiguous domains.
#'
#' @return A data.table containing the summarized information for each organism.
#'         Columns include the organism name, presence/absence of MSH6 and its domains,
#'         and the minimum distance to the nearest domain.
#'
#' @examples
#' # Assuming you have the necessary data set up:
#' result <- make_MSH6_organism_table("path_to_file.txt", 1000)
#'
make_MSH6_organism_table <- function(file, dist_threshold) {
  msh6<-fread(file)

  msh6_organism <-
    msh6[, .(
      MutS = sum(MutS),
      Tudor_domain = sum(Tudor_domain + Tudor_interpro),
      Tudor_dist = min(tudor_tblastn, tudor_neighbor),
      PWWP_domain = sum(PWWP_domain+ PWWP_interpro),
      PWWP_dist = min(PWWP_tblastn, PWWP_neighbor)
    ), by = Organism]

  msh6_organism$Tudor_call<-apply(msh6_organism, 1, function(row){
    domain=as.numeric(row["Tudor_domain"])
    dist=as.numeric(row["Tudor_dist"])
    if(domain>0) return("Present")
    if(dist<dist_threshold) return("Ambiguous")
    return("Absent")
  })

  msh6_organism$PWWP_call<-apply(msh6_organism, 1, function(row){
    domain=as.numeric(row["PWWP_domain"])
    dist=as.numeric(row["PWWP_dist"])
    if(domain>0) return("Present")
    if(dist<dist_threshold) return("Ambiguous")
    return("Absent")
  })

  ambiguous_msh6<-nrow(msh6_organism[MutS<1])
  ambiguous_Tudor<-nrow(msh6_organism[Tudor_call=="Ambiguous"])
  ambiguous_PWWP<-nrow(msh6_organism[PWWP_call=="Ambiguous"])

  message("removing ", ambiguous_msh6, " organisms with ambiguous MSH6 ortholog (lacking interpro support)")
  message("removing ", ambiguous_Tudor, " organisms with ambiguous Tudor domain calls (possible neighbor)")
  message("removing ", ambiguous_PWWP, " organisms with ambiguous PWWP domain calls (possible neighbor)")

  final_msh6<-msh6_organism[MutS>0 & Tudor_call!="Ambiguous" & PWWP_call!="Ambiguous"]
  message(nrow(final_msh6)," organisms with definitive MSH6 and domains")

  return(final_msh6)
}

#MSH6 evolution project


