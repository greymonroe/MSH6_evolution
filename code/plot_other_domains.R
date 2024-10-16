#' Plot Other Domains in a Phylogenetic Tree
#'
#' This function visualizes domains in a given phylogenetic tree based on annotations from InterProScan.
#' It returns a list containing plots of the domain and the tree, as well as the organisms associated with the domain.
#'
#' @param interpro A data frame or data.table read from InterProScan annotations of proteins.
#' @param tree A phylogenetic tree object, typically representing the NCBI taxonomy.
#' @param proteins A data.table containing Organism and Protein columns.
#' @param domain A character string specifying the domain name, e.g., "zinc-finger of acetyl-transferase ESCO".
#' @param interprosource A character string specifying the analysis source from InterProScan, e.g., "Pfam".
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{domain_plot}: A ggplot object visualizing the domain.
#'   \item \code{tree_plot}: A ggplot object visualizing the phylogenetic tree.
#'   \item \code{Organisms}: A character vector of organism names associated with the domain.
#' }
#'
#' @examples
#' # Assuming you have interpro annotations 'my_interpro', a tree object 'my_tree', and MSH6 annotations 'my_msh6':
#' # plot_other_domains(my_interpro, my_tree, my_msh6, "zinc-finger of acetyl-transferase ESCO", "Pfam")
#'
#' @importFrom ggplot2 ggplot aes geom_segment
#' @importFrom ggtree ggtree geom_nodelab geom_tippoint geom_tiplab
#'
plot_other_domains <- function(interprosource, domain, interpro, tree, proteins) {
  domains<-interpro[Analysis==interprosource &  Signature_Description==domain]
  domain_proteins<-interpro[Analysis==interprosource & Protein %in% domains$Protein]
  domain_proteins$Organism<-proteins$Organism[match(domain_proteins$Protein, proteins$Protein)]
  domain_plot<-ggplot(domain_proteins, aes(x=Start_Location, xend=Stop_Location, y=Protein, yend=Protein, col=Signature_Description))+
    geom_segment()
  Organisms<-unique(domain_proteins$Organism)
  taxonomy<-lapply(Organisms, function(o) get_ancestral_nodes(tree,o))
  taxonomy <- taxonomy[sapply(taxonomy, function(x) !is.null(x))]

  common<-common_taxonomy(Organisms, tree, 1)

  subtree<-extract.clade(tree, common)

  subtree<-drop.tip(subtree, subtree$tip.label[!subtree$tip.label %in% proteins$Organism])
  proteins$tip<-(proteins$Organism %in% Organisms)
  proteins<-proteins[Organism %in% subtree$tip.label]

  tree_plot<-ggtree(subtree, size=0.1, layout="circular")
  tree_plot<-tree_plot %<+% as.data.frame(proteins) +geom_nodelab(size=1)+ geom_tippoint(aes(col=tip), size=0.01)+geom_tiplab(aes(col=tip), size = 0.5)

  return(list(domain_plot, tree_plot, Organisms))
}

out<-plot_other_domains("Pfam", "zinc-finger of acetyl-transferase ESCO", interpro, tree, msh6)
#MSH6 Evolution

