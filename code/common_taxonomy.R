
#' Find the Common Node in a Phylogenetic Tree
#'
#' This function identifies the common node among a set of organisms in a given phylogenetic tree.
#' It then returns the node that is a specified depth below the common node.
#'
#' @param Organisms A character vector of organism names from NCBI.
#' @param tree A phylogenetic tree object, typically read using the `read.tree` function from the `ape` package.
#' @param depth An integer specifying how many levels below the common node to go.
#'
#' @return The node that is a specified depth below the common node. If no common node is found, the function returns `NULL`.
#'
#' @examples
#' # Assuming you have a tree object named 'my_tree' and a vector of organisms named 'my_organisms':
#' # common_node(my_organisms, my_tree, 2)
#'
#' @importFrom stats setNames
#'

common_taxonomy <- function(Organisms, tree, depth) {

  taxonomy<-lapply(Organisms, function(o) get_ancestral_nodes(tree,o))
  taxonomy <- taxonomy[sapply(taxonomy, function(x) !is.null(x))]

  common_elements <- taxonomy[[1]]

  for (i in 2:length(taxonomy)) {
    common_elements <- intersect(common_elements, taxonomy[[i]])
    if (length(common_elements) == 0) {
      return(NULL)
    }
  }

  return(common_elements[depth])
}

#MSH6 evolution project
