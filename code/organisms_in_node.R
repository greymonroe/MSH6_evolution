organisms_in_node <- function(tree, Organisms, node, exclude = FALSE) {

  # Initialize a vector to store tip labels
  all_tips <- c()

  # Loop through each node to gather or exclude tip labels
  for (n in node) {
    subtree <- ape::extract.clade(tree, node = n)
    if (exclude) {
      tree <- ape::drop.tip(tree, subtree$tip.label)
    } else {
      all_tips <- c(all_tips, subtree$tip.label)
    }
  }

  # If exclude mode, get the remaining tip labels from the tree
  if (exclude) {
    all_tips <- tree$tip.label
  }

  # Filter the Organisms data.table based on the gathered tip labels
  sub_organisms <- Organisms[Organism %in% all_tips]

  return(sub_organisms)
}

#MSH6 evolution project
