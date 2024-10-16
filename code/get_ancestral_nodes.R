

get_ancestral_nodes <- function(tree, tip_label) {
  if (tip_label %in% tree$tip.label) {
    tip_num <- which(tree$tip.label == tip_label)
    node <- tree$edge[which(tree$edge[, 2] == tip_num), 1]
    nodes <- numeric(0)

    while (!is.na(node) && node > length(tree$tip.label)) {
      nodes <- c(nodes, node)
      node <- ifelse(length(which(tree$edge[, 2] == node)) > 0,
                     tree$edge[which(tree$edge[, 2] == node), 1],
                     NA)
    }

    # Check if node labels exist and retrieve them
    if (!is.null(tree$node.label)) {
      node_names <- tree$node.label[nodes - length(tree$tip.label)]
      return(node_names)
    } else {
      return(nodes)
    }

  } else {
    warning("The provided tip label does not exist in the tree.")
    return(NULL)
  }
}

#MSH6 evolution project
