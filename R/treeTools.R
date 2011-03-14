descendants <- function(x, tipsToRoot){
  ## Identify terminal nodes (leaves) associated with each internal
  ## node of a tree.
  ##
  ## Input
  ## =====
  ##
  ## * x - an object of class phylo, eg x <- nj(ddist)
  ## * tipsToRoot - output of call to (private) ape function seq_root2tip
  ##
  ## Value
  ## =====
  ##
  ## A list indexed by node, each element consisting of a numeric
  ## vector of dependent leaves

  ## x <- read.tree('/home/nhoffman/seqtax/bv/mltree/output/RAxML_result.bv_amplicon')

  ## a list defining sequence from the root to each edge.
  if(missing(tipsToRoot)){
    tipsToRoot <- .Call("seq_root2tip", x$edge, length(x$tip.label), x$Nnode, PACKAGE = "ape")
  }

  ## data structure to hold leaves for each node
  nodes <- lapply(seq(max(x$edge)), function(i) numeric(0))

  ## TODO: uses a nested loop: find a more efficient way,
  ## perhaps having something to do with subtrees
  for(lineage in tipsToRoot){
    nedge <- length(lineage)
    leaf <- lineage[nedge]
    for(node in lineage[-nedge]){
      nodes[[node]] <- c(nodes[[node]], leaf)
    }
  }

  ## reduce each set of leaves
  lapply(nodes, unique)

}

getTips <- function(x){
  ## Returns a vector of tip indices in tree order
  ## * x - an object of class phylo, eg x <- nj(ddist)

  Ntip <- length(x$tip.label)
  x$edge[x$edge[, 2] <= Ntip, 2]
}

getNodes <- function(x){
  ## Returns a vector of internal node indices in tree order
  ## * x - an object of class phylo, eg x <- nj(ddist)

  Ntip <- length(x$tip.label)
  x$edge[x$edge[, 2] > Ntip, 2]
}
