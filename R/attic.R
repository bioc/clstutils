placeTree <- function(file){

  ## Read numbered reference tree from a .place file;
  ## return an object of class phylo

  tag <- '# numbered reference tree: '

  line <- grep(tag,
               readLines(file, n=50), ## n=50 should be safe to capture header info
               value=TRUE)
  treestr <- substr(line, nchar(tag) + 1, nchar(line))

  ape::read.tree(text=treestr)
}

unplace <- function(ptree){
  ## Return a new phylo object with leading "N@" labels from
  ## stripped from tip.labels.

  ptree$tip.label <- sapply(strsplit(ptree$tip.label, split='@'), '[', 2)
  return(ptree)
}

.pplacerTreeExplained <- function(tre){

  nodes <- as.numeric(tre$node.label)
  tips <- as.numeric(sapply(strsplit(tre$tip.label, split='@'), '[[', 1))

}

placeMap <- function(ptree){

  ## Describes the mapping between node numbers as assigned by
  ## phylo{ape} and pplacer.
  ##
  ## Input
  ## -----
  ##
  ## * ptree - an object of class phylo corresponding to the "numbered
  ##   reference tree" line in a .place file
  ##
  ## Output
  ## ------
  ##
  ## A data.frame with the following columns:
  ##
  ## * node - corresponds to ptree$edge[2,] (numeric)
  ## * parent - corresponds to ptree$edge[1,] (numeric)
  ## * type - 'node' or 'tip' (factor)
  ## * pplacer - node or leaf number as defined in pplacer tree (numeric)
  ## * name - name of corresponding leaf in reference tree,
  ##   or pplacer edge number (character)
  ## * label - edge or leaf label in pplacer tree (character)
  ##
  ## Details
  ## -------
  ##
  ## Note that nrow(m) == length(ptree$node.label) + length(ptree$tip.label) - 1
  ## and that there will be no row corresponding to the root node, hence
  ## setdiff(seq(length(ptree$node.label) + length(ptree$tip.label)), pm$node) \
  ##   == length(ptree$tip.label) + 1
  ##
  ## Examples
  ## --------
  ##
  ## m <- placeMap(ptree)
  ## with(m, node[match(891,pplacer)]) ## maps pplacer numbers to ape numbers
  ## with(m, pplacer[match(450,node)]) ## maps ape numbers to pplacer numbers

  ## Notes:
  ## tip.label: a vector of mode character giving the names of the tips; the
  ##           order of the names in this vector corresponds to the
  ##           (positive) number in ‘edge’.

  ## edge: a two-column matrix of mode numeric where each row represents
  ##       an edge of the tree; the nodes and the tips are symbolized
  ##       with numbers; the tips are numbered 1, 2, ..., and the nodes
  ##       are numbered after the tips. For each row, the first column
  ##       gives the ancestor.

  ## http://lists.r-forge.r-project.org/pipermail/phylobase-devl/2008-August/000206.html
  ## For a tree bar with 5 taxa, node i has label bar$node.label[i - 5],
  ## the same is done with the tip.label, but you do not have to
  ## substract the number of taxa.

  ntip <- length(ptree$tip.label)

  spl <- strsplit(ptree$tip.label, split='@')
  leaves <- sapply(spl, '[', 1)
  leaf.names <- sapply(spl, '[', 2)

  nodes <- ptree$node.label

  ## need to add ancestry of the root node
  edges <- rbind(c(NA, length(ptree$tip.label)+1), ptree$edge)

  ## ii <- ptree$edge[,2]
  ii <- edges[,2]
  isNode <- ii > ntip
  nodeii <- ifelse(isNode, ii-ntip, NA)

  data.frame(
             node = ii,
             parent = edges[,1],
             type = factor(ifelse(isNode, 'node', 'tip')),
             pplacer = as.numeric(ifelse(isNode, nodes[nodeii], leaves[ii])),
             name = ifelse(isNode, nodes[nodeii], leaf.names[ii]),
             label = ifelse(isNode, nodes[nodeii], ptree$tip.label[ii]),
             stringsAsFactors=FALSE)
}

