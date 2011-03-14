classifyPlacements <- function(taxdata, treedists, placetab, ..., verbose=FALSE, debug=FALSE){

  ## Perform classification of placements by pplacer on a phylogenetic
  ## tree using clst.
  ##
  ## Input
  ## -----
  ##
  
  vcat <- function(msg){if(verbose){cat(msg)}}
  vprint <- function(msg){if(verbose){print(msg)}}
  
  ## dmat is the square matrix of distances between tips in the
  ## reference tree
  dmat <- treedists$dmat
  
  ## taxonomic data
  taxNames <- taxdata$taxNames
  taxTab <- taxdata$taxTab
  ## TODO: define leftmost column as root in clst::classifyIter
  taxTab$tax_id <- NULL
  
  ## distances: dvects is a rectangular matrix containing distances
  ## between nodes (rows) and tips (columns)
  requiredCols <- c('at','edge','branch')
  if(!all(requiredCols %in% colnames(placetab))){
    stop(cat(paste('placetab must have all of the columns', requiredCols, '\n')))
  }

  dvects <- with(placetab, {
    treedists$dists[at+1,,drop=FALSE] + treedists$paths[at+1,,drop=FALSE]*edge + branch
  })
  
  ## compatibility checks for input data:
  if(!nrow(taxTab) == nrow(dmat) | !nrow(taxTab) == ncol(dmat)){
    stop('treedists$dmat must have an edge length equal to the number of rows in taxdata$taxTab')
  }
  
  stopifnot(all(rownames(taxdata$taxTab) == rownames(treedists$dmat)))
  
  ## perform the classification.  each list element contains results
  ## from a single atgroup; perform the classification and summarize the
  ## results

  ## classified <- lapply(1:10,
  classified <- lapply(seq(nrow(dvects)),
                       function(i){
                         
                         vcat(gettextf('*** i = %s ***\n',i))

                         cc <- clst::classifyIter(dmat, taxTab, dvects[i,], ..., verbose=debug)

                         last <- cc[[length(cc)]]
                         d <- with(last, {
                           if(length(matches) > 0){
                             tids <- matches
                           }else{
                             ## if no match, save the tally for the taxon
                             ## with the lowest median distance to the query
                             tids <- rownames(tally)[which.min(tally$median)]
                           }
                           
                           cbind(data.frame(tax_id=tids,
                                            tax_name=taxNames[tids],
                                            rank=names(rank)),
                                 tally[tids,,drop=FALSE])

                         })

                         d$d <- last$thresh$D
                         d$at <- placetab$at[i]
                         d$atgroup <- placetab$atgroup[i]

                         ## d$tally <- tab$tally[i]
                         vprint(d)
                         return(d)
                       })

  cdata <- do.call(rbind, classified)
  rownames(cdata) <- NULL  
    
  return(cdata)  
}


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

placeMat <- function(file){

  ## Parses a lower triangular matrix of between-edge distances
  ## produced by placeutil --distmat
  ##
  ## Input
  ## -----
  ##
  ## * file - name of input file
  ##
  ## Output
  ## ------
  ##
  ## A list(ddist,paths) of two square matrices
  ##
  ## * ddist - distances between distal ends of pairs of edges.
  ## * paths - Describes the path between edge i,j as parallel (-1) or
  ##   serial (1).

  ## some notes from Erick:
  ##   a placement on an edge looks like this:

  ## proximal
  ##  |
  ##  |   d_p
  ##  |
  ##  |---- x
  ##  |
  ##  |   d_d
  ##  |
  ##  |
  ## distal

  ## d_p is the distance from the placement x to the proximal side of the
  ## edge, and d_d the distance to the distal side.

  ## If the distance from x to a leaf y is an S-distance Q, then the
  ## path from x to y will go through the distal side of the edge and
  ## we will need to add d_d to Q to get the distance from x to y.  If
  ## the distance from x to a leaf y is a P-distance Q, then the path
  ## from x to y will go through the proximal side of the edge, and we
  ## will need to subtract off d_d from Q to get the distance from x
  ## to y. In either case, we always need to add the
  ## length of the pendant edge, which is the second column.

  ## To review, say the values of the two leftmost columns are a and b
  ## for a given placement x, and that it is on an edge i.  We are
  ## interested in the distance of x to a leaf y, which is on edge j.
  ## We look at the distance matrix, entry (i,j), and say it is an
  ## S-distance Q. Then our distance is Q+a+b.  If it is a P-distance
  ## Q, then the distance is Q-a+b.

  ## The distances between leaves should always be P-distances, and there
  ## we need no trickery.

  v <- scan(file, what="character", quiet=TRUE)
  
  N <- ceiling(sqrt(length(v)*2))

  stopifnot(N*(N-1)/2 == length(v))

  m <- matrix(nrow=N, ncol=N)
  m[upper.tri(m)] <- v

  ## numeric values
  ddist <- substring(m,2)
  mode(ddist) <- "numeric"

  lt <- lower.tri(ddist)
  ddist[lt] <- t(ddist)[lt]  

  diag(ddist) <- 0

  ## textual values (P,S)
  pstrs <- substring(m,1,1)
  
  ## these values will be applied to the position of placement on an
  ## edge; P (parallel) leaves are corrected by multiplying by -1, S
  ## (serial) by 1
  paths <- ifelse(pstrs == 'P', -1L, 1L)

  paths[lt] <- t(paths)[lt]
  
  diag(paths) <- 1L
  rownames(paths) <- NULL

  list(ddist=ddist, paths=paths)
}

## placeMat <- function(file){

##   ## Parses a lower triangular matrix of between-edge distances
##   ## produced by placeutil --distmat
##   ##
##   ## Input
##   ## -----
##   ##
##   ## * file - name of input file
##   ##
##   ## Output
##   ## ------
##   ##
##   ## A list(ddist,paths) of two square matrices
##   ##
##   ## * ddist - distances between distal ends of pairs of edges.
##   ## * paths - Describes the path between edge i,j as parallel (-1) or
##   ##   serial (1).

##   ## some notes from Erick:
##   ##   a placement on an edge looks like this:

##   ## proximal
##   ##  |
##   ##  |   d_p
##   ##  |
##   ##  |---- x
##   ##  |
##   ##  |   d_d
##   ##  |
##   ##  |
##   ## distal

##   ## d_p is the distance from the placement x to the proximal side of the
##   ## edge, and d_d the distance to the distal side.

##   ## If the distance from x to a leaf y is an S-distance Q, then the
##   ## path from x to y will go through the distal side of the edge and
##   ## we will need to add d_d to Q to get the distance from x to y.  If
##   ## the distance from x to a leaf y is a P-distance Q, then the path
##   ## from x to y will go through the proximal side of the edge, and we
##   ## will need to subtract off d_d from Q to get the distance from x
##   ## to y. In either case, we always need to add the
##   ## length of the pendant edge, which is the second column.

##   ## To review, say the values of the two leftmost columns are a and b
##   ## for a given placement x, and that it is on an edge i.  We are
##   ## interested in the distance of x to a leaf y, which is on edge j.
##   ## We look at the distance matrix, entry (i,j), and say it is an
##   ## S-distance Q. Then our distance is Q+a+b.  If it is a P-distance
##   ## Q, then the distance is Q-a+b.

##   ## The distances between leaves should always be P-distances, and there
##   ## we need no trickery.

##   ## N <- max(count.fields(file))+1
##   N <- length(readLines(file)) + 1

##   x <- scan(file, what=as.list(rep('character',N)), fill=TRUE, quiet=TRUE)
##   m <- do.call(cbind, x)
##   m <- rbind(rep("",N), m)

##   ## numeric values
##   ddist <- substring(m,2)
##   mode(ddist) <- "numeric"
##   ut <- upper.tri(ddist)
##   ddist[ut] <- t(ddist)[ut]
##   diag(ddist) <- 0

##   ## textual values (P,S)
##   pstrs <- substring(m,1,1)
  
##   ## these values will be applied to the position of placement on an
##   ## edge; P (parallel) leaves are corrected by multiplying by -1, S
##   ## (serial) by 1
##   paths <- ifelse(pstrs == 'P', -1L, 1L)

##   paths[ut] <- t(paths)[ut]
##   diag(paths) <- 1L
##   rownames(paths) <- NULL

##   list(ddist=ddist, paths=paths)
## }

treeDists <- function(distfile, placefile, pmat, pmap){

  ## Provides objects (dists, paths) that can be used to calculate
  ## vectors of distances between an internal node and each leaf
  ## node. Also returns a square matrix of distances between leaf
  ## nodes.
  ##
  ## Input
  ## -----
  ##
  ## * distfile - path to output of 'placeutil --distmat'
  ## * placefile - path to pplacer output
  ## * pmat - square matrix (output of clstutils:::placeMat); may be
  ##   provided in the place of distfile
  ## * pmap - output of clstutils:::placeMap (data.frame); may be
  ##   provided in the place of pmap

  ##
  ## Value
  ## -----
  ##
  ## A list with named elements $dists, $paths, and $dmat
  ##
  ##  * $dists, $paths - rectangular matrices with rows corresponding
  ##    to all nodes in pplacer order, and columns corresponding to
  ##    tips in phylo order.
  ##  * $dmat - square matrix containing distances between pairs of tips.
  ##

  if(missing(pmat)){
    if(missing(distfile)){
      stop('"distfile" is required if "pmat" is missing.')
    }
    pmat <- placeMat(distfile)
  }

  if(missing(pmap)){
    if(missing(placefile)){
      stop('"placefile" is required if "pmap" is missing.')
    }
    tre <- placeTree(placefile)
    pmap <- placeMap(tre)
  }

  for(colname in c('type','parent')){
    if(!colname %in% colnames(pmap)){
      stop(paste('pmap must contain a column named',colname))
    }
  }
  
  list(
       dists = .rectangle(pmat$ddist, pmap),
       paths = .rectangle(pmat$paths, pmap),
       dmat = with(subset(pmap, type == 'tip'),{
         m <- pmat$ddist[pplacer+1, pplacer+1]
         colnames(m) <- rownames(m) <- name
         m })
       )
}

.rectangle <- function(m, pmap){

  ## Rehsape a square matrix for the purpose of calculating distances
  ## between a given node and each tip.
  ##
  ## Input
  ## -----
  ##
  ## * m - square matrix
  ## * pmap - output of placeMap (data.frame)
  ##
  ## Value
  ## -----
  ##
  ## A rectangular matrix with rows corresponding to all nodes in
  ## pplacer order, and columns corresponding to tips in phylo order.

  ## columns: tips in phylo order
  ## tips is in the same order as tre$tip.label:
  ## ie, all(unplace(tre)$tip.label == tips$name)
  tips <- subset(pmap, pmap$type == 'tip')
  out <- m[,tips$pplacer+1] ## convert pplacer numbering to 1-index
  colnames(out) <- tips$name

  ## rows: all nodes in pplacer order
  rownames(out) <- with(subset(pmap, !is.na(parent)), label[order(pplacer)])

  return(out)
}

placeData <- function(file){

  ## Input
  ## -----
  ##
  ## * file - output of pplacer (.place file); may be the name of a
  ##   file or an open connection.
  ##
  ## Output
  ## ------
  ##
  ## * a data.frame with columns
  ##   'name','hit','at','mlwr','ppost','mlll','bml','edge','branch'
  ##   (see Details)
  ##
  ## Details
  ## -------
  ##
  ## Columns in the data.frame output include 'name' (name of the
  ## sequence), 'hit' (integer corresponding to each row of output for
  ## a given sequence in the pplacer file); plus the following columns
  ## as described in
  ## (http://matsen.fhcrc.org/pplacer/manual.html#placefile):
  ##
  ## 1. Placement edge number, numbered as in a postorder traversal
  ## 2. ML likelihood weight ratio (i.e. the normalized ML likelihood values)
  ## 3. Posterior probability (or a dash if the -p option was not set)
  ## 4. ML log likelihood
  ## 5. Bayes marginal likelihood (or a dash if the -p option was not set)
  ## 6. The ML distance from the distal (farthest from the root) side of the edge
  ## 7. The ML pendant branch length

  cnames <- c('at','mlwr','ppost','mlll','bml','edge','branch','tax_id','x')

  lines <- readLines(file)
  commentIdx <- grep('^#', lines, value=FALSE, perl=TRUE)
  nameIdx <- grep('^>', lines, value=FALSE, perl=TRUE)

  names <- c('none', substring(lines[nameIdx], 2))

  ## remove lines with comments, names and sequences
  discard <- c(commentIdx, nameIdx, nameIdx+1)
  lines[discard] <- paste(rep('\t',length(cnames)-1),collapse='')

  con <- textConnection(lines)
  tab <- read.table(con, sep='\t')
  close(con)
  colnames(tab) <- cnames

  ## rows corresponding to each name extend from here....
  first <- c(1, nameIdx)
  ## ... to here
  last <- c(nameIdx - 1, length(lines))

  each <- last - first + 1

  ## add sequence names; no idea why tab$name <- rep(names, each=last -
  ## first + 1) doesn't give the expected behavior here
  tab$name <- unlist(sapply(seq_along(names), function(i) rep(names[i],each[i])))
  ## subtract 2 from each value of each: first two lines are name, sequence
  tab$hit <- as.integer(unlist(sapply(seq_along(each), function(i) seq(each[i]))) - 2)

  ## delete unneeded rows and rearrange columns
  tab <- tab[-discard, c('name','hit', setdiff(cnames,'x'))]
  tab$tax_id <- substring(tab$tax_id, 3)

  return(tab)
}


