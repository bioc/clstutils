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

  stopifnot(setequal(rownames(taxTab), rownames(dmat)))

  ## taxTab and margins of dmat must be in the same order
  taxTab <- taxTab[match(rownames(dmat),rownames(taxTab)),]

  stopifnot(all(rownames(taxTab) == rownames(dmat)))

  
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

edgeMap <- function(placefile){

  treestr <- fromJSON(file=placefile)$tree

  spl <- strsplit(treestr, '[,(]')[[1]]

  edgelist <- lapply(spl, function(s){
    data.frame(
               name = strsplit(s, ':')[[1]][1],
               edge = as.numeric(strsplit(s, '(\\]|\\[)')[[1]][2])
               )
  })[grepl(':[^[]+\\[',spl)]
  
  edgemap <- do.call(rbind, edgelist)

  ## i corresponds to margins of placeMat output
  edgemap$i <- edgemap$edge + 1

  return(edgemap)
}


treeDists <- function(placefile, distfile){

  ## Provides objects (dists, paths) that can be used to calculate
  ## vectors of distances between an internal node and each leaf
  ## node. Also returns a square matrix of distances between leaf
  ## nodes.
  ##
  ## Input
  ## -----
  ##
  ## * placefile - path to pplacer output
  ## * distfile - path to output of 'placeutil --distmat'

  ##
  ## Value
  ## -----
  ##
  ## A list with named elements $dists, $paths, and $dmat
  ##
  ##  * $dists, $paths - rectangular matrices with rows corresponding
  ##    to all nodes, and named columns corresponding to tips.
  ##  * $dmat - square matrix containing distances between pairs of tips.
  ##


  edgemap <- edgeMap(placefile)

  ## TODO: if distfile is missing, create using `guppy distmat`

  pmat <- placeMat(distfile)

  leaves <- edgemap$i
  
  dists <- pmat$ddist[,leaves]
  paths <- pmat$paths[,leaves]
  colnames(dists) <- colnames(paths) <- edgemap$name

  dmat <- pmat$ddist[leaves,leaves]
  colnames(dmat) <- rownames(dmat) <- edgemap$name
  
  list(dists=dists, paths=paths, dmat=dmat)
       
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


