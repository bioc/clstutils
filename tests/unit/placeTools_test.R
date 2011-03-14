source('unit/setup.R')
VERBOSE <- TRUE

library(clstutils)
library(ape)

placeMap <- clstutils:::placeMap
placeMat <- clstutils:::placeMat
placeTree <- clstutils:::placeTree

test_placefile <- function(){
  checkTrue(file.exists(placefile))
}

test_distfile <- function(){
  checkTrue(file.exists(distfile))
}

test_placeTree01 <- function(){
  tre <- placeTree(placefile)
  checkTrue(class(tre) == 'phylo')
}

test_placeTree02 <- function(){
  ## works on compressed file?
  checkTrue(file.exists(placefilez))
  tre <- placeTree(placefilez)
  checkTrue(class(tre) == 'phylo')
}

test_placeMat01 <- function(){

  ## distfile contains a triangular matrix
  pmat <- placeMat(distfile)
  lines <- readLines(file(distfile))
  checkTrue(nrow(pmat$ddist) == length(lines) + 1)
  checkTrue(nrow(pmat$paths) == length(lines) + 1)
}

test_placeMat02 <- function(){
  checkTrue(file.exists(distfilez))
  
  ## distfile contains a triangular matrix
  pmat <- placeMat(distfilez)
  lines <- readLines(file(distfilez))
  checkTrue(nrow(pmat$ddist) == length(lines) + 1)
  checkTrue(nrow(pmat$paths) == length(lines) + 1)
}

test_placeMap01 <- function(){

  ptree <- placeTree(placefile)
  pmap <- placeMap(ptree)

  checkTrue(with(pmap, all(grepl('@', subset(label, type=='tip')))))
  checkTrue(with(pmap, !all(grepl('@', subset(label, type=='node')))))  
}

test_treeDists01 <- function(){
  checkException(treeDists(placefile=placefile))
}

test_treeDists02 <- function(){
  checkException(treeDists(distfile=distfile))
}

test_treeDists03 <- function(){
  x <- treeDists(distfile=distfile, placefile=placefile)
  checkTrue(setequal(c('dists','paths','dmat'), names(x)))
}

test_placeData01 <- function(){
  dat <- placeData(placefile)
  checkTrue(all(colnames(dat) == c('name','hit','at','mlwr','ppost','mlll','bml','edge','branch','tax_id')))
}
