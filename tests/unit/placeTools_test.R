source('unit/setup.R')
VERBOSE <- TRUE

library(clstutils)
library(ape)

test_placefile <- function(){
  checkTrue(file.exists(jsonfile))
}

test_distfile <- function(){
  checkTrue(file.exists(distfile))
}

test_placeMat01 <- function(){
  
  ## distfile contains a triangular matrix
  ## read uncompressed file

  pmat <- clstutils:::placeMat(distfile)
  lines <- readLines(file(distfile))
  checkTrue(nrow(pmat$ddist) == length(lines) + 1)
  checkTrue(nrow(pmat$paths) == length(lines) + 1)
}

test_placeMat02 <- function(){
  checkTrue(file.exists(distfilez))
  
  ## distfile contains a triangular matrix
  ## read compressed file
  
  pmat <- clstutils:::placeMat(distfilez)
  lines <- readLines(file(distfilez))
  checkTrue(nrow(pmat$ddist) == length(lines) + 1)
  checkTrue(nrow(pmat$paths) == length(lines) + 1)
}

test_treeDists01 <- function(){
  x <- treeDists(distfile=distfile, placefile=jsonfile)
  checkTrue(setequal(c('dists','paths','dmat'), names(x)))
  checkTrue(sum(diag(x$dmat)) == 0)  
}

## test_placeData01 <- function(){
##   dat <- placeData(placefile)
##   checkTrue(all(colnames(dat) == c('name','hit','at','mlwr','ppost','mlll','bml','edge','branch','tax_id')))
## }

test_edgeMap01 <- function(){

  emap <- clstutils:::edgeMap(jsonfile)
  ## confirm that the number of sequences is the same as the refpkg
  refpkg <- refpkgContents(vag_refpkg)
  fafile <- refpkg$files$aln_fasta
  library(ape)
  seqs <- read.dna(fafile, format='fasta')
  checkTrue(nrow(emap) == nrow(seqs))
  
  ## sequence names all the same?
  checkTrue(setequal(emap$name, rownames(seqs)))
  
}

test_classifyPlacements01 <- function(){
  taxdata <- taxonomyFromRefpkg(vag_refpkg)
  treedists <- treeDists(distfile=distfile, placefile=jsonfile)

  save(taxdata, treedists, file='unit/test_classifyPlacements01.rda')

  
  placetab <- data.frame(name=c('lcrisp'),
                         at=c(860),
                         edge=c(0.000010),
                         branch=c(0.030468))

  cc <- classifyPlacements(taxdata, treedists, placetab, debug=FALSE)##, dStart=0.89)
  print(cc)

  ## this doesn't really test anything meaningful - TODO: come up with
  ## something better
  ## checkTrue(all(cc$at == c(358,342)))
}
