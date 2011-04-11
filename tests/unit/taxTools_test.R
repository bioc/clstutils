source('unit/setup.R')
VERBOSE <- TRUE

library(ape)
library(clstutils)

data(seqs)
data(seqdat)



# [1] "Enterococcus avium"    "Enterococcus faecalis" "Enterococcus faecium"

test_abspath <- function(){

  files <- c(
             '.',
             '..',
             'runfile.sh',
             'unit'
           )

  paths <- clstutils:::.abspath(files)
  checkTrue(all(file.exists(paths)))  
}

test_refpkgContents <- function(){

  contents <- refpkgContents(vag_refpkg)
  checkTrue(all(c('files','md5','metadata') %in% names(contents)))
  checkTrue(setequal(names(contents$files),names(contents$md5)))  
}

test_taxonomyFromRefpkg01 <- function(){

  taxdata <- taxonomyFromRefpkg(vag_refpkg)

  ## all columns should be character vectors
  for(col in colnames(taxdata$taxTab)){
    checkTrue(typeof(taxdata$taxTab[[col]]) == 'character')    
  }
  
  with(taxdata, checkTrue(all(taxTab$tax_id %in% names(taxNames))))

  ## confirm that the order of taxTab reflects that of seq_info
  contents <- refpkgContents(vag_refpkg)
  seq_info <- read.csv(contents$files$seq_info, colClasses='character')

  checkTrue(all(rownames(taxdata$taxTab) == seq_info$seqname))
  
}


test_findOutliers01 <- function(){

  taxon <- "Enterococcus faecium"
    
  someseqs <- seqs[seqdat$tax_name == taxon,]
  somedat <- subset(seqdat, seqdat$tax_name == taxon)
  checkTrue(all(rownames(someseqs) == somedat$seqname))
  
  dmat <- dist.dna(someseqs, pairwise.deletion=TRUE, as.matrix=TRUE, model='raw')

  outliers <- findOutliers(dmat, cutoff=0.015)
  
}

test_findOutliers02 <- function(){

  taxon <- "Enterococcus faecium"
    
  someseqs <- seqs[seqdat$tax_name == taxon,]
  somedat <- subset(seqdat, seqdat$tax_name == taxon)
  checkTrue(all(rownames(someseqs) == somedat$seqname))
  
  dmat <- dist.dna(someseqs, pairwise.deletion=TRUE, as.matrix=TRUE, model='raw')

  ## values of quant closer to 0 are more stringent
  outliers1 <- findOutliers(dmat, quant=0)
  outliers2 <- findOutliers(dmat, quant=.5)
  outliers3 <- findOutliers(dmat, quant=1)

  checkTrue(sum(outliers1) > sum(outliers2))
  checkTrue(sum(outliers2) > sum(outliers3))   
}


test_maxDists01 <- function(){

  taxon <- "Enterococcus faecium"
    
  someseqs <- seqs[seqdat$tax_name == taxon,]
  somedat <- subset(seqdat, seqdat$tax_name == taxon)
  checkTrue(all(rownames(someseqs) == somedat$seqname))
  
  dmat <- dist.dna(someseqs, pairwise.deletion=TRUE, as.matrix=TRUE, model='raw')
  selected <- maxDists(dmat, N=1)
  checkTrue(length(selected) == 1)
}

test_maxDists02 <- function(){

  taxon <- "Enterococcus faecium"
    
  someseqs <- seqs[seqdat$tax_name == taxon,]
  somedat <- subset(seqdat, seqdat$tax_name == taxon)
  checkTrue(all(rownames(someseqs) == somedat$seqname))
  
  dmat <- dist.dna(someseqs, pairwise.deletion=TRUE, as.matrix=TRUE, model='raw')
  selected <- maxDists(dmat, N=10)
  checkTrue(length(selected) == 10)
}

test_maxDists03 <- function(){

  taxon <- "Enterococcus faecium"
    
  someseqs <- seqs[seqdat$tax_name == taxon,]
  somedat <- subset(seqdat, seqdat$tax_name == taxon)
  checkTrue(all(rownames(someseqs) == somedat$seqname))
  
  dmat <- dist.dna(someseqs, pairwise.deletion=TRUE, as.matrix=TRUE, model='raw')
  outliers <- findOutliers(dmat, cutoff=0.015)
  selected <- maxDists(dmat, N=10, exclude=outliers)

  checkTrue(length(selected) == 10)
  checkTrue(all(!selected %in% seq(nrow(dmat))[outliers]))
}

