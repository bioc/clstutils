source('unit/setup.R')
VERBOSE <- TRUE

library(clstutils)
library(ape)

test_seqs <- function(){

  data(seqs)
  checkTrue(is.matrix(seqs))

}

test_seqdata <- function(){

  data(seqs)
  data(seqdat)
  checkTrue(is.data.frame(seqdat))
  checkTrue(nrow(seqs) == nrow(seqdat))  
  checkTrue(all(rownames(seqs) == seqdat$seqname))
}


