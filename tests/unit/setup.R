unit_output <- 'unit_output'

## create destination for output files
dir.create(unit_output, showWarnings=FALSE)

VERBOSE <- TRUE

unitFuncName <- function(){
  unlist(sapply(sys.frames(), function(f)f$funcName))[1]
}

pdfName <- function(){
  gettextf('%s.pdf',file.path(unit_output,unitFuncName()))
}

expand <- function(fname, destdir, package='clstutils'){

  archive <- system.file('extdata', fname, package=package)  

  if(grepl('.tar', fname, fixed=TRUE)){
    unarch <- sub('\\.tar\\..+$', '', fname)
    system(gettextf('tar -xf %s --directory=%s', archive, destdir))
  }else if(grepl('.bz2', fname, fixed=TRUE)){
    unarch <- sub('\\.bz2$', '', fname)
    system(gettextf('cp %s %s', archive, destdir))
    system(gettextf('bzip2 -df %s', file.path(destdir, fname)))    
  }else{
    stop(paste('problem expanding',fname))
  }

  file.path(destdir, unarch)
}

vag_refpkg <- expand('vaginal_16s.refpkg.tar.bz2', unit_output)
distfile <- expand('merged.distmat.bz2', unit_output)

jsonfile <- system.file('extdata', 'merged.json', package='clstutils')
distfilez <- system.file('extdata', 'merged.distmat.bz2', package='clstutils')

