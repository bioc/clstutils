prettyTree <- function(x, groups, fill, X, O, indices,
                       labels,
                       show=rep(TRUE, length(x)),
                       largs=list(cex=0.75, bty='n', yjust=0.5),
                       parargs=list(mar=c(bottom=5,
                                        left=2,
                                        top=2,
                                        right=ifelse(is.null(largs),2,8)),
                           xpd=NA),
                       pointargs=list(),
                       glyphs,
                       shuffleGlyphs=NA,
                       ...){

  ##
  ## Input
  ## -----
  ##
  ## * x - an object of class phylo, eg x <- nj(ddist)
  ## * groups - a factor (or object coercible) to a factor assigning group identity
  ##   to leaf nodes in x
  ## * fill - vector (logical or indices) of points to fill
  ## * X - vector of points to mark with an X
  ## * O - vector of points to mark with a circle
  ## * indices - label points with indices (all points if 'yes', or a subset indicated
  ##   by a vector)
  ## * labels - character vector of tip labels in the same order as x$tip.label
  ## * show - boolean vector of points to show
  ## * largs - arguments controlling appearance of the legend or NULL
  ##   for no legend
  ## * parargs - arguments to pass par()
  ## * pointargs - arguments to pass points() (other than pch, col, bg)
  ## * glyphs - a data.frame with columns named $col and $pch corresponding
  ##   to elements of unique(groups) - eg, output of getGlyphs
  ## * shuffleGlyphs - NA or an integer (argument to set.seed)  
  ## * ... passed to plot.tree
  ##
  ## Value
  ## -----
  ##
  ## Plots to the active device; no return value.
  ##
  ## Details
  ## -------
  ##
  ## Vectors specifying annotation should be in the order of row or column
  ## labels of the distance matrix used to generate x.
  
  groups <- factor(groups)

  if(missing(glyphs)){
    glyphs <- clst:::getGlyphs(length(unique(groups)),
                        shuffleGlyphs=shuffleGlyphs,
                        shapes=c(cir=21, sq=22, diam=23, tri=24, utri=25)) 
  }

  col <- as.character(glyphs$col[groups])
  pch <- glyphs$pch[groups]

  pch[!show] <- NA

  tips <- getTips(x)
  nTips <- length(tips)
  ##   par(mar=c(bottom=3, left=2, top=2, right=8),
  ##       xpd=NA)

  do.call(par, parargs)
  plot(x, show.tip.label = FALSE,
       ...)
  nje <- get("last_plot.phylo", .PlotPhyloEnv)
  ## axisPhylo()

  bg <- rep("transparent", nTips)
  if (!(missing(fill))) {
    bg[fill] <- col[fill]
  }

  xx <- nje$xx[tips]
  yy <- nje$yy[tips]

  ## points(x=xx,y=yy,pch=pch[tips],col=col[tips],
  ##        bg=bg[tips])

  do.call(points,
          c(list(x=xx,y=yy,pch=pch[tips],col=col[tips],bg=bg[tips]),
            pointargs)
          )

  if(!missing(indices)){
    iLabels <- rep(NA,nTips)
    iLabels[indices] <- seq_along(xx)[indices]
    text(xx,yy,labels=iLabels[tips],pos=4)
  }

  if(!missing(labels)){
    text(xx,yy,labels=labels,pos=4)
  }

  if(!missing(X)){
    if(is.numeric(X)) X <- seq(nTips) %in% X
    X <- X[tips]
    points(x=xx[X],y=yy[X], pch=c(x=4), cex=1.5)
  }

  if(!missing(O)){
    if(is.numeric(O)) O <- seq(nTips) %in% O
    O <- O[tips]
    points(x=xx[O], y=yy[O], pch=c(o=21), cex=2)
  }

  if(!is.null(largs)){
    groupTab <- table(groups)
    do.call(legend,
            c(list(
                   x=max(xx)*1.02,
                   y=mean(yy),
                   legend=gettextf('%s (%s)', names(groupTab), groupTab),
                   pch=glyphs$pch,
                   col=as.character(glyphs$col)
                   ),
              largs)
            )
  }
}
