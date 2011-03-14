svgTree <- function(x, file='tree.svg', groups,
                    desc1=groups,
                    desc2,
                    title='',
                    cex=rep(1, length(groups)),
                    border=rep(FALSE, length(groups)),
                    ...){
  ## * x - an object of class phylo, eg x <- nj(ddist)
  ## * file - name of output file
  ## * groups - a factor (or object coercible to a factor) assigning group identity
  ##   to leaf nodes in x
  ## * desc1 - optional first line of labels for toolTips;
  ##   if not provided, groups will be used.
  ## * desc2 - optional second line of labels for toolTips
  ## * title - title for the plot
  ## * cex - relative size of points
  ## * border - logical vector where TRUE indicates points
  ##   with a black border (ie col='black').
  ## * ... - additional arguments passed to plot.tree

  library(RSVGTipsDevice)
  
  if(missing(groups)){
    groups <- rep('leaf', length(x$tip.label))
  }
  
  groups <- factor(groups)

  ## TODO: make conditional desc2 cleaner
  ## check if groups, desc2 is length of groups
  if(missing(desc2)){
    desc2 <- rep(' ', length(groups))
  }
  
  glyphs <- clst:::getGlyphs(length(unique(groups)))

  col <- as.character(glyphs$col[groups])
  pch <- glyphs$pch[groups]

  tips <- getTips(x)

  devSVGTips(file=file,
             toolTipMode=2,
             title=title,
             width=10,
             height=10)

  ## par(xpd=FALSE)
  par(mai=c(bottom=0.5, left=0.5, top=0.5, right=2))
  plot(x, show.tip.label = FALSE, ...)
  nje <- get("last_plot.phylo", .PlotPhyloEnv)

  invisible(
            sapply(
                   tips,
                   function(i){
                     setSVGShapeToolTip(
                                        title=x$tip.label[i],
                                        desc1=as.character(desc1[i]),
                                        desc2=as.character(desc2[i])
                                        )
                     points(x=nje$xx[i],
                            y=nje$yy[i],
                            pch=pch[i],
                            col=ifelse(border, 'black', col[i]),
                            bg=col[i],
                            cex=cex[i])

                   }
                   )
            )
  invisible(dev.off())

}

svgTree2 <- function(x, file='tree.svg',
                     leaf.groups=rep('leaf',length(x$tip.label)),
                     leaves=list(
                         groups=leaf.groups,
                         desc1=leaf.groups,
                         desc2=rep('leaf', length(x$tip.label)),
                         cex=rep(1, length(x$tip.label))
                         ),
                     node.groups=rep('node',length(x$node.label)),
                     nodes=list(
                         groups=node.groups,
                         desc1=node.groups,
                         desc2=rep('node', length(x$node.label)),
                         cex=rep(1, length(x$node.label))
                         ),
                     title='',
                     ...){
  
  ## * x - an object of class phylo, eg x <- nj(ddist)
  ## * file - name of output file
  ## * groups - a factor (or object coercible to a factor) assigning group identity
  ##   to leaf nodes in x
  ## * desc1 - optional first line of labels for toolTips;
  ##   if not provided, groups will be used.
  ## * desc2 - optional second line of labels for toolTips
  ## * title - title for the plot
  ## * cex - relative size of points
  ## * ... - additional arguments passed to plot.tree

  library(RSVGTipsDevice)
  
  groups <- factor(c(leaves$groups,nodes$groups))
  leaves$groups <- factor(leaves$groups, levels=levels(groups))
  nodes$groups <- factor(nodes$groups, levels=levels(groups))
  
  glyphs <- clst:::getGlyphs(length(levels(groups)))

  col <- as.character(glyphs$col[groups])
  pch <- glyphs$pch[groups]

  devSVGTips(file=file,
             toolTipMode=2,
             title=title,
             width=10,
             height=10)

  ## par(xpd=FALSE)
  par(mai=c(bottom=0.5, left=0.5, top=0.5, right=2))
  plot(x, show.tip.label = FALSE, ...)
  nje <- get("last_plot.phylo", .PlotPhyloEnv)

  with(nodes,{
    invisible(
              sapply(
                     getNodes(x),
                     function(i){
                       ## convert to node index
                       n <- i - length(x$tip.label)                       
                       setSVGShapeToolTip(
                                          title=x$node.label[n],
                                          desc1=as.character(desc1[n]),
                                          desc2=as.character(desc2[n])
                                          )                       
                       points(x=nje$xx[i],
                              y=nje$yy[i],
                              pch=pch[i],
                              ## col=col[i],
                              col='black',
                              bg=col[i],
                              cex=cex[n]
                              )
                     }
                     )
              )
  })
  
  with(leaves,{
    invisible(
              sapply(
                     getTips(x),
                     function(i){
                       setSVGShapeToolTip(
                                          title=x$tip.label[i],
                                          desc1=as.character(desc1[i]),
                                          desc2=as.character(desc2[i])
                                          )
                       points(x=nje$xx[i],
                              y=nje$yy[i],
                              pch=pch[i],
                              ## col=col[i],
                              col='black',
                              bg=col[i],
                              cex=cex[i])
                     }
                     )
              )
  })
  
  invisible(dev.off())

}
