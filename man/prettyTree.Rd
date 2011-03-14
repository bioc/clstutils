\name{prettyTree}
\alias{prettyTree}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Draw an annotated phylogenetic tree. 
}

\description{ Extends \code{\link[ape]{plot.phylo}} to draw a
phylogenetic tree with additional annotation.}

\usage{

prettyTree(x, groups, fill,
           X, O, indices, labels,
           show = rep(TRUE, length(x)),
           largs = list(cex = 0.75, bty = "n", yjust = 0.5),
           parargs = list(mar = c(bottom = 5, left = 2, top = 2,
                              right = ifelse(is.null(largs), 2, 8)),
               xpd = NA),
           pointargs = list(), glyphs,
           shuffleGlyphs = NA, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{x}{an object of class phylo, eg \code{x <- nj(ddist)}}

\item{groups}{a factor (or object coercible) to a factor assigning group
identity to leaf nodes in x}

\item{fill}{vector (logical or indices) of points to fill}

\item{X}{vector of points to mark with an X}

\item{O}{vector of points to mark with a circle}

\item{indices}{label points with indices (all points if 'yes', or a
subset indicated by a vector)}

\item{labels}{character vector of tip labels in the same order as
\code{x$tip.label}}

\item{show}{boolean vector of points to show}

\item{largs}{arguments controlling appearance of the legend or NULL for
no legend}

\item{parargs}{arguments to pass par()}

\item{pointargs}{arguments to pass points() (other than pch, col, bg)}

\item{glyphs}{a \code{data.frame} with columns named 'col' and 'pch'
corresponding to elements of \code{unique(groups)}}

\item{shuffleGlyphs}{NA or an integer (argument to \code{set.seed})}

\item{\ldots}{passed to \code{\link[ape]{plot.phylo}}}

}

\details{
  
\code{prettyTree} adds to a plot drawn by \code{\link[ape]{plot.phylo}}

Vectors specifying annotation should be in the order of row or column
labels of the distance matrix used to generate \code{x}.

}

\value{
Plots to the active device; no return value.
}

% \references{
% %% ~put references to the literature/web site here ~
% }

\author{
Noah Hoffman
}

\note{
See package vignette for additional examples.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[ape]{plot.phylo}}
}
\examples{
library(ape)
data(seqs)
data(seqdat)
prettyTree(nj(dist.dna(seqs)), groups=seqdat$tax_name)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{aplot}

