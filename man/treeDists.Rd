\name{treeDists}
\alias{treeDists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Provide objects for determining distances among nodes of a reference tree.  
}
\description{

Provides objects (dists, paths) that can be used to calculate
vectors of distances between an internal node and each leaf
node. Also returns a square matrix of distances between leaf
nodes.

}
\usage{
treeDists(distfile, placefile, pmat, pmap)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

\item{distfile}{path to output of \command{placeutil --distmat}}

\item{placefile}{path to \command{pplacer} output}
  
\item{pmat}{square matrix (output of \code{clstutils:::placeMat}); may be
provided in the place of \code{pmat}}

\item{pmap}{output of \code{clstutils:::placeMap} (data.frame); may be
provided in the place of \code{pmap}}

}
\details{

 A placement on an edge looks like this:
 \preformatted{
 proximal
  |
  |   d_p
  |
  |---- x
  |
  |   d_d
  |
  |
 distal}
 d_p is the distance from the placement x to the proximal side of the
 edge, and d_d the distance to the distal side.

 If the distance from x to a leaf y is an S-distance Q, then the
 path from x to y will go through the distal side of the edge and
 we will need to add d_d to Q to get the distance from x to y.  If
 the distance from x to a leaf y is a P-distance Q, then the path
 from x to y will go through the proximal side of the edge, and we
 will need to subtract off d_d from Q to get the distance from x
 to y. In either case, we always need to add the
 length of the pendant edge, which is the second column.

 To review, say the values of the two leftmost columns are a and b
 for a given placement x, and that it is on an edge i.  We are
 interested in the distance of x to a leaf y, which is on edge j.
 We look at the distance matrix, entry (i,j), and say it is an
 S-distance Q. Then our distance is Q+a+b.  If it is a P-distance
 Q, then the distance is Q-a+b.

 The distances between leaves should always be P-distances, and there
 we need no trickery.

 (thanks to Erick Matsen for this description)
 
}
\value{

A list with the following elements:

\item{dists}{rectangular matrix of distances with rows corresponding to
  all nodes in pplacer order, and columns corresponding to tips in the
  order of the corresponding \code{phylo\{ape\}} object.}

\item{paths}{rectangular matrix in the same configuration as
  \code{dists} with values of 1 or -1 if the path between nodes is
  serial or parallel, respectively (see Details)}

\item{dmat}{square matrix containing distances between pairs of tips.}

}
\references{
  
Documentation for \command{pplacer} can be found here:
\url{http://matsen.fhcrc.org/pplacer/}

}

\author{
Noah Hoffman
}

\note{
The output of this function is required for \code{\link{classifyPlacements}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{classifyPlacements}}  
}
\examples{

placefile <- system.file('extdata','qalign.place.bz2', package='clstutils')
distfile <- system.file('extdata','qalign.distmat.bz2', package='clstutils')
treedists <- treeDists(distfile=distfile, placefile=placefile)

## coordinates of a single placement
placetab <- data.frame(at=49, edge=5.14909e-07, branch=5.14909e-07)

## dvects is a matrix in which each row corresponds to a vector of
## distances between a single placement along the edge of the reference
## tree used to generate 'distfile', and each column correspons to a
## reference sequence (ie, a terminal node).

dvects <- with(placetab, {
treedists$dists[at+1,,drop=FALSE] + treedists$paths[at+1,,drop=FALSE]*edge + branch
})

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

