\name{classifyPlacements}
\alias{classifyPlacements}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Taxonomic classification by phylogenetic placement.
}

\description{

  Given taxonomic information from a reference package and inter-node
  distances from a reference tree, perform classification of one or more
  placements provided by \command{pplacer}.

}

\usage{
classifyPlacements(taxdata, treedists, placetab, ...,
    verbose = FALSE, debug = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{taxdata}{\code{data.frame}, output of \code{\link{taxonomyFromRefpkg}} }

  \item{treedists}{output of \code{\link{treeDists}}}

  \item{placetab}{a data.frame with columns \code{at}, \code{edge}, and
\code{branch}}

  \item{\dots}{extra arguments passed to \code{\link[clst]{classifyIter}}}

  \item{verbose}{writes progress messages to terminal if TRUE}
  
  \item{debug}{be very verbose if TRUE}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{

The output is a \code{data.frame} describing the taxonomic assignment,
along with a description of the confidence of the classification. See
the man page for \code{\link[clst]{classify}} for details on the output.

}
% \references{

% }
\author{
Noah Hoffman
}
% \note{
% %%  ~~further notes~~
% }

\seealso{
\code{\link{treeDists}}, \code{\link{taxonomyFromRefpkg}}, \code{\link{placeData}}
}

\examples{
placefile <- system.file('extdata', 'qalign.place.bz2', package='clstutils')
distfile <- system.file('extdata', 'qalign.distmat.bz2', package='clstutils')
refpkgz <- system.file('extdata', 'vaginal_16s.refpkg.tar.bz2', package='clstutils')

tmpdir <- tempdir()
system(paste('tar -xf', refpkgz, '-C', tmpdir))
refpkg <- file.path(tmpdir, 'vaginal_16s.refpkg')

treedists <- treeDists(distfile=distfile, placefile=placefile)
taxdata <- taxonomyFromRefpkg(refpkg, seqnames=rownames(treedists$dmat), lowest_rank='species')
placetab <- data.frame(at=49, edge=5.14909e-07, branch=5.14909e-07)
result <- classifyPlacements(taxdata, treedists, placetab)
result
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

