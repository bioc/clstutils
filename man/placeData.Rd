\name{placeData}
\alias{placeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read output file produced by \command{pplacer}
}
\description{
Given a file produced by \command{pplacer}, return a \code{data.frame}
of placement location and statistics.
}

\usage{
placeData(file)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{Name of the input file.}

}
\details{

 Columns in the data.frame output include 'name' (name of the
 sequence), 'hit' (integer corresponding to each row of output for a
 given sequence in the pplacer file); plus the following columns as
 described in (http://matsen.fhcrc.org/pplacer/manual.html#placefile):

 \enumerate{
 \item Placement edge number, numbered as in a postorder traversal
 \item ML likelihood weight ratio (i.e. the normalized ML likelihood values)
 \item Posterior probability (or a dash if the -p option was not set)
 \item ML log likelihood
 \item Bayes marginal likelihood (or a dash if the -p option was not set)
 \item The ML distance from the distal (farthest from the root) side of the edge
 \item The ML pendant branch length
 \item The taxon id assigned by pplacer
}
}

\value{ A data.frame with columns
  'name','hit','at','mlwr','ppost','mlll','bml','edge','branch','tax_id' (see
  Details) }

\references{
See the manual for \command{pplacer} at \url{http://matsen.fhcrc.org/pplacer/}
}

\author{
Noah Hoffman
}

% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }

\examples{
placefile <- system.file('extdata','qalign.place.bz2', package='clstutils')
str(placedata <- placeData(placefile))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

