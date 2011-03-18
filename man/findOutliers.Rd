\name{findOutliers}
\alias{findOutliers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify outlier objects given a square distance matrix.
}

\description{
  Outliers are defined as elements with edge length to the
  centermost element > cutoff. The distance threshold (cutoff) can be either
  specified, or calculated as a quantile of all pairwise distances in
  the matrix.
}

\usage{
findOutliers(mat, quant, cutoff)
}

\arguments{
  \item{mat}{square matrix of distances}

  \item{quant}{given all pairwise distances x, calculate distance
  threshold as quantile(x, quant).  Values closer to 0 are more
  stringent.}

  \item{cutoff}{an absolute cutoff overriding quant}

}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{

  Returns a boolean vector corresponding to margin of mat; outliers have
  a value of TRUE.
  
}
% \references{
% %% ~put references to the literature/web site here ~
% }
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
library(ape)
data(seqs)
data(seqdat)
dmat <- ape::dist.dna(seqs[seqdat$tax_name == 'Enterococcus faecium',],
        pairwise.deletion=TRUE, as.matrix=TRUE, model='raw')
summary(dmat[lower.tri(dmat)])
outliers <- findOutliers(dmat, cutoff=0.015)
table(outliers)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}
