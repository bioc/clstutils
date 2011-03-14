\name{maxDists}
\alias{maxDists}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Select a maximally diverse set of items given a distance matrix.
}

\description{ Given a square matrix of pairwise distances, return
indices of N objects with a maximal sum of pairwise distances. }

\usage{
maxDists(mat, idx = NA, N = 1,
         exclude = rep(FALSE, nrow(mat)),
         include.center = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{mat}{square distance matrix}

\item{idx}{starting indices; if missing, starts with the object with the
maximum median distance to all other objects.}

\item{N}{total number of selections; length of idx is subtracted.}

\item{exclude}{boolean vector indicating elements to exclude from the
  calculation.}

\item{include.center}{includes the "most central" element (ie, the one with
  the smallest median of pairwise distances to all other elements) if TRUE}
}

% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
A vector of indices corresponding to the margin of \code{mat}.
}

% \references{
% %% ~put references to the literature/web site here ~
% }

\author{
Noah Hoffman
}

\note{
Note that it is important to evaluate if the candidate sequences contain
outliers (for example, mislabeled sequences), because these will
assuredly be included in a maximally diverse set of elements! 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{findOutliers}}
}

\examples{
library(ape)
library(clstutils)
data(seqs)
data(seqdat)
efaecium <- seqdat$tax_name == 'Enterococcus faecium'
seqdat <- subset(seqdat, efaecium)
seqs <- seqs[efaecium,]
dmat <- ape::dist.dna(seqs, pairwise.deletion=TRUE, as.matrix=TRUE, model='raw')

## find a maximally diverse set without first identifying outliers
picked <- maxDists(dmat, N=10)
picked
prettyTree(nj(dmat), groups=ifelse(1:nrow(dmat) \%in\% picked,'picked','not picked'))

## restrict selected elements to non-outliers
outliers <- findOutliers(dmat, cutoff=0.015)
picked <- maxDists(dmat, N=10, exclude=outliers)
picked
prettyTree(nj(dmat), groups=ifelse(1:nrow(dmat) \%in\% picked,'picked','not picked'),
X = outliers)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

