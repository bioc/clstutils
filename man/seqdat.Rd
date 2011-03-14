\name{seqdat}
\alias{seqdat}
\docType{data}
\title{
Annotation for the Enterococcus sequence data set.
}

\description{ Provides annotation for \code{link{seqs}}, an aligned 16S
rRNA sequences representing three Enterococcus species.  }

\usage{data(seqdat)}
\format{
  A data frame with 200 observations on the following 5 variables.
  \describe{
    \item{\code{seqname}}{a character vector}
    \item{\code{accession}}{a character vector containing GenBank
      accession numbers.}
    \item{\code{tax_id}}{a character vector}
    \item{\code{tax_name}}{a character vector}
    \item{\code{isType}}{a logical vector indicating if the sequence is
      from a type strain.}
  }
}

% \details{

% }

\source{
These sequences were downloaded from the Ribosomal Database Project
website \url{http://rdp.cme.msu.edu/}
}

% \references{
% %%  ~~ possibly secondary sources and usages ~~
% }
\examples{
data(seqdat)
with(seqdat,{
table(tax_name, isType)
})
}
\keyword{datasets}
