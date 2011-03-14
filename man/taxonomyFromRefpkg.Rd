\name{taxonomyFromRefpkg}
\alias{taxonomyFromRefpkg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract taxonomic information from a refpkg.
}
\description{

  Construct a data.frame providing the lineage of each sequence
  represented in the reference package.

}
\usage{
taxonomyFromRefpkg(path, seqnames, lowest_rank = NA)
}


\arguments{

\item{path}{path to a refpkg directory}

\item{seqnames}{optional character vector of sequence names. If
provided, determines the order of rows in \code{$taxTab}}

\item{lowest_rank}{name of the most specific (ie, rightmost) rank to
 include. Default is the name of the rightmost column in
 \code{refpkg_contents$taxonomy}}

}

% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }

\value{
A list with the following elements:
  
\item{taxNames}{a named character vector of taxonomic names (names are tax_ids)}

\item{taxTab}{a \code{data.frame} in which each row corresponds to a
  reference sequence and contains a tax_id followed by the corresponding
  lineage (columns are "root"...lowest_rank)}

}

\references{

    The decsription and specification for a reference package can be found
  in the project repository in github:
  \url{https://github.com/fhcrc/taxtastic}

  Scripts and tools for creating reference packages are provided in the
  python package \command{taxonomy}, also available from the taxtastic
  project site.

}
\author{

Noah Hoffman
  
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{refpkgContents}}
}
\examples{
archive <- 'vaginal_16s.refpkg.tar.bz2'
destdir <- tempdir()
system(gettextf('tar -xjf \%s --directory="\%s"',
                system.file('extdata',archive,package='clstutils'),
                destdir))
refpkg <- file.path(destdir, sub('.tar.bz2','',archive))
reftax <- taxonomyFromRefpkg(refpkg)
str(reftax)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

