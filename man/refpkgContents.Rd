\name{refpkgContents}
\alias{refpkgContents}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Read the contents of a collection of reference sequences ("refpkg").
}

\description{
Read the manifest file from a refpackage and
return a list containing the package contents.
}

\usage{
refpkgContents(path, manifest = "CONTENTS.json")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{path to a refpkg directory}

  \item{manifest}{name of the manifest file}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
Returns a list of lists. Run \code{example(refpkgContents)} for details.
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
\code{\link{taxonomyFromRefpkg}}
}
\examples{
archive <- 'vaginal_16s.refpkg.tar.gz'
destdir <- tempdir()
system(gettextf('tar -xzf \%s --directory="\%s"',
                system.file('extdata',archive,package='clstutils'),
                destdir))
refpkg <- file.path(destdir, sub('.tar.gz','',archive))
contents <- refpkgContents(refpkg)
str(refpkg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

