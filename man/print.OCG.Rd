\name{print.OCG}
\alias{print.OCG}

\title{
Print a Summary of an \code{OCG} Object
}
\description{
This function prints summary statistics for an \code{OCG} object to the screen.
}
\usage{
\method{print}{OCG}(x, ...)
}

\arguments{
  \item{x}{
An object of class \code{OCG}.
}
  \item{\dots}{
Further arguments passed to or from other methods.
}
}

\value{
Prints summary data to the screen.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\examples{
## Generate graph and extract OCG communities.
g <- swiss[,3:4]
oc <- getOCG.clusters(g)

## Print summary statistics to the screen.
print(oc)
}

