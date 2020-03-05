\name{getAllNestedComm}
\alias{getAllNestedComm}

\title{
Find Nested Communities
}
\description{
This function returns communities of nodes that are entirely nested within other larger communities of nodes.
}
\usage{
getAllNestedComm(x, verbose = FALSE, plot = FALSE)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{verbose}{
Logical, whether to print to the screen a warning that individual community IDs are not clustered in any other communities. Defaults to FALSE.
}
  \item{plot}{
Logical, whether to plot graphs of the nested communities. Defaults to FALSE.
}
}
\details{
Nested community structures may reveal interesting relationships among sets of nodes.
}
\value{
A named list of integer vectors; names are integers referring to nested communities, and the integer vectors are the communities that the named community is nested in.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{getNestedHierarchies}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Find nested communities.
getAllNestedComm(lc)
}

