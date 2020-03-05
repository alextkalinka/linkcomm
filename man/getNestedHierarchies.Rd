\name{getNestedHierarchies}
\alias{getNestedHierarchies}

\title{
Find Nested Structures in Communities
}
\description{
This function determines whether a particular community is nested within any other communities.
}
\usage{
getNestedHierarchies(x, clusid = 1, verbose = TRUE, plot = TRUE, ids = FALSE)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{clusid}{
An integer value indicating the community ID whose nesting structure will be tested. Defaults to 1.
}
  \item{verbose}{
Logical, whether to display a warning that a particular community is not nested in any other communities on the screen. Defaults to FALSE.
}
  \item{plot}{
Logical, whether to plot a graph layout of the nested community.
}
  \item{ids}{
Logical, whether to return only the community IDs that the community is nested in, or the node names also. Defaults to FALSE.
}
}

\value{
Either a list of character vectors, each giving the nodes that the community is nested in, or an integer vector of community IDs that the community is nested in.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{getAllNestedComm}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Determine if community 1 is nested in any other communities.
getNestedHierarchies(lc, clusid = 1)
}

