\name{getEdgesIn}
\alias{getEdgesIn}

\title{
Extract Edge Indices from Communities
}
\description{
This function returns edge indices that belong to certain communities or that are incident upon certain nodes.
}
\usage{
getEdgesIn(x, clusterids = 1, nodes = NULL, all = FALSE)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm} or \code{OCG}.
}
  \item{clusterids}{
An integer vector of community IDs. Defaults to community 1.
}
  \item{nodes}{
A character vector specifying node(s) for which edge indices should be returned. Overrides \code{clusterids}. Defaults to NULL.
}
  \item{all}{
Logical, whether the edges for all communities to which the named nodes belong should be returned. Will have an effect only if \code{nodes} is not NULL. If FALSE, edges that are directly incident upon the named nodes will be returned. Defaults to FALSE.
}
}

\value{
An integer vector of edge indices.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Get edges from community 1.
getEdgesIn(lc)
}

