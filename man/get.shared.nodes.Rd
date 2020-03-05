\name{get.shared.nodes}
\alias{get.shared.nodes}

\title{
Get Nodes Shared by Communities
}
\description{
This function returns the nodes that are shared by specific sets of communities.
}
\usage{
get.shared.nodes(x, comms)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm} or \code{OCG}.
}
  \item{comms}{
An integer vector giving the community IDs for which an intersecting set of shared nodes should be returned.
}
}

\value{
A character vector giving the shared node names.
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

## Get shared nodes for communities 1 and 7.
get.shared.nodes(lc, comms = c(1,7))
}

