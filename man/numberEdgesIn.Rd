\name{numberEdgesIn}
\alias{numberEdgesIn}

\title{
Extract Node Community Membership by Edges
}
\description{
This function returns the number of edges node(s) contain in each community.
}
\usage{
numberEdgesIn(x, clusterids = 1:x$numbers[3], nodes)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm} or \code{OCG}.
}
  \item{clusterids}{
An integer vector of community IDs. Defaults to all communities.
}
  \item{nodes}{
A character vector specifying node(s) for which edge membership should be returned.
}
}

\value{
A named list of named integer vectors specifying the number of edges in each community a node belongs in. Names of the integer vectors are community IDs, and names of the list are node names.
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

## Get edges from community 1.
numberEdgesIn(oc, nodes = 1)
}

