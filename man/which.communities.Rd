\name{which.communities}
\alias{which.communities}

\title{
Extract Community Membership for Nodes
}
\description{
This function returns the community IDs of the communities to which one or more nodes belong.
}
\usage{
which.communities(x, nodes)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm} or \code{OCG}.
}
  \item{nodes}{
A character vector specifying the nodes.
}
}

\value{
An integer vector of community IDs.
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
which.communities(oc, nodes = 1)
}

