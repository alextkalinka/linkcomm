\name{get.community.overlaps}
\alias{get.community.overlaps}

\title{
Find Node Overlaps between Communities
}
\description{
This function returns lists of communities that share at least one node with each other.
}
\usage{
get.community.overlaps(x)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm} or \code{OCG}.
}
}

\value{
A list of length equal to the number of communities. Each element contains an integer vector giving the community IDs for communities sharing at least one node with each community. NAs indicate that a community shares no nodes with any other communities.
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

## Get list of overlapping communities.
ov <- get.community.overlaps(lc)
}

