\name{LinkDensities}
\alias{LinkDensities}

\title{
Calculate Link Community Link Densities
}
\description{
This function calculates link densities for link communities.
}
\usage{
LinkDensities(x, clusterids = 1:x$numbers[3])
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{clusterids}{
An integer vector of community IDs. Defaults to all communities.
}
}
\details{
The link density of community i is
\deqn{D_{i}=\frac{e_{i}-n_{i}+1}{(n_{i}(n_{i}-1)/2)-n_{i}+1}}{D_{i} = (e_{i}-n_{i}+1)/((n_{i}(n_{i}-1)/2)-n_{i}+1)}
where \eqn{e_{i}} is the number of edges in community i and \eqn{n_{i}} is the number of nodes in community i.
}
\value{
A named numerical vector, where the names are community IDs and the numbers are link densities.
}
\references{
Ahn, Y.Y., Bagrow, J.P., and Lehmann, S. (2010). Link communities reveal multiscale complexity in networks. \emph{Nature} \bold{466}, 761-764.

Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{plot.linkcomm}}, \code{\link{plotLinkCommSummComm}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Calculate link densities.
ld <- LinkDensities(lc)
}

