\name{orderCommunities}
\alias{orderCommunities}

\title{
Order Link Communities According to the Dendrogram
}
\description{
This function returns link communities in the same order as in the hierarchical clustering dendrogram.
}
\usage{
orderCommunities(x, clusterids = 1:x$numbers[3], verbose = TRUE)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{clusterids}{
An integer vector of community IDs. Defaults to all communities.
}
  \item{verbose}{
Logical, whether to print progress of the calculation to the screen. Defaults to TRUE.
}
}
\details{
Ordering link communities according to the dendrogram can aid in visualization when plotting them as a Spencer circle because it minimses crossing over between links.
}
\value{
A list with the following components:
\item{ordered }{A list of integer vectors. These are the ordered communities of links.}
\item{clusids }{An integer vector of community IDs in their new order.}
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{plot.linkcomm}}, \code{\link{plotLinkCommGraph}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Order communities according to the dendrogram.
orderCommunities(lc)
}

