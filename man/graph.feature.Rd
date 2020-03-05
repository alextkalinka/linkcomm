\name{graph.feature}
\alias{graph.feature}

\title{
Make Node or Edge Graph Features
}
\description{
This function returns vectors of node sizes or edge widths for use in \code{plot.linkcomm}.
}
\usage{
graph.feature(x, type = "nodes", clusterids = 1:length(x$clusters), 
              nodes = NULL, indices, features, default = 15, showall = FALSE)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{type}{
A character string specifying either \code{"nodes"} or \code{"edges"}.
}
  \item{clusterids}{
An integer vector of community IDs that will be plotted. Defaults to all communities.
}
  \item{nodes}{
A character vector specifying node(s) that will be plotted. Overrides \code{clusterids}. Defaults to NULL.
}
  \item{indices}{
An integer vector specifying the indices of the nodes or edges that will be given specific size or width values. See \code{\link{getNodesIn}} and \code{\link{getEdgesIn}} for ways to generate these indices. Also see examples in \code{vignette(topic = "linkcomm", package = "linkcomm")}.
}
  \item{features}{
An integer vector specifying the node or edge sizes for the nodes or edges that are to be changed. If there is a single value then this will be applied to all nodes or edges specified in \code{indices}, otherwise the \code{features} vector must be the same length as the \code{indices} vector and the values will be matched to each other.
}
  \item{default}{
An integer value specifying the node size or edge width that all nodes or edges not specified by \code{indices} will take. Defaults to 15.
}
  \item{showall}{
Logical, whether edges that don't belong to communities will also be plotted or not. Defaults to FALSE.
}
}

\value{
A named integer vector of node sizes or edge widths. The names will be either node names or edge indices.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{plotLinkCommGraph}}, \code{\link{getNodesIn}}, \code{\link{getEdgesIn}}, \code{vignette(topic = "linkcomm", package = "linkcomm")}.
}

\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Make node size vector for all nodes that belong to community 1.
graph.feature(lc, indices = getNodesIn(lc, type = "indices"), features = 20, default = 5)

## Make edge width vector for all edges that belong to community 1.
graph.feature(lc, type = "edges", indices = getEdgesIn(lc), features = 5, default = 1)
}

