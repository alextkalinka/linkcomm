\name{plot.linkcomm}
\alias{plot.linkcomm}

\title{
The \code{linkcomm} Plotting Function
}
\description{
This function plots various different \code{linkcomm} graphs.
}
\usage{
\method{plot}{linkcomm}(x, type = "", ...)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{type}{
A character string specifying the type of plot. Can be one of \code{"summary"}, \code{"members"}, \code{"graph"}, \code{"commsumm"}, and \code{"dend"}. See Details below.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}
\details{
\code{"summary"} plots the dendrogram and partition density plot side-by-side;\cr \code{"members"} plots a community membership matrix;\cr \code{"graph"} plots a graph layout of the network with coloured link communities;\cr \code{"commsumm"} plots a bar graph or pie chart summarising community modularity or connectedness for each community;\cr \code{"dend"} plots a dendrogram with coloured link communities.\cr\cr See the individual plotting functions for details of arguments that can be passed to \code{plot.linkcomm}: \code{\link{plotLinkCommSumm}}, \code{\link{plotLinkCommMembers}}, \code{\link{plotLinkCommGraph}}, \code{\link{plotLinkCommSummComm}}, and \code{\link{plotLinkCommDend}}.
}
\value{
Plots to the current device.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{plotLinkCommSumm}}, \code{\link{plotLinkCommMembers}}, \code{\link{plotLinkCommGraph}}, \code{\link{plotLinkCommSummComm}}, \code{\link{plotLinkCommDend}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Plot a graph of link communities.
plot(lc, type = "graph")
}

