\name{plot.OCG}
\alias{plot.OCG}

\title{
The \code{OCG} Plotting Function
}
\description{
This function plots various different \code{OCG} graphs.
}
\usage{
\method{plot}{OCG}(x, type = "", ...)
}

\arguments{
  \item{x}{
An object of class \code{OCG}.
}
  \item{type}{
A character string specifying the type of plot. Can be one of \code{"members"} or \code{"graph"}. See Details below.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}
\details{
\code{"members"} plots a community membership matrix;\cr \code{"graph"} plots a graph layout of the network with coloured link communities.\cr\cr See the \code{OCG} plotting function for details of arguments that can be passed to \code{plot.OCG}: \code{\link{plotOCGraph}}, \code{\link{plotLinkCommMembers}}.
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
\code{\link{plotOCGraph}}, \code{\link{plotLinkCommMembers}}
}
\examples{
## Generate graph and extract OCG communities.
g <- swiss[,3:4]
oc <- getOCG.clusters(g)

## Plot a graph of OCG clusters.
plot(oc, type = "graph")
}

