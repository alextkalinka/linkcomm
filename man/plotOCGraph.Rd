\name{plotOCGraph}
\alias{plotOCGraph}

\title{
Plot a Graph Layout of OCG Communities
}
\description{
This function is called by \code{plot.OCG} to plot a graph layout of the OCG communities.
}
\usage{
plotOCGraph(x, clusterids = 1:x$numbers[3], nodes = NULL, pie.local = TRUE, 
        incident = TRUE, layout = layout.fruchterman.reingold, 
        vertex.radius = 0.03, scale.vertices = 0.05, edge.color = "grey",
        vertex.label.color = "black", vertex.label.cex = 0.8,
        pal = brewer.pal(7,"Set2"), shownodesin = 0, vlabel = TRUE, 
        random = TRUE, ...)
}

\arguments{
  \item{x}{
An object of class \code{OCG}.
}
  \item{clusterids}{
An integer vector of community IDs. Defaults to all communities.
}
  \item{nodes}{
A character vector naming the nodes to be plotted. If NULL, then community IDs are used instead. Defaults to NULL.
}
  \item{pie.local}{
Logical, whether to include pie segments for just the communities to which the chosen nodes belong, or for all communities. In the latter case, communities that are not present among the chosen nodes will appear as empty segments in the node pies. Defaults to TRUE.
}
  \item{incident}{
Logical, whether to include just the communities of the named node(s), or the community membership of all nodes that interact with the named node(s). Defaults to TRUE.
}
  \item{layout}{
A character string or function identifying the layout algorithm to be used for positioning nodes in the graph. Defaults to \code{layout.fruchterman.reingold}. See details for alternative layouts.
}
\item{vertex.radius}{
A real number specifying the vertex radius. Defaults to 0.03.
}
\item{scale.vertices}{
A real number specifying the fraction of increase in vertex radius for each community membership. If NULL then all vertices are the same size. Defaults to 0.05.
}
\item{edge.color}{
A character string specifying the colour of edges. Defaults to \code{"grey"}.
}
  \item{vertex.label.color}{
A character string specifying the color of node labels. Defaults to \code{"black"}.
}
  \item{vertex.label.cex}{
A numerical value specifying the size of the node labels. Defaults to 0.8.
}
  \item{pal}{
A character vector describing a colour palette to be used for colouring the link communities in the graph. Defaults to \code{brewer.pal(7, "Set2")}.
}
  \item{shownodesin}{
An integer value specifying the number of communities a node must belong to before it will be displayed. If 0 then all nodes are displayed. Defaults to 0.
}
  \item{vlabel}{
Logical, whether node labels are to be added. Defaults to TRUE.
}
  \item{random}{
Logical, whether to randomise the link colours. Defaults to TRUE.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}
\details{
Here we describe the parameters for plotting OCG community graphs using:\cr
\code{plot(x, type = "graph", layout = layout)}\cr\cr
Various graph layouts are available:
\enumerate{
\item \code{layout.random}
\item \code{layout.circle}
\item \code{layout.sphere}
\item \code{layout.fruchterman.reingold}
\item \code{layout.kamada.kawai}
\item \code{layout.spring}
\item \code{layout.reingold.tilford}
\item \code{layout.fruchterman.reingold.grid}
\item \code{layout.lgl}
\item \code{layout.graphopt}
\item \code{layout.mds}
\item \code{layout.svd}
\item \code{layout.norm}
}
All of these are described in more detail in the \code{\link[igraph]{igraph}} package.
}
\value{
A graph plot.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{plot.OCG}}, \code{\link[igraph]{igraph.plotting}}
}
\examples{
## Generate graph and extract OCG communities.
g <- swiss[,3:4]
oc <- getOCG.clusters(g)

## Plot a graph of OCG communities.
plot(oc, type = "graph")
}

