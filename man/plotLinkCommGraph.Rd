\name{plotLinkCommGraph}
\alias{plotLinkCommGraph}

\title{
Plot a Graph Layout of Link Communities
}
\description{
This function is called by \code{plot.linkcomm} to plot a graph layout of the link communities.
}
\usage{
plotLinkCommGraph(x, clusterids = 1:length(x$clusters), nodes = NULL,
        layout = layout.fruchterman.reingold, pal = brewer.pal(7, "Set2"),
        random = TRUE, node.pies = TRUE, pie.local = TRUE, vertex.radius = 0.03,
        scale.vertices = 0.05, edge.color = NULL, vshape = "none", vsize = 15,
        ewidth = 3, margin = 0, vlabel.cex = 0.8, vlabel.color = "black",
        vlabel.family = "Helvetica", vertex.color = "palegoldenrod",
        vlabel = TRUE, col.nonclusters = "black", jitter = 0.2, circle = TRUE,
        printids = TRUE, cid.cex = 1, shownodesin = 0, showall = FALSE,
        verbose = TRUE, ...)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{clusterids}{
An integer vector of community IDs. Defaults to all communities.
}
  \item{nodes}{
A character vector naming the nodes to be plotted. If NULL, then community IDs are used instead. Defaults to NULL.
}
  \item{layout}{
A character string or function identifying the layout algorithm to be used for positioning nodes in the graph. Defaults to \code{layout.fruchterman.reingold}. See details for alternative layouts.
}
  \item{pal}{
A character vector describing a colour palette to be used for colouring the link communities in the graph. Defaults to \code{brewer.pal(7, "Set2")}.
}
  \item{random}{
Logical, whether to randomise the link colours. Defaults to TRUE.
}
  \item{node.pies}{
Logical, whether to plot node pies showing as a pie chart the fraction of each node's edges which belong to each of its communities. Defaults to TRUE.
}
  \item{pie.local}{
Logical, whether to include pie segments for just the communities to which the chosen nodes belong, or for all communities. In the latter case, communities that are not present among the chosen nodes will appear as empty segments in the node pies. Will be ignored if node.pies is FALSE. Defaults to TRUE.
}
\item{vertex.radius}{
A real number specifying the vertex radius. Defaults to 0.03. Will be ignored if node.pies is FALSE.
}
\item{scale.vertices}{
A real number specifying the fraction of increase in vertex radius for each community membership. If NULL then all vertices are the same size. Defaults to 0.05. Will be ignored if node.pies is FALSE.
}
\item{edge.color}{
A character string specifying the colour of edges. Defaults to \code{"grey"}. Will be ignored if node.pies is FALSE.
}
  \item{vshape}{
A character string specifying the shape of the nodes. Can be one of \code{"none"}, \code{"circle"}, \code{"square"}, \code{"csquare"}, \code{"rectangle"}, \code{"crectangle"}, and \code{"vrectangle"}. Defaults to \code{"none"}.
}
  \item{vsize}{
An integer vector of node sizes. If there is a single value this will be used for all nodes. If there are multiple values, it must be the same length as the number of nodes in the network to be visualized. See \code{\link{graph.feature}} for generating these values. This argument only has an effect when \code{vshape} is not set to \code{"none"}. Defaults to 15.
}
  \item{ewidth}{
An integer vector of edge widths. If there is a single value this will be used for all edges. If there are multiple values, it must be the same length as the number of edges in the network to be visualized. See \code{\link{graph.feature}} for generating these values. Defaults to 3.
}
  \item{margin}{
A numerical value specifying the amount of empty space around the graph. Negative values will zoom into the graph. Defaults to 0.
}
  \item{vlabel.cex}{
A numerical value specifying the size of node labels. Defaults to 0.8.
}
  \item{vlabel.color}{
A character string specifying the color of node labels. Defaults to \code{"black"}.
}
  \item{vlabel.family}{
A character string specifying the font family for node labels. Defaults to \code{"Helvetica"}.
}
  \item{vertex.color}{
A character string specifying the colour of nodes. If this is a character vector then the colours will be recycled. Defaults to \code{"palegoldenrod"}.
}
  \item{vlabel}{
Logical, whether node labels are to be added. Defaults to TRUE.
}
  \item{col.nonclusters}{
A character string specifying the colour of edges that do not belong to any communities. Will only have an effect if \code{showall} is TRUE. Defaults to \code{"black"}.
}
  \item{jitter}{
A numerical value specifying the range (negative to positive) of random noise that will be added to nodes that have identical x-y coordinates. Defaults to 0.2. Only used for Spencer circle layouts.
}
  \item{circle}{
Logical, whether to display a circle for a Spencer circle layout. Defaults to TRUE.
}
  \item{printids}{
Logical, whether to display community IDs at their anchor points around the Spencer circle. Defaults to TRUE.
}
  \item{cid.cex}{
A numerical value specifying the size of community IDs around the Spencer circle. Defaults to 1.
}
  \item{shownodesin}{
An integer value specifying the number of communities a node must belong to before it will be displayed. If 0 then all nodes are displayed. Defaults to 0.
}
  \item{showall}{
Logical, whether to display all links in the network regardless of whether they belong to communities or not. Defaults to FALSE.
}
  \item{verbose}{
Logical, whether to print the progress of the calculation to the screen. Defaults to TRUE.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}
\details{
Here we describe the parameters for plotting link community graphs using:\cr
\code{plot(x, type = "graph", layout = layout)}\cr\cr
Various graph layouts are available:
\enumerate{
\item \code{"spencer.circle"}
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
All of these, except the \code{"spencer.circle"}, are described in more detail in the \code{\link{igraph}} package. The \code{"spencer.circle"} is described in \code{\link{layout.spencer.circle}}.
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
\code{\link{plot.linkcomm}}, \code{\link{layout.spencer.circle}}, \code{\link{graph.feature}}, \code{\link{igraph.plotting}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Plot a graph of link communities.
plot(lc, type = "graph")

## Plot a graph of link communities using a Spencer circle layout.
plot(lc, type = "graph", layout = "spencer.circle")
}

