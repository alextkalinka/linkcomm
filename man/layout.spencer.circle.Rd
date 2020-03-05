\name{layout.spencer.circle}
\alias{layout.spencer.circle}

\title{
Calculate Node Coordinates for a Spencer Circle
}
\description{
This function returns the x-y coordinates for nodes in a Spencer circle together with community anchor positions.
}
\usage{
layout.spencer.circle(x, clusterids = 1:x$numbers[3], verbose = TRUE, 
                      jitter = 0.2)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{clusterids}{
An integer vector of community IDs. Defaults to all communities.
}
  \item{verbose}{
Logical, whether to print the progress of the calculation to the screen. Defaults to TRUE.
}
  \item{jitter}{
A positive numerical value specifying the range (negative to positive) of random, uniformly distributed noise that will be added to nodes that have identical x-y coordinates. Defaults to 0.2.
}
}
\details{
This algorithm anchors communities evenly around the circumference of a circle in their dendrogram order (to minimise crossing over of links) and positions nodes within the circle according to how many links they possess in each of the communities (Spencer, 2010). Thus, nodes that have links to a lot of communities will get pushed into the centre of the circle making this method well suited for representing ego networks where one or a small number of nodes belong to multiple communities.
}
\value{
A list with the following components:
\item{nodes }{A numerical matrix with nodes as rows and with 2 columns; the first contains the x coordinates and the second the y coordinates.}
\item{anchors }{A numerical matrix with communities as rows and with 2 columns of x and y coordinates.}
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.

Spencer, R. (2010). \url{http://scaledinnovation.com/analytics/communities/comlinks.html}
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

## Extract x-y coordinates for nodes in a Spencer circle.
layout.spencer.circle(lc)
}

