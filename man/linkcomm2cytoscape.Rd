\name{linkcomm2cytoscape}
\alias{linkcomm2cytoscape}

\title{
Write an Edge Attribute File for Cytoscape
}
\description{
This function writes out an edge attribute file for visualising the link communities in Cytoscape.
}
\usage{
linkcomm2cytoscape(x, interaction = "pp", ea = "temp.ea")
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{interaction}{
A character string indicating the type of interaction between nodes. Defaults to \code{"pp"} for protein-protein interaction.
}
  \item{ea}{
A character string indicating the file for writing the edge attributes. Defaults to \code{"temp.ea"}.
}
}
\details{
Cytoscape is an open source platform for complex-network analysis and visualization (Shannon et al. 2003).
}
\value{
Used for its side-effect of writing an edge attribute file to disk.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.

Shannon, P. et al. (2003) Cytoscape: A software environment for integrated models of biomolecular interaction networks. \emph{Genome Research} \bold{13}, 2498-2504.

}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Write an edge attribute file to disk.
linkcomm2cytoscape(lc)
}

