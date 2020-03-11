\name{linkcomm2clustnsee}
\alias{linkcomm2clustnsee}

\title{
Write a Partition File for Clust&See
}
\description{
This function writes out a partition file which can be imported into the Cytoscape plug-in Clust&See.
}
\usage{
linkcomm2clustnsee(x, file = "temp.cns", network.name = NULL)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm} or \code{OCG}.
}
  \item{file}{
A character string naming a Clust&See partition file (.cns extension). Defaults to \code{"temp.cns"}.
}
  \item{network.name}{
A character string providing a name for the network. This name must correspond to the file name of the network that will be imported into Cytoscape. If NULL, the object name, \code{x}, is used. Defaults to NULL.
}
}
\details{
Cytoscape is an open source platform for complex-network analysis and visualization, and Clust&See (Spinelli et al. 2013) is a Cytoscape plug-in used for visualizing and manipulating the clusters produced by various network clustering algorithms.
}
\value{
Used for its side-effect of writing a Clust&See partition file to disk.
}
\references{
Becker, E. et al. (2012) Multifunctional proteins revealed by overlapping clustering in protein interaction network. \emph{Bioinformatics} \bold{28}, 84-90.

Gambette, P. and Guenoche, A. (2011) Bootstrap clustering for graph partitioning. \emph{RAIRO-Operations Research} \bold{45}, 339-352.

Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.

Shannon, P. et al. (2003) Cytoscape: A software environment for integrated models of biomolecular interaction networks. \emph{Genome Research} \bold{13}, 2498-2504.

Spinelli, L. et al. (2013) Clust&See: a Cytoscape plugin for the identification, visualization, and manipulation of network clusters. \emph{BioSystems} \bold{113}, 91-95.

}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Write a partition file to disk.
\dontrun{
linkcomm2clustnsee(lc)
}

}

