\name{meta.communities}
\alias{meta.communities}

\title{
Produce a Set of Meta-Communities
}
\description{
This function returns meta-communities by hierarchically clustering link communities based on their number of shared nodes.
}
\usage{
meta.communities(x, hcmethod = "ward.D", deepSplit = FALSE)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm} or \code{OCG}.
}
  \item{hcmethod}{
A character string naming the hierarchical clustering method to use. Can be one of \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"}, or \code{"centroid"}. Defaults to \code{"ward.D"}.
}
  \item{deepSplit}{
Logical or integer value from 0 to 4 indicating how fine-grained the clusters should be with larger values giving increased cluster structure. Defaults to FALSE.
}
}
\details{
Extracting meta-communities allows the user to explore community relatedness and structure at higher levels. Community relatedness is calculated using the Jaccard coefficient and the number of nodes that community i and j share:
\deqn{S(i,j)=\frac{|n_{i}\cap n_{j}|}{|n_{i}\cup n_{j}|}}{S(i,j)=|intersect(i,j)|/|union(i,j)|}
\cr
A hierarchical clustering dendrogram is generated based on the above distance metric and then an algorithm is used to automatically extract clusters (Langfelder et al. 2008). This function works best on large networks with a large number of link communities.
}

\value{
An object of class \code{linkcomm}, which is a list containing the following components:
\item{numbers }{An integer vector with the number of edges, nodes, and communities.}
\item{hclust }{An object of class \code{\link{hclust}}, which contains information about the hierarchical clustering of links.}
\item{pdmax }{A numerical value indicating the height of the dendrogram at which the partition density is maximised.}
\item{pdens }{A numerical matrix with 2 columns; the first is the heights at which clusters appear and the second is the partition density.}
\item{nodeclusters }{A data frame consisting of 2 columns; the first contains node names, and the second contains single community IDs for each node. All communities and their nodes are represented, but not necessarily all nodes.}
\item{clusters }{A list of integer vectors containing the link IDs that belong to each community. Community IDs are the numerical position of the communities in the list.}
\item{edges }{A data frame with 3 columns; the first two contain nodes that interact with each other, and the third is an integer vector of community IDs indicating community membership for each link.}
\item{numclusters }{A named integer vector. Names are node names and integer values are the number of communities to which each node belongs.}
\item{clustsizes }{A named integer vector. Names are community IDs and integer values indicate the number of nodes that belong in each community.}
\item{igraph }{An object of class \code{\link{igraph}}. The network is represented here as an \code{igraph} object.}
\item{edgelist }{A character matrix with 2 columns containing the nodes that interact with each other.}
\item{directed }{Logical indicating whether the network is directed.}
\item{bipartite }{Logical indicating whether the network is bi-partite.}
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.

Langfelder, P., Zhang, B., and Horvath, S. (2008) Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. \emph{Bioinformatics} \bold{24}, 719-720.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}
\seealso{
\code{\link{getClusterRelatedness}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Produce meta-communities.
\dontrun{mc <- meta.communities(lc)}
}

