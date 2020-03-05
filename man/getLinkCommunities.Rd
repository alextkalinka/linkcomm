\name{getLinkCommunities}
\alias{getLinkCommunities}

\title{
Extract Link Communities from a Network
}
\description{
This function extracts link communities from networks of arbitrary size and type.
}
\usage{
getLinkCommunities(network, hcmethod = "average", use.all.edges = FALSE,
                   edglim = 10^4, directed = FALSE, dirweight = 0.5,
                   bipartite = FALSE, dist = NULL, plot = TRUE, 
                   check.duplicates = TRUE, removetrivial = TRUE, 
                   verbose = TRUE)
}

\arguments{
  \item{network}{
An edge list, which is a matrix or data frame with 2 or 3 columns. The first 2 columns contain the nodes that interact with each other, which can be character strings or integer values. The optional third column is a numerical vector of weights for each edge. Can also be a character string naming a file containing an edge list.
}
  \item{hcmethod}{
A character string naming the hierarchical clustering method to use. Can be one of \code{"ward"}, \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"}, or \code{"centroid"}. Defaults to \code{"average"} (if the number of edges is greater than \code{edglim} then \code{"single"} is used).
}
  \item{use.all.edges}{
Logical, indicating whether edge similarities should be calculated for all pairs of edges (TRUE), or only for edge pairs that share a node (FALSE) as in the original Ahn et al. (2010) algorithm. Defaults to FALSE. If TRUE, networks are treated as undirected.
}
  \item{edglim}{
An integer value indicating the largest number of edges permissible for the hierarchical clustering to be handled in memory. Above this value the upper triangular dissimilarity matrix will be written to disk and read and written as clustering proceeds until the file size is 0 bytes (see Details below). Defaults to \eqn{10^{4}}.
}
  \item{directed}{
Logical, whether the network is directed. Defaults to FALSE.
}
  \item{dirweight}{
A numerical value between 1 and 0 inclusive indicating the weight that will be attached to edges that share a node but are in the opposite orientation. Defaults to 0.5. Will be ignored if \code{directed} = FALSE.
}
  \item{bipartite}{
Logical, whether the input network is bi-partite. See Details for an explanation of how bi-partite networks are handled. Defaults to FALSE.
}
  \item{dist}{
An object of class \code{"dist"} representing a user-defined distance matrix for the network. Note, this must be the lower triangular matrix of an n*n distance matrix, where n is the number of edges in the network (make sure duplicated edges are removed). If NULL, then the distance matrix is calculated by the algorithm. Defaults to NULL.
}
  \item{plot}{
Logical, whether to plot summary output from the algorithm (dendrogram and partition density plot). Defaults to TRUE. Note, if there are more than 1500 but less than edglim edges then the dendrogram will be plotted without colour and in a separate panel from the partition density to avoid lengthy rendering times; when there are more than edglim edges then only the partition density will be plotted.
}
  \item{check.duplicates}{
Logical, whether to check for and remove loops, duplicate edges, and bi-directional edges. Defaults to TRUE. Note, if you wish to avoid this step by setting this parameter to FALSE then you must be certain that there are no duplicate edges in the network.
}
  \item{removetrivial}{
Logical, whether to remove trivial community clusters that contain 2 edges. Defaults to TRUE.
}
  \item{verbose}{
Logical, whether to display the progress of the algorithm on the screen. Defaults to TRUE.
}
}
\details{
This is the main algorithm used for extracting link communities from networks of arbitrary size and type. Input networks may be directed, weighted, both directed and weighted, or neither. The algorithm used is the one outlined by Ahn et al. (2010). The similarity between links, \eqn{e_{ik}} and \eqn{e_{jk}}, that share a node, \eqn{k}{k}, is calculated using the Jaccard coefficient
\deqn{S(e_{ik},e_{jk})=\frac{|n_{+}(i)\cap n_{+}(j)|}{|n_{+}(i)\cup n_{+}(j)|}}{S(e_{ik},e_{jk}) = |intersect(n_{+}(i), n_{+}(j))|/union(|n_{+}(i), n_{+}(j))|}
where \eqn{n_{+}(i)} refers to the first-order node neighbourhood of node \eqn{i}{i}, which includes node \eqn{i}{i} itself (inclusive neighbour set). After assigning pairwise similarities to all of the links in the network, the links are hierarchically clustered using single-linkage clustering, and the resulting dendrogram is cut at a point that maximises the density of links within the clusters normalising against the maximum and minimum numbers of links possible in each cluster, known as the partition density.
For directed and weighted networks, the Tanimoto coefficient is used for assigning similarity between links
\deqn{S(e_{ik},e_{jk})=\frac{\mathbf{a}_{i}.\mathbf{a}_{j}}{|\mathbf{a}_{i}|^{2}+|\mathbf{a}_{j}|^{2}-\mathbf{a}_{i}.\mathbf{a}_{j}}}{S(e_{ik},e_{jk}) = (a_{i}.a_{j})/(|a_{i}|^{2}+|a_{j}|^{2}-a_{i}.a_{j})}
where \eqn{\mathbf{a}_{i}}{a_{i}} refers to a vector describing the weights of links between node \eqn{i}{i} and the nodes in the first-order neighbourhoods of both nodes \eqn{i}{i} and \eqn{j}{j} (equal to 0 in the event of an absent link). For directed networks, links to nodes shared by both node \eqn{i}{i} and \eqn{j}{j} are given a user-defined weight below 1 if they are in the opposite orientation.\cr\cr\cr
For bi-partite networks, the set of neighbours (instead of the inclusive neighbour set) is used to count nodes for the edge similarity metric because node \eqn{i}{i} and node \eqn{j}{j} cannot share an edge in a bi-partite network. The partition density for bi-partite networks is calculated as:
\deqn{D_{c} = \frac{2}{M}\sum_{c}m_{c}\frac{m_{c}+1-n_{c}}{2n_{c0}n_{c1}-2(n_{c}-1)}}{ (2/M)*sum_c m_c*(m_c+1-n_c)/(2*n_c0*n_c1-2*(n_c-1)) }
where \eqn{M}{M} is the total number of edges, \eqn{m_c}{m_c} is the number of edges in subset \eqn{c}{c}, \eqn{n_c}{n_c} is the number of nodes in subset \eqn{c}{c}, \eqn{n_{c0}}{n_c0} is the number of nodes in partition 0, and \eqn{n_{c1}}{n_c1} is the number of nodes in partition 1.

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
Ahn, Y.Y., Bagrow, J.P., and Lehmann, S. (2010). Link communities reveal multiscale complexity in networks. \emph{Nature} \bold{466}, 761-764.

Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}
\note{
When the number of links is less than \code{edglim} the hierarchical clustering will be handled in memory. Above this value the upper triangular dissimilarity matrix will be compressed and written to disk and read and written as clustering proceeds until the file size is 0 bytes using a compiled C++ function. In this case the hierarchical clustering method will always be \code{"single"} to enhance performance for large networks. The size of \code{edglim} can be modified to suit the computer resources available to the user. As a guide, a network with \eqn{10^{4}} links will require \eqn{((10^{4})^{2})*8 = 800} MB to be handled in an uncompressed format in the memory.\cr\cr\cr
For directed networks, a pair of bidirectional interactions between two nodes cannot be assigned similarities and the edge that appears lower in the edge list for the network will be discarded.\cr\cr\cr
When \code{use.all.edges} is TRUE, the algorithm may be slow as all pairs of edges will be compared (\eqn{n^2}{n^2} comparisons, where \eqn{n}{n} is the number of edges).
}

\seealso{
\code{\link{plot.linkcomm}}, \code{\link{newLinkCommsAt}}, \code{\link{meta.communities}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Extract communities by writing a temporary file to disk.
lc <- getLinkCommunities(g, edglim = 10)

## Use similarities between all pairs of edges.
lc <- getLinkCommunities(g, use.all.edges = TRUE)

## Directed network.
lc <- getLinkCommunities(g, directed = TRUE, dirweight = 0.8)

## Weighted network.
g <- cbind(swiss[,3:4], runif(nrow(swiss[,3:4])))
lc <- getLinkCommunities(g)

## Directed and weighted network.
lc <- getLinkCommunities(g, directed = TRUE, dirweight = 0.8)

}

