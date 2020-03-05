\name{getClusterRelatedness}
\alias{getClusterRelatedness}

\title{
Hierarchichal Clustering of Link Communities
}
\description{
This function hierarchically clusters the link communities themselves and returns an object of class \code{hclust}.
}
\usage{
getClusterRelatedness(x, clusterids = 1:x$numbers[3], hcmethod = "ward.D", 
        cluster = TRUE, plot = TRUE, cutat = NULL, col = TRUE, 
        pal = brewer.pal(11, "Spectral"), labels = FALSE, plotcut = TRUE, 
        right = TRUE, verbose = TRUE, ...)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{clusterids}{
An integer vector of community IDs. Defaults to all communities.
}
  \item{hcmethod}{
A character string naming the hierarchical clustering method to use. Can be one of \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"}, or \code{"centroid"}. Defaults to \code{"ward.D"}.
}
  \item{cluster}{
Logical, whether to cluster the communities. If FALSE, the function returns the upper triangular dissimilarity matrix as a vector. Defaults to TRUE.
}
  \item{plot}{
Logical, whether to plot the cluster dendrogram.
}
  \item{cutat}{
A numerical value at which to cut the dendrogram. If NULL, the dendrogram is not cut and meta-communities are not returned. Defaults to NULL.
}
  \item{col}{
Logical, whether to colour the dendrogram. Defaults to TRUE.
}
  \item{pal}{
A character vector describing a colour palette to be used for colouring the meta-communites in the dendrogram plot. Defaults to \code{brewer.pal(11, "Spectral")}.
}
  \item{labels}{
Logical, whether to add labels to the dendrogram plot.
}
  \item{plotcut}{
Logical, whether to display a horizontal line where the dendrogram is cut. Defaults to TRUE.
}
  \item{right}{
Logical, whether to orient the dendrogram to the right. Defaults to TRUE.
}
  \item{verbose}{
Logical, whether to display the progress of the calculation on the screen. Defaults to TRUE.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}
\details{
Extracting meta-communities allows the user to explore community relatedness and structure at higher levels. Community relatedness is calculated using the Jaccard coefficient and the number of nodes that community i and j share:
\deqn{S(i,j)=\frac{|n_{i}\cap n_{j}|}{|n_{i}\cup n_{j}|}}{S(i,j)=|intersect(i,j)|/|union(i,j)|}
}
\value{
Either a numerical vector (the upper triangular dissimilarity matrix - if \code{cluster = FALSE}), a list of integer vectors (the meta-communities - if \code{cutat} is not NULL), or an object of class \code{hclust} (if \code{cluster} is TRUE and \code{cutat} is NULL).
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{meta.communities}}, \code{\link{cutDendrogramAt}}, \code{\link{hclust}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Cluster the link communities.
getClusterRelatedness(lc)

## Cluster the link communities, cut the dendrogram, and return the meta-communities.
getClusterRelatedness(lc, cutat = 1)
}

