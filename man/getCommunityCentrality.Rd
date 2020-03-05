\name{getCommunityCentrality}
\alias{getCommunityCentrality}

\title{
Calculate Community Centrality Measures for Nodes
}
\description{
This function returns community-based node centrality measures.
}
\usage{
getCommunityCentrality(x, nodes = names(x$numclusters), type = "commweight", 
                       normalise = TRUE)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{nodes}{
A character vector giving the names of nodes for calculating community centrality scores. Defaults to all nodes.
}
  \item{type}{
A character string naming the community centrality measure. Can be one of \code{"commweight"} or \code{"commconn"} (see Details below). Defaults to \code{"commweight"}.
}
  \item{normalise}{
Logical, whether to normalise community connectedness for \code{"commconn"}. Defaults to TRUE. Will be ignored for \code{"commweight"}.
}
}
\details{
Community-based measures of node centrality provide an alternative to classic measures of node centrality. \code{"commweight"} weights each community that a node belongs to by how similar that community is to each of the other communities to which the node also belongs. For node i the community centrality is
\deqn{C_{C}(i)=\sum_{i\in j}^{N}\left(1-\frac{1}{m}\sum_{i\in j\cap k}^{m}S(j,k)\right)}{Cc(i) = sum_{i in j}^{N}(1 - (1/m)*sum_{i in j&k}^{m}(S(j,k)))}
where the main sum is over the \eqn{N} communities to which node \eqn{i} belongs, and \eqn{S(j,k)} refers to the similarity between community \eqn{j} and \eqn{k}, calculated as the Jaccard coefficient for the number of shared nodes between each community pair, and this is averaged over the \eqn{m} communities paired with community \eqn{j} and in which node \eqn{i} jointly belongs.
\code{"commconn"} weights each community that a node belongs to by how many connections the community forms outside of itself relative to how many connections the community has within itself (the inverse of modularity), so that nodes that belong to more highly connecting communitites will receive a higher community centrality score. For node i the community centrality is
\deqn{C_{C}(i)=\sum_{i\in j}^{N}e_{ij}\frac{\check{e}_{B(j)}}{\check{e}_{W(j)}}}{Cc(i) = sum_{i in j}^{N}e_{ij}.(\check{e}_{B(j)}/\check{e}_{W(j)})}
where \eqn{e_{ij}} is the number of edges node i has in community j, \eqn{\check{e}_{B(j)}=\frac{e_{B(j)}}{n_{j}\bar{d}}} is the number of edges community j makes outside of itself normalised by the number of nodes in community j multiplied by the average degree in the network, and \eqn{\check{e}_{W(j)}=\frac{e_{W(j)}}{n(n-1)/2}} is the number of edges within community j normalised by the total number possible.
}
\value{
A named numerical vector where the names are node names and the numbers are community centrality measures.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{getCommunityConnectedness}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Calculate community centrality.
cc <- getCommunityCentrality(lc)

## Calculate community centrality using "commconn" measure.
cc <- getCommunityCentrality(lc, type = "commconn")
}

