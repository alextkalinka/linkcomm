\name{getCommunityConnectedness}
\alias{getCommunityConnectedness}

\title{
Calculate Community Connectedness or Modularity
}
\description{
This function returns a measure of how relatively outwardly or inwardly connected a community is.
}
\usage{
getCommunityConnectedness(x, clusterids = 1:x$numbers[3], conn = "conn", 
                          normalise = TRUE, verbose = FALSE)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{clusterids}{
An integer vector of community IDs. Defaults to all communities.
}
  \item{conn}{
A character string naming the connectedness measure to use. Can be one of \code{"conn"} or \code{"mod"} (see Details below). Defaults to \code{"conn"}.
}
  \item{normalise}{
Logical, whether to normalise community connectedness measures by the number of nodes in individual communities. Defaults to TRUE.
}
  \item{verbose}{
Logical, whether to display the progress of the calculation on the screen. Defaults to FALSE.
}
}
\details{
The connectedness and modularity of different communities indicates whether a particular community is bridging several other communities, or existing as a relatively isolated module. The modularity of community i is
\deqn{M_{i}=\left(\frac{e_{w}(i)}{n_{i}(n_{i}-1)/2}\right).\left(\frac{e_{b}(i)}{n_{i}\hat{d}}\right)^{-1}}{M(i) = (e_{w}(i)/n_{i}(n_{i}-1)/2).(e_{b}(i)/n_{i}*d)^-1}
where \eqn{e_{w}(i)} is the number of edges within community i, \eqn{e_{b}(i)} is the number of edges community i makes to other communities, \eqn{n_{i}} is the number of nodes in community i, and \eqn{\hat{d}} is the average degree in the network. Community connectedness is the inverse of this value.
}
\value{
A named numerical vector, where the names are community IDs and the numbers are community connectedness or modularity scores.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{getCommunityCentrality}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Get community connectedness.
cc <- getCommunityConnectedness(lc, normalise = FALSE)
}

