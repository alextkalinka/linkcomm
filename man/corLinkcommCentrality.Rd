\name{corLinkcommCentrality}
\alias{corLinkcommCentrality}

\title{
Correlation of Community Centrality with Classic Centrality
}
\description{
This function calculates the correlation between the community centrality and classic centrality measures for a set of nodes in a network, and plots a scatterplot of the relationship together with a fitted straight line.
}
\usage{
corLinkcommCentrality(x, centrality = "degree", type = "commweight", 
                      method = "spearman", plot = TRUE, pch = 20, ...)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{centrality}{
A character string naming the classic centrality measure. Can be one of \code{"degree"}, \code{"betweenness"}, \code{"closeness"}, and \code{"constraint"}. Defaults to \code{"degree"}.
}
  \item{type}{
A character string naming the type of community centrality. Can be \code{"commweight"} or \code{"commconn"}, defaults to \code{"commweight"}.
}
  \item{method}{
A character string naming the correlation method. Can be one of \code{"spearman"}, \code{"pearson"}, or \code{"kendall"}. Defaults to \code{"spearman"}.
}
  \item{plot}{
Logical, whether to plot a scatterplot of the relationship, defaults to TRUE.
}
  \item{pch}{
An integer specifying the plot symbol (see \code{\link{par}}). Defaults to 20.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}
\details{
The correlation between community centrality and classic centrality measures, such as degree or betweenness, may reveal discrepancies, thereby indicating that community centrality scores provide a unique reflection of node importance.
}
\value{
A correlation coefficient.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getCommunityCentrality}} 
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Correlate community centrality with degree centrality.
corLinkcommCentrality(lc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% __ONLY ONE__ keyword per line
