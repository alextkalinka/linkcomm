\name{plotLinkCommDend}
\alias{plotLinkCommDend}

\title{
Plot a Coloured Dendrogram of Link Communities
}
\description{
This function is called by \code{plot.linkcomm} to plot a dendrogram of coloured link communities.
}
\usage{
plotLinkCommDend(x, col = TRUE, pal = brewer.pal(9, "Set1"), 
        height = x$pdmax, right = FALSE, labels = FALSE, plotcut = TRUE, 
        droptrivial = TRUE, leaflab = "none", verbose = TRUE, ...)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{col}{
Logical, whether to add community-specific colours. Defaults to TRUE.
}
  \item{pal}{
A character vector describing a colour palette to be used for colouring the communites in the dendrogram plot. Defaults to \code{brewer.pal(9, "Set1")}.
}
  \item{height}{
A numerical value specifying the height at which the dendrogram is cut. Defaults to the maximum partition density height.
}
  \item{right}{
Logical, whether to orient the dendrogram to the right. Defaults to FALSE.
}
  \item{labels}{
Logical, whether to include labels in the dendrogram. Defaults to FALSE.
}
  \item{plotcut}{
Logical, whether to display a horizontal line where the dendrogram is cut. Defaults to TRUE.
}
  \item{droptrivial}{
Logical, whether to not colour communities of size 2. Defaults to TRUE.
}
  \item{leaflab}{
A character string describing the leaf labels on the dendrogram. Can be one of \code{"none"}, \code{"perpendicular"}, or \code{"textlike"}. Defaults to \code{"none"}.
}
  \item{verbose}{
Logical, whether to display the progress of colouring the dendrogram on screen. Defaults to TRUE.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}
\details{
Here we describe the parameters for plotting coloured dendrograms using:\cr
\code{plot(x, type = "dend")}
}
\value{
A dendrogram plot.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{plot.linkcomm}}
}
\examples{
## Generate graph and extract link communities.
g <- swiss[,3:4]
lc <- getLinkCommunities(g)

## Plot a coloured dendrogram.
plot(lc, type = "dend")
}

