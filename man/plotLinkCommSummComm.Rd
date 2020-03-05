\name{plotLinkCommSummComm}
\alias{plotLinkCommSummComm}

\title{
Plot a Summary of the Link Communities
}
\description{
This function is called by \code{plot.linkcomm} to plot either connectedness or modularity of individual link communities.
}
\usage{
plotLinkCommSummComm(x, clusterids = 1:x$numbers[3], summary = "conn", 
        pie = FALSE, col = TRUE, pal = brewer.pal(11, "Spectral"), 
        random = FALSE, verbose = TRUE, ...)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{clusterids}{
An integer vector of community IDs. Defaults to all communities.
}
  \item{summary}{
A character string specifying the community summary. Can be one of \code{"conn"}, \code{"mod"}, \code{"ld"} for connectedness, modularity, and link densities respectively. Defaults to \code{"conn"}.
}
  \item{pie}{
Logical, whether to plot a pie graph. If FALSE, a bar plot is plotted. Defaults to FALSE.
}
  \item{col}{
Logical, whether to colour each community differently. Defaults to TRUE.
}
  \item{pal}{
A character vector describing a colour palette to be used for colouring the link communities. Defaults to \code{brewer.pal(11, "Spectral")}.
}
  \item{random}{
Logical, whether to randomise the link colours. Defaults to FALSE.
}
  \item{verbose}{
Logical, whether to print the progress of the calculation to the screen. Defaults to TRUE.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}
\details{
Here we describe the parameters for plotting link community summaries using:\cr
\code{plot(x, type = "commsumm", type = "mod")}
}
\value{
A bar graph or pie chart.
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

\dontshow{dev.off()}
## Plot the modularity of the link communities.
plot(lc, type = "commsumm", summary = "mod")
}

