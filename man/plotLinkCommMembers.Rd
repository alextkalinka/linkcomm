\name{plotLinkCommMembers}
\alias{plotLinkCommMembers}

\title{
Plot a Community Membership Matrix for Link Communities
}
\description{
This function is called by \code{plot.linkcomm} to plot a community membership matrix for the link communities.
}
\usage{
plotLinkCommMembers(x, nodes = head(names(x$numclusters), 10), 
        pal = brewer.pal(11, "Spectral"), shape = "rect", total = TRUE, 
        fontsize = 11, nspace = 3.5, maxclusters = 20)
}

\arguments{
  \item{x}{
An object of class \code{linkcomm}.
}
  \item{nodes}{
A character vector specifying the node names that will be included in the plot. Defaults to the 10 nodes that belong to the most communities.
}
  \item{pal}{
A character vector describing a colour palette to be used for community-specific colouring. Defaults to \code{brewer.pal(11, "Spectral")}.
}
  \item{shape}{
A character string specifying the shape of matrix entries. Can be one of \code{"rect"} or \code{"circle"}. Defaults to \code{"rect"}.
}
  \item{total}{
Logical, whether to display the number of communities each node belongs to and the number of nodes in each community. Defaults to TRUE.
}
  \item{fontsize}{
A numerical value specifying font size for the node names. Defaults to 11.
}
  \item{nspace}{
A numerical value specifying how much space to leave at the left for fitting in node names. Defaults to 3.5.
}
  \item{maxclusters}{
An integer value specifying the maximum number of communities to display. Defaults to 20.
}
}

\value{
A community membership matrix plot.
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

## Plot a community membership matrix.
plot(lc, type = "members")
}

