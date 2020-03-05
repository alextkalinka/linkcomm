\name{edge.duplicates}
\alias{edge.duplicates}

\title{
Find and Remove Edge Loops and Duplicates
}
\description{
This function finds and removes loops, edge duplicates, and bi-directional edges.
}
\usage{
edge.duplicates(network, verbose = TRUE)
}

\arguments{
  \item{network}{
An edge list, which is either a data frame or a character or integer matrix with two columns.
}
  \item{verbose}{
Logical, whether to display the progress of the function on the screen. Defaults to TRUE.
}
}
\details{
This function removes loops, duplicate and bi-directional edges; the edges that occur closer to the end of the edge list will be removed.
}
\value{
A list with the following elements: \code{edges} - a character matrix of the edges in the network with any loops or duplicate edges removed; \code{inds} - an integer vector of the edge indices of any loop or duplicate edges in the original network.
}
\references{
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. \emph{Bioinformatics} \bold{27}, 2011-2012.
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\examples{
## Generate graph and remove loops and duplicate edges.
g <- swiss[,3:4]
g.dup <- edge.duplicates(g)

}

