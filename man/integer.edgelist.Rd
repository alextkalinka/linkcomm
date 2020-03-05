\name{integer.edgelist}
\alias{integer.edgelist}

\title{
Convert A Network to an Integer Edgelist
}
\description{
This function converts a character string edgelist into an integer edgelist.
}
\usage{
integer.edgelist(network)
}

\arguments{
  \item{network}{
An edge list, which is a matrix or data frame with 2 or 3 columns. The first 2 columns contain the nodes that interact with each other, which can be character strings or integer values. The optional third column is a numerical vector of weights for each edge, which is stripped from the output.
}

}

\value{
A list containing the following components:
\item{edges }{A matrix with two columns containing the integer edgelist.}
\item{nodes }{A named integer vector mapping node integer IDs to their character string equivalents.}

}

\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}


\seealso{
\code{\link{getLinkCommunities}} 
}
\examples{
## Generate graph and convert into an integer edgelist.
g <- data.frame(letters[1:10],letters[6:15])
gi <- integer.edgelist(g)
}

