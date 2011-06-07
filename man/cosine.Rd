\name{cosine}
\alias{cosine}
\title{Cosine Measure}

\description{Calculates the cosine measure between two vectors or between all
  column vectors of a matrix.}

\usage{cosine(x, y=NULL)}

\arguments{
  \item{x}{A vector or a matrix (e.g., a document-term matrix).}
  \item{y}{Optional: a vector with compatible dimensions to ‘x’. If
    `NULL', all column vectors of ‘x’ are correlated.}
}


\value{Returns a n x n similarity matrix of cosine values, comparing all n
  column vectors against each other. Executed on two vectors, their
  cosine similarity value is returned.}

\seealso{\code{\link[lsa]{cosine}} from which this function was copied.}

\author{Fridolin Wild <email: f.wild@open.ac.uk>}

\examples{

     ## the cosine measure between two vectors

     vec1 = c( 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
     vec2 = c( 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0 )
     cosine(vec1,vec2)




     ## the cosine measure for all vectors of a matrix

     vec3 = c( 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0 )
     matrix = cbind(vec1,vec2, vec3)
     cosine(matrix)

}