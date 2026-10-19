\name{PowerSet}
\alias{PowerSet}
\title{
	Build the power set of the values \{0,1,2,...,\code{iK-1}\}.
}
\description{
	Build the power set of the values \{0,1,2,...,\code{iK-1}\}.
}
\usage{
    PowerSet(iK)
}
\arguments{
\item{iK}{
  \code{numeric} an integer between 0 and 20 giving the number of elements in the series \{0,1,2,...,\code{iK-1}\}.
}
}
\details{
	The function returns a \code{list} of \code{numeric} vectors with the indices representing all the 2^iK subsets.
	The empty subset \{\} is represented by the \code{numeric(0)} vector.
	}
\value{
An object of the class \code{list}.
}
\author{ Leopoldo Catania }
\examples{

library(eDMA)
PowerSet(5)

}
