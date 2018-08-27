\name{PowerSet}
\alias{PowerSet}
\title{
	Build the power set of the values \{0,1,2,...,\code{iK}\}.
}
\description{
	Build the power set of the values \{0,1,2,...,\code{iK}\}.
}
\usage{
    PowerSet(iK)
}
\arguments{
\item{iK}{
  \code{numeric} an integer value indicating the end of the series \{0,1,2,...,\code{iK}\}.
}
}
\details{
	The function returns a \code{list} of \code{numeric} vectors with the indices representing all the 2^iK subsets.
	The empty subset \{\} is represented by the \code{numeric(0)} vector.
	}
\value{, 4
An object of the class \code{list}.
}
\author{ Leopoldo Catania & Nima Nonejad }
\examples{

library(eDMA)
PowerSet(5)

}
