\name{Lag}
\alias{Lag}
\title{
	Lag a vector or a matrix of observations
}
\description{
	Lag a vector or a matrix of observations by \code{iK} periods.
}
\usage{
    Lag(mY, iK)
}
\arguments{
\item{mY}{
a vector or a matrix of observations.
}
	\item{iK}{
	An integer indicating the number of lag.
	}
}
\details{
	The function returns a \code{numeric} vector or a \code{matrix} depending on the class of \code{mY}. The dimension of the object is preserved and NA's are introduced for the missing values.
}
\value{
An object of the class \code{numeric} or \code{matrix}.
}
\author{ Leopoldo Catania & Nima Nonejad }
\examples{
\dontrun{
#  Code chunk of Catania and Nonejad (2016) Dynamic Model Averaging
#  for Practitioners in Economics and Finance: The eDMA Package
library(eDMA)

## load data
data("USData")

UDData_lagged <- Lag(USData, 1)

}
}
