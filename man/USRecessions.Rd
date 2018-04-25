\name{USRecessions}
\docType{data}
\alias{USRecessions}
\title{data: Dates of U.S. recessions as inferred by GDP-based recession indicator (JHDUSRGDPBR).}
\description{
Dates of U.S. recessions as inferred by GDP-based recession indicator (JHDUSRGDPBR) downloaded from FRED. \url{https://fred.stlouisfed.org/series/JHDUSRGDPBR#0}.

From the FRED website:

The series assigns dates to U.S. recessions based on a mathematical model of the way that recessions differ from expansions. Whereas the NBER business cycle dates are based on a subjective assessment of a variety of indicators, the dates here are entirely mechanical and are calculated solely from historically reported GDP data. Whenever the GDP-based recession indicator index rises above 67\%, the economy is determined to be in a recession. The date that the recession is determined to have begun is the first quarter prior to that date for which the inference from the mathematical model using all data available at that date would have been above 50\%. The next time the GDP-based recession indicator index falls below 33\%, the recession is determined to be over, and the last quarter of the recession is the first quarter for which the inference from the mathematical model using all available data at that date would have been below 50\%.
}
\usage{data(USData)}
\format{A \link{zoo} object containing 206 observations at quarterly frequency from 1960-Q1 to 2011-Q4.}
\examples{
data("USRecessions")
}

\references{
Hamilton, James, Dates of U.S. recessions as inferred by GDP-based recession indicator [JHDUSRGDPBR], retrieved from FRED, Federal Reserve Bank of St. Louis; \url{https://fred.stlouisfed.org/series/JHDUSRGDPBR}, February 4, 2018.
}
\keyword{datasets}
