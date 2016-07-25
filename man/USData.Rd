\name{USData}
\docType{data}
\alias{USData}
\title{data: Quarterly US inflation and associated predictors}
\description{
This is the dataset used in Groen et al. (2013) and is downloadable from \url{http://www.tandfonline.com/doi/suppl/10.1080/07350015.2012.727718}.\cr
The variables are:

GDPDEF : Quarterly log changes in the Gross Domestic Product implicit price deflator.\cr
ROUTP  : Real GDP in volume terms.\cr
RCONS  : Real durable personal consumption expenditures in volume terms\cr
RINVR  : Real residential investment in volume terms\cr
PIMP   : Import deflator\cr
UNEMP  : Unemployment ratio\cr
NFPR   : Non-farm payrolls data on employment\cr
HSTS   : Housing starts\cr
OIL    : Real spot price of oil\cr
FOOD   : Real food commodities price index\cr
RAW    : Real raw material commodities price index\cr
M2     : The M2 monetary aggregate\cr
YL     : The level factor\cr
TS     : The slope factor\cr
CS     : Curvature factor\cr
MS     : One-year ahead inflation expectations that come from the Reuters/Michigan Survey of Consumers.

}
\usage{data(USData)}
\format{A \link{xts} object containing 206x16 observations from from 1960-01-01 to 2011-04-01.}
\references{
Groen, Jan JJ, Richard Paap, and Francesco Ravazzolo. (2013) Real--time inflation forecasting in a changing world. \emph{Journal of Business & Economic Statistics}, \bold{31.1} : 29--44.
}
\keyword{datasets}
