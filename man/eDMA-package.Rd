\name{eDMA-package}
\alias{eDMA-package}
\alias{eDMA}
\docType{package}
\title{
Dynamic Model Averaging with Modifications
}
\description{
Perform Dynamic Model Averaging with modifications introduced in Dangl and Halling (2012) using parallel computing.
}
\details{
\tabular{ll}{
Package: \tab eDMA\cr
Type: \tab Package\cr
Version: \tab 1.4-0\cr
Date: \tab 2016-11-24\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Leopoldo Catania & Nima Nonejad

Maintainer: Leopoldo Catania <leopoldo.catania@uniroma2.it>
}
\references{
Raftery, Adrian E., Miroslav Karny, and Pavel Ettler. "Online prediction under model uncertainty via dynamic model averaging: Application to a cold rolling mill." Technometrics 52.1 (2010): 52-66.

Dangl, Thomas, and Michael Halling. "Predictive regressions with time-varying coefficients." Journal of Financial Economics 106.1 (2012): 157-181.

Raftery, Adrian E., David Madigan, and Jennifer A. Hoeting. "Bayesian model averaging for linear regression models." Journal of the American Statistical Association 92.437 (1997): 179-191.

Harrison, Jeff, and Mike West. Bayesian Forecasting & Dynamic Models. Springer, 1999.
}
\keyword{ package }
\seealso{ \link{DMA}
}
\examples{
library(eDMA)

## load data
data("USData")

## do DMA, keep the first three predictors fixed and the intercept
Fit <- DMA(GDPDEF ~ Lag(GDPDEF, 1) + Lag(GDPDEF, 2) + Lag(GDPDEF, 3) +
            Lag(ROUTP, 1) + Lag(UNEMP, 1), data = USData, vDelta = c(0.9,0.95,0.99),
             vKeep = c(1, 2, 3, 4))

}
