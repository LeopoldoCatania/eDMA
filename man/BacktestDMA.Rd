\name{BacktestDMA}
\alias{BacktestDMA}
\title{
	Backtest measures for Dynamic Model Averaging and comparison with Dynamic Model Selection
}
\description{
	Backtest measures for Dynamic Model Averaging and comparison with Dynamic Model Selection. This function evaluates the out of sample performance of DMA and compare it with DMS.
}
\usage{
    BacktestDMA(object, iBurnPeriod = NULL)
}
\arguments{
\item{object}{
an object of the class \link{DMA-class}, created using the function \link{DMA}.
}
	\item{iBurnPeriod}{
	An integer indicating the length of the burn-in period. By default \code{iBurnPeriod = NULL}. If \code{iBurnPeriod = NULL} then one third of the total sample is used as the burn-in in period and a warning is printed.
	}
}
\details{
	The function returns a \code{matrix} with Mean Squared Error (MSE), Mean Absolute Error (MAD) and Predictive Likelihood for DMA and DMS using the predictions during the out-of-sample period.
}
\value{
An object of the class \code{matrix}.
}
\author{ Leopoldo Catania & Nima Nonejad }
\examples{

library(eDMA)

## load data
data("USData")

## do DMA, keep the first three predictors fixed and the intercept
Fit <- DMA(GDPDEF ~ Lag(GDPDEF, 1) + Lag(GDPDEF, 2) + Lag(GDPDEF, 3) +
            Lag(ROUTP, 1) + Lag(UNEMP, 1), data = USData, vDelta = c(0.9, 0.95, 0.99),
            vKeep = c(1, 2, 3))

BacktestDMA(Fit, iBurnPeriod = 32)

}
