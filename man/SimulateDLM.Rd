\name{SimulateDLM}
\alias{SimulateDLM}
\title{
	Simulate from DLM of West and Harrison (1999).
}
\description{
		Simulate from DLM of West and Harrison (1999), as in Section 2 of Catania and Nonejad (2016).
}
\usage{
    SimulateDLM(iT, mX, vBeta0, mW, dV, dPhi)
}
\arguments{
\item{iT}{
\code{numeric}, number of observation to simulate
}
\item{mX}{
\code{matrix} of dimension \code{iT} x N m where N is the number of covariates.
}
\item{vBeta0}{
\code{numeric} vector with initial value for the regressor coefficients.
}
\item{mW}{
\code{matrix} covariance matrix of the state equation.
}
\item{dV}{
\code{numeric} variance of the observation (measurement equation).
}
\item{dPhi}{
\code{numeric} value for the autoregressive parameter of the regressors. It imposes that all the regressors have the same autoregressive parameters, if \code{dPhi = 1}, then the regressors evolve as random-walks.
}


}
\details{
	The function returns a \code{list} of two elements: \code{vY} and \code{mBeta}. \code{vY} is a \code{iT} x 1 \code{numeric} vector of simulated dependent variables. \code{mBeta} is a \code{matrix} of dimension \code{iT x ncol(mX)} of regressor coefficients.
}
\value{
An object of the class \code{list}.
}
\author{ Leopoldo Catania & Nima Nonejad }
\references{
Catania, Leopoldo, and Nima Nonejad. "Dynamic Model Averaging for Practitioners in Economics and Finance: The eDMA Package." arXiv preprint arXiv:1606.05656 (2016).

West, Mike. Bayesian forecasting. John Wiley & Sons, Inc., 1999.
}
\examples{
\dontrun{


set.seed(7892)

iT = 500
iK = 3

dV = 0.1
mW = diag(iK + 1) * 0.01
dPhi = 1

vBeta0 = rep(0, iK + 1)
mX = cbind(1, matrix(rnorm(iT * (iK)), iT, iK))

lOut = SimulateDLM(iT, mX, vBeta0, mW, dV, dPhi)
vY = lOut$vY

}
}
