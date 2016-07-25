\name{SimData}
\docType{data}
\alias{SimData}
\title{data: Simulated data from DLM of West and Harrison (1999).}
\description{
This is the simulated dataset used in Section 4.1 of Catania and Nonejad (2016).
}
\usage{data(USData)}
\format{A \link{matrix} object containing 500 x 6 simulated observations.}
\examples{
\dontrun{

#the data set has been generated as:

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

mX = mX[, -1]

iK_Add = 2

mX_add = matrix(rnorm(iT * iK_Add), iT, iK_Add)

SimData = cbind(y = vY, mX, mX_add)
colnames(SimData) = c("y", paste("x", 2:(iK + iK_Add + 1), sep = ""))

}
}

\references{
Catania, Leopoldo, and Nima Nonejad. "Dynamic Model Averaging for Practitioners in Economics and Finance: The eDMA Package." arXiv preprint arXiv:1606.05656 (2016).\cr
West, Mike. Bayesian forecasting. John Wiley & Sons, Inc., 1999.
}
\keyword{datasets}
