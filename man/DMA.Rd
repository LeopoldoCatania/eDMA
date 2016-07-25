\name{DMA}
\alias{DMA}
\title{
	Perform Dynamic Model Averaging
}
\description{
	Implements the Dynamic Model Averaging procedure with the possibility of different valued of the instability parameter.
}
\usage{
DMA(formula, data, vDelta = c(0.9, 0.95, 0.99), dAlpha = 0.99,
    vKeep = NULL, bZellnerPrior = FALSE, dG = 100, bParallelize = TRUE, iCores = NULL)
}
\arguments{
\item{formula}{
an object of class \link{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted
}
	\item{data}{
	an object of the class \link{data.frame}, (or object coercible by as.data.frame to a data frame) containing the variables in the model. It can also be an object of the classes \code{ts}, \code{xts} or \code{zoo}. If this is the case, the time information is used in the graphical representation of the results.
	}
	\item{vDelta}{
   	D x 1 numeric vector representing the \bold{\eqn{\delta}} parameter. By default \code{vDelta = c(0.9, 0.95, 0.99)}.
	}
	\item{dAlpha}{
		numeric variable representing \eqn{\alpha}. By default \code{dAlpha = 0.99}.
	}
	\item{vKeep}{
		\code{numeric} vector of indexes representing the predictors that must be always included in the models. The combinations of predictors that do not include the variables declared in \code{vKeep} are automatically discarded. The indexes must be consistent with the model description given in \code{formula}, i.e., if the first and the fourth variables always have to be included, then we must set vKeep=c(1, 4). Note that, the intercept (if not removed from \code{formula}) is always in the first position. It can also be a character vector indicating the names of the predictors if these are consistent with the provided \code{formula}. If \code{vKeep = "KS"} the "Kitchen Sink" formulation is adopted, i.e., all the predictors are always included, see, e.g., Paye (2012).	By default all the combinations are considered, i.e. \code{vKeep = NULL}.
	}
		\item{bZellnerPrior}{
	Boolean variable indicating whether the Zellner prior should be used for the coefficients at time t=0. Default = \code{FALSE}.
	}
		\item{dG}{
		numeric variable equal to 100 by default. If \code{bZellnerPrior = TRUE} this represent the variable 'g' in Eq. (4) of Dangl Halling (2012). Otherwise, if \code{bZellnerPrior = FALSE} it represents the scaling factor for the variance covariance matrix of the normal prior for \bold{\eqn{\theta_0}}, i.e. \bold{\eqn{\theta_0}}~N(0,dG*I) where I is the identity matrix.
	}
		\item{bParallelize}{
		 Boolean variable indicating whether to use multiple processors to speed up the computations. By default \code{bParallelize = TRUE}.
	}
		\item{iCores}{
	integer indicating the number of cores to use if \code{bParallelize = TRUE}. By default all but one cores are used. The number of cores is guessed using the
        \code{detectCores()} function from the \code{parallel} package
	}
}
\details{
See Catania and Nonejad (2016) for further details.
}
\value{
An object of the class \code{DMA}, see \link{DMA-class}.
}
\references{
Dangl, T., & Halling, M. (2012). Predictive regressions with time--varying coefficients. \emph{Journal of Financial Economics}, \bold{106(1)}, 157--181.  \url{doi:10.1016/j.jfineco.2012.04.003}.\cr

Catania, Leopoldo, and Nima Nonejad. "Dynamic Model Averaging for Practitioners in Economics and Finance: The eDMA Package." arXiv preprint arXiv:1606.05656 (2016).\cr

Paye, B.S. (2012). 'Deja vol': Predictive Regressions for Aggregate Stock Market Volatility
Using Macroeconomic Variables.\emph{Journal of Financial Economics}, \bold{106(3)}, 527-546. ISSN 0304-405X.
\url{doi:10.1016/j.jfineco.2012.06.005}. URL \url{http://www.sciencedirect.com/science/article/pii/S0304405X12001316}.
}
\author{ Leopoldo Catania & Nima Nonejad }
\examples{
\dontrun{
#  Code chunk of Catania and Nonejad (2016) Fast Dynamic Model Averaging
#  for Practitioners in Economics and Finance: The eDMA Package
library(eDMA)

## load data
data("USData")

## do DMA, keep the first three predictors fixed and the intercept
Fit = DMA(GDPDEF ~ Lag(GDPDEF, 1) + Lag(GDPDEF, 2) + Lag(GDPDEF, 3) +
            Lag(ROUTP, 1) + Lag(UNEMP, 1), data = USData, vDelta = c(0.9,0.95,0.99),
             vKeep = c(1, 2, 3, 4))

Fit
}
}
