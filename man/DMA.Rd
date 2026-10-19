\name{DMA}
\alias{DMA}
\title{
	Perform Dynamic Model Averaging
}
\description{
	Implements the Dynamic Model Averaging procedure with the possibility of also performing averaging over a grid of forgetting factor values
}
\usage{
DMA(formula, data = NULL, vDelta = c(0.9, 0.95, 0.99), dAlpha = 0.99,
    vKeep = NULL, bZellnerPrior = FALSE, dG = 100, bParallelize = TRUE,
    iCores = NULL, dBeta = 1.0, dS0 = 1.0, dN0 = 2.0)
}
\arguments{
\item{formula}{
an object of class \link{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
}
	\item{data}{
	an object of the class \link{data.frame}, (or object coercible by as.data.frame to a data frame) containing the variables in the model. It can also be an object of the classes \code{ts}, \code{xts} or \code{zoo}. If this is the case, the time information is used in the graphical representation of the results. The last observation of the dependent variable can be equal to \code{NA}. This is the case when a series of length \eqn{T} is available but the user wants to have predictions for time \eqn{T+1}, see Details.
	}
	\item{vDelta}{
   	D x 1 numeric vector representing the the grid of values of \bold{\eqn{\delta}}. By default \code{vDelta = c(0.9, 0.95, 0.99)}.
	}
	\item{dAlpha}{
		numeric variable representing \eqn{\alpha}. By default \code{dAlpha = 0.99}.
	}
	\item{vKeep}{
		\code{numeric} vector of indices representing the predictors that must be always included in the models. The combinations of predictors that do not include the variables declared in \code{vKeep} are automatically discarded. The indices must be consistent with the model description given in \code{formula}, i.e., if the first and the fourth variables always have to be included, then we must set vKeep=c(1, 4). Note that, the intercept (if not removed from \code{formula}) is always in the first position. It can also be a character vector indicating the names of the predictors if these are consistent with the provided \code{formula}. If \code{vKeep = "KS"} the "Kitchen Sink" formulation is adopted, i.e., all the predictors are always included, see, e.g., Paye (2012).	By default all the combinations are considered, i.e. \code{vKeep = NULL}.
	}
		\item{bZellnerPrior}{
	Boolean variable indicating whether the Zellner prior should be used on the coefficients at time t=0. Default = \code{FALSE}.
	}
		\item{dG}{
		numeric variable equal to 100 by default. If \code{bZellnerPrior = TRUE} this represent the variable 'g' in Eq. (4) of Dangl Halling (2012). Otherwise, if \code{bZellnerPrior = FALSE} it represents the scaling factor for the variance covariance matrix of the normal prior for \bold{\eqn{\theta_0}}, i.e. \bold{\eqn{\theta_0}}~N(0,dG*I) where I is the identity matrix.
	}
		\item{bParallelize}{
		 Boolean variable indicating whether to use multiple processors to speed up the computations. By default \code{bParallelize = TRUE}.
	}
		\item{iCores}{
	integer indicating the number of cores to use if \code{bParallelize = TRUE}. By default all but one cores are used. The number of cores is guessed using the
        \code{detectCores()} function from the \code{parallel} package.
	}
			\item{dBeta}{
	numeric scalar in (0, 1] indicating the forgetting factor for the measurement variance, see Prado and West (2010) pp. 132 and Beckmann and Schussler (2014).
	}
\item{dS0}{Positive initial observation scale, fixed before observing the data. Defaults to 1.}
\item{dN0}{Positive initial degrees-of-freedom accumulator. Defaults to 2.}

}
\details{
All discounts must lie in (0, 1]. The model space consists of nonempty
subsets of the columns of the model matrix, including factor expansions and
formula transformations. Indices in \code{vKeep} refer to these columns;
\code{"KS"} includes all columns. At most 20 optional columns are allowed.
An estimated working-memory limit of 1 GiB is checked before allocating
model/delta state arrays. Formula offsets are not supported.

At each step both the conditional model probabilities and delta probabilities
are raised to \code{dAlpha} and separately normalized. These prediction
probabilities are used consistently for point forecasts, predictive densities,
and all predictive variance components. After observing the response, Bayes'
rule updates both probability layers using the same conditional predictive
likelihoods. Filtered coefficients use these posterior probabilities, with
absent coefficients represented by zeros. Normalization is performed in log
space. This coherent discounted prediction convention is an explicit change
from the lagged-weight forecast formulas in the published paper; see NEWS.

The coefficient prior has mean zero and covariance \code{dG * I} unless
\code{bZellnerPrior = TRUE}. The first row is processed through the same
prediction and update equations as subsequent rows, with no coefficient
covariance inflation before that first observation. Observation scale starts
at \code{dS0} and the degrees-of-freedom accumulator starts at \code{dN0}.
Predictive Student-t densities retain the package's convention
\eqn{\nu_t = \beta n_{t-1}+1}. The scale is
\eqn{Q_t = S_{t-1}+x_t'R_t x_t}. Observation and coefficient contributions
to predictive variance multiply their scale contributions by
\eqn{\nu_t/(\nu_t-2)} when \eqn{\nu_t>2}; otherwise the predictive
variance is infinite. Between-model and between-delta contributions use
prediction means and prediction probabilities. The coefficient covariance
update remains the plug-in Kalman recursion, without a subsequent scale-ratio
rescaling; it is not an exact Normal--inverse-gamma filter.

The Zellner covariance uses least squares on all supplied observed responses,
with residual sum of squares divided by \eqn{T-1}, and is reused across deltas.
A terminal forecast-only row is excluded. A residual-scale floor of
\code{dS0 * 1e-12} prevents a zero prior covariance for an exact fit.
Rank-deficient or ill-conditioned designs are rejected. Because this prior
uses the complete observed sample, historical forecasts with the Zellner
option are not real-time forecasts. Use the fixed prior for causal filtering.

DMS forecasts retain the nested selection rule: select the largest prediction
probability for delta, then the largest conditional model probability within
that delta. The model-size and maximum/top-ten-percent probability diagnostics
instead describe the marginal posterior distribution over models.

Rows with incomplete data are omitted, except that a final missing response
with complete predictors is retained for prediction only. Nonfinite observed
data are rejected. Original time indexes are preserved for ts, zoo and xts
inputs. No artificial response is imputed.



There might be situations when the practitioner desires to predict \eqn{T+1} conditional on observation till time \eqn{T} in a true out-of-sample fashion. In such circumstances the user can substitute the future value of the dependent variable with an \code{NA}. This way, the code treats the last observation as missing and does not perform backtesting or updating of the coefficients. However, the filter provides us with the necessary quantities to perform prediction. The predicted value \eqn{\hat{y_{T+1}} = E[y_T+1 | F_T]} as well as the predicted variance decomposition can then be extracted using the \link{getLastForecast} method. The other quantities that can be extracted, for example via the \link{as.data.frame} method, will ignore the presence of the last \code{NA} and report results only for the first \eqn{T} observations.

See Catania and Nonejad (2016) for further details.
}
\value{
An object of the class \code{DMA}, see \link{DMA-class}.
}
\references{

Beckmann, J., & Schussler, R. (2014). Forecasting Equity Premia using Bayesian Dynamic Model Averaging (No. 2914). Center for Quantitative Economics (CQE), University of Muenster. \cr

Dangl, T., & Halling, M. (2012). Predictive regressions with time--varying coefficients. \emph{Journal of Financial Economics}, \bold{106(1)}, 157--181.  \doi{doi:10.1016/j.jfineco.2012.04.003}.\cr

Catania, Leopoldo, and Nima Nonejad (2018). "Dynamic Model Averaging for Practitioners in Economics and Finance: The eDMA Package." Journal of Statistical Software, 84(11), 1-39. \doi{10.18637/jss.v084.i11}.\cr

Paye, B.S. (2012). 'Deja vol': Predictive Regressions for Aggregate Stock Market Volatility
Using Macroeconomic Variables.\emph{Journal of Financial Economics}, \bold{106(3)}, 527-546. ISSN 0304-405X.
\doi{doi:10.1016/j.jfineco.2012.06.005}. URL \url{http://www.sciencedirect.com/science/article/pii/S0304405X12001316}.

Prado, R., & West, M. (2010). Time series: modeling, computation, and inference. CRC Press.
}
\author{ Leopoldo Catania }
\examples{

#  Code chunk of Catania and Nonejad (2016) Fast Dynamic Model Averaging
#  for Practitioners in Economics and Finance: The eDMA Package
library(eDMA)

## load data
data("USData")

## do DMA, keep the first three predictors fixed and the intercept
Fit <- DMA(GDPDEF ~ Lag(GDPDEF, 1) + Lag(GDPDEF, 2) + Lag(GDPDEF, 3) +
            Lag(ROUTP, 1) + Lag(UNEMP, 1), data = USData, vDelta = c(0.9,0.95,0.99),
             vKeep = c(1, 2, 3, 4))

Fit

}
