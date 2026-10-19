\name{DMA-class}
\Rdversion{0.1}
\docType{class}
\alias{DMA-class}
\alias{show,DMA-method}
\alias{plot,DMA,missing-method}
\alias{as.data.frame,DMA-method}
\alias{coef,DMA-method}
\alias{residuals,DMA-method}
\alias{summary,DMA-method}
\alias{inclusion.prob}
\alias{inclusion.prob,DMA-method}
\alias{pred.like}
\alias{pred.like,DMA-method}
\alias{getLastForecast}
\alias{getLastForecast,DMA-method}
\title{class: Class for the DMA class}
\description{
Class for the DMA estimate.
}
\section{Objects from the Class}{Objects are created by \code{DMA()}.}
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"list"}
    Contains information about the DMA specification.}
    \item{\code{data}:}{Object of class \code{"list"}
    Contains the data given to the \code{\link{DMA}} function.}
     \item{\code{Est}:}{Object of class \code{"list"}
    Contains the estimated quantities.}
  }
}
\section{Methods}{
\describe{
\item{as.data.frame}{\code{signature(object = "DMA")}:
Extracts estimated quantities, (see note). }
\item{plot}{\code{signature(x = "DMA", y = "missing")}:
Plots estimated quantities. }
\item{show}{\code{signature(object = "DMA")}}.
\item{summary}{\code{signature(object = "DMA")}:
Print a summary of the estimated model. This method accepts the additional argument \code{iBurnPeriod} corresponding at the length of the burn-in period. By default \code{iBurnPeriod = NULL}, i.e. all the sample is used.
}
\item{coef}{\code{signature(object = "DMA")}:
Extract the filtered regressor coefficients. This method accepts the additional argument \code{iBurnPeriod} corresponding at the length of the burn-in period. By default \code{iBurnPeriod = NULL}, i.e. all the sample is used.
}
\item{residuals}{\code{signature(object = "DMA")}:
Extract the residuals of the model. This method accepts the additional argument \code{iBurnPeriod} corresponding at the length of the burn-in period. By default \code{iBurnPeriod = NULL}, i.e. all the sample is used. The additional Boolean argument \code{standardize} controls if standardize residuals should be returned. By default \code{standardize = FALSE}. The additional argument \code{type} permits to choose between residuals evaluated using DMA or DMS. By default \code{type = "DMA"}.
}
\item{inclusion.prob}{\code{signature(object = "DMA")}:
Extract the inclusion probabilities of the regressors. This method accepts the additional argument \code{iBurnPeriod} corresponding at the length of the burn-in period. By default \code{iBurnPeriod = NULL}, i.e. all the sample is used.
}
\item{pred.like}{\code{signature(object = "DMA")}:
Extract the predictive log-likelihood series. This method accepts the additional argument \code{iBurnPeriod} corresponding at the length of the burn-in period. By default \code{iBurnPeriod = NULL}, i.e. all the sample is used. The additional argument \code{type} permits to choose between predictive likelihood evaluated using DMA or DMS. By default \code{type = "DMA"}.
}
\item{getLastForecast}{\code{signature(object = "DMA")}:
If the last observation of the dependent variable was \code{NA}, i.e. the practitioner desidered to predict \eqn{Y_{T+1}} having a sample of length \eqn{T} (without backtesting the result), this method can be used to extract the predicted value \eqn{\hat{y_T+1} = E[y_{T+1} | F_T]} as well as the predicted variance decomposition using the prediction weights and Student-t variance convention documented in \code{DMA}.
}

}
}
\note{
For compatibility, \code{as.data.frame()} returns a matrix for ordinary data,
an xts object for Date/POSIX indexes, or a zoo object for other time indexes.
It does not return a base data frame. Dimensions are preserved for a single
column or a single retained observation. Burn-in must be an integer from
zero to one less than the observed sample size. All rows, including the
first, now contain defined prediction and filtering results.

Residual standardization divides by the sample standard deviation over the
retained residuals; it is not conditional-volatility standardization.
\code{type} must be \code{"DMA"} or \code{"DMS"}. Plotting preserves the
current graphics device and restores its parameters. In noninteractive use,
the default plot is the point forecast. Nonfinite variances are omitted from
plot lines. \code{getLastForecast()} returns \code{PointForecast},
\code{VarianceDecomposition} and \code{ForecastDate} (NULL without an index),
and gives an informative error if there is no forecast-only row.

Additional extraction names are \code{mpmt_pred} for prediction probabilities
over delta, \code{mLpdf} for conditional predictive log likelihoods,
\code{vdf} for predictive degrees of freedom, and \code{mcovhat} for diagonal
mixture variances of the plug-in filtered coefficient distributions, including
between-component mean dispersion and zero coefficients in absent models.
The latter are not exact conjugate Student-t posterior variances.
Predictive variance components refer to the same distribution as the forecast
and are calculated before observing the current response. They include the
Student-t variance multiplier; total variance can be infinite.

DMS point forecasts use nested delta/model selection. The marginal posterior
model diagnostics need not identify that same model. The top-ten-percent
statistic sums the largest \eqn{\lceil 0.1M\rceil} marginal probabilities.


The \code{as.data.frame()} method permits to extract several estimated quantities. It accepts the two additional arguments: \code{which} with possible values:
\itemize{
  \item{\code{mincpmt}}: Posterior inclusion probabilities of the predictors.
  \item{\code{vsize}}: Expected number of predictors (average size).
  \item{\code{mtheta}}: Filtered estimates of the regression coefficients.
  \item{\code{mpmt}}: Posterior probability of the degree of instability.
  \item{\code{vyhat}}: Point forecasts.
  \item{\code{vLpdfhat}}: Predictive log-likelihood.
  \item{\code{vdeltahat}}: Posterior weighted average of delta.
  \item{\code{mvdec}}: representing the y_t variance decomposition. The function returns a T x 5 matrix whose columns contains the variables.
  \itemize{
    \item{\code{vtotal}}: total variance.
    \item{\code{vobs}}: Observational variance.
    \item{\code{vcoeff}}: Variance due to errors in the estimation of the coefficients.
    \item{\code{vmod}}: Variance due to model uncertainty.
    \item{\code{vtvp}}: Variance due to uncertainty with respect to the choice of the degrees of time--variation in the regression coefficients.
   }
 \item \code{vhighmp_DMS}: Highest posterior model probability.
  \item \code{vhighmpTop01_DMS}: Sum of the 10\% highest posterior model probabilities.
}
and \code{iBurnPeriod} which is an integer indicating the length of the burn-in period. For instance, if \code{iBurnPeriod = 50} the first 50 observations are removed from the output. By default \code{iBurnPeriod = NULL} meaning that all the observations are returned.
}
\references{
Catania, Leopoldo, and Nima Nonejad (2018). "Dynamic Model Averaging for Practitioners in Economics and Finance: The eDMA Package." Journal of Statistical Software, 84(11), 1-39. \doi{10.18637/jss.v084.i11}.\cr
}
\author{ Leopoldo Catania }
\keyword{classes}
