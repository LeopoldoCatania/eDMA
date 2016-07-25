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
\title{class: Class for the DMA class}
\description{
Class for the DMA estimate.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
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
\item{show}{\code{signature(object = "DMA")}}
\item{summary}{\code{signature(object = "DMA")}:
Print a summary of the estimated model. This method accepts the additional argument \code{iBurnPeriod} corresponding at the length of the burn-in period. By default \code{iBurnPeriod = NULL}, i.e. all the sample is used.
}
\item{coef}{\code{signature(object = "DMA")}:
Extract the filtered regressor coefficients. This method accepts the additional argument \code{iBurnPeriod} corresponding at the length of the burn-in period. By default \code{iBurnPeriod = NULL}, i.e. all the sample is used.
}
\item{residuals}{\code{signature(object = "DMA")}:
Extract the residuals of the model. This method accepts the additional argument \code{iBurnPeriod} corresponding at the length of the burn-in period. By default \code{iBurnPeriod = NULL}, i.e. all the sample is used. The additional Boolean argument \code{standardize} controls if standardize residuals should be returned. By default \code{standardize = FALSE}. The additional argument \code{Type} permits to choose between residuals evaluated using DMA or DMS. By default \code{Type = "DMA"}.
}
\item{inclusion.prob}{\code{signature(object = "DMA")}:
Extract the inclusion probabilities of the regressors. This method accepts the additional argument \code{iBurnPeriod} corresponding at the length of the burn-in period. By default \code{iBurnPeriod = NULL}, i.e. all the sample is used.
}
\item{pred.like}{\code{signature(object = "DMA")}:
Extract the predictive log likelihood series. This method accepts the additional argument \code{iBurnPeriod} corresponding at the length of the burn-in period. By default \code{iBurnPeriod = NULL}, i.e. all the sample is used. The additional argument \code{Type} permits to choose between predictive likelihood evaluated using DMA or DMS. By default \code{Type = "DMA"}.
}
}
}
\note{
The \code{as.data.frame()} method permits to extract several estimated quantities. It accepts the two additional arguments: \code{which} with possible values:
\itemize{
  \item{\code{mincpmt}}: Posterior inclusion probability of the predictors
  \item{\code{vsize}}: Expected number of predictors (average size)
  \item{\code{mtheta}}: Filtered estimates of the regression coefficients
  \item{\code{mpmt}}: Posterior probability of the degree of instability
  \item{\code{vyhat}}: Point forecasts
  \item{\code{vLpdfhat}}: Predictive likelihood
  \item{\code{vdeltahat}}: Posterior weighted average of delta
  \item{\code{mvdec}}: representing the y_t variance decomposition. The function returns a T x 5 matrix whose columns contains the variables
  \itemize{
    \item{\code{vtotal}}: total variance
    \item{\code{vobs}}: Observational variance
    \item{\code{vcoeff}}: Variance due to errors in the estimation of the coefficients
    \item{\code{vmod}}: Variance due to model uncertainty
    \item{\code{vtvp}}: Variance due to uncertainty with respect to the choice of the degrees of time--variation in the regression coefficients.
   }
 \item \code{vhighmp_DMS}: Highest posterior model probability.
  \item \code{vhighmpTop01_DMS}: Sum of the 10\% highest posterior model probabilities.
}
and \code{iBurnPeriod} which is an integer indicating the length of the burn-in period. For instance, if \code{iBurnPeriod = 50} the first 50 observations are removed from the output. By default \code{iBurnPeriod = NULL} meaning that all the observations are returned.
}
\author{ Leopoldo Catania & Nima Nonejad }
\keyword{classes}
