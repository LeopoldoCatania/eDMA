setClass("DMA", representation(model = "list", data = "list", Est = "list"))

setMethod("show", "DMA", function(object) {

  iK = ncol(object@data$mF)
  iT = nrow(object@data$mF)
  iM = object@Est$iM
  vDelta = object@model$vDelta
  iD = length(vDelta)
  dBeta = object@model$dBeta
  dAlpha = object@model$dAlpha
  bZellnerPrior = object@model$bZellnerPrior
  dG = object@model$dG
  elapsedTime = object@model$elapsedTime
  vKeep = object@model$vKeep
  FixedVar = object@model$FixedVar

  vNames = colnames(object@data$mF)
  if (is.null(vNames))
    vNames = paste("var", 1:iK, sep = ".")

  cat(paste("\n------------------------------------------"))
  cat(paste("\n-        Dynamic Model Averaging         -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification\t")
  cat(paste("\nT     =", iT))
  cat(paste("\nn     =", iK))
  cat(paste("\nd     =", iD))
  cat(paste("\nAlpha =", dAlpha))
  cat(paste("\nBeta  =", dBeta))
  cat(paste("\nModel combinations =", iM))
  cat(paste("\nModel combinations including averaging over delta =", iM*iD))
  cat(paste("\n------------------------------------------"))
  if (bZellnerPrior) {
    cat(paste("\nPrior : Zellner's with degree of shrinkage, dG, equal to", dG))
  } else {
    cat(paste("\nPrior : Multivariate Gaussian with mean vector 0 and covariance matrix equal to: ", dG, " x diag(",iK,")",sep = "" ))
  }
  if (FixedVar) {
    cat(paste("\n\nVariables always included :", paste(vNames[vKeep], collapse = ", ")))
  }
  cat(paste("\n------------------------------------------"))
  cat(paste("\nThe grid for delta:\n"))
  cat(paste("\nDelta = ", paste(vDelta, collapse = ", ")))
  cat(paste("\n------------------------------------------"))
  cat(paste("\n\nElapsed time\t:", round(as.double(elapsedTime, units = "secs"), 2), "secs"))
})

setMethod("summary", "DMA", function(object, iBurnPeriod = NULL) {

  Call = object@model$Call
  vRes = residuals(object, iBurnPeriod = iBurnPeriod)
  mTheta = coef(object, iBurnPeriod = iBurnPeriod)
  mProb  = inclusion.prob(object, iBurnPeriod)

  ETheta   = apply(mTheta, 2, mean)
  SDTheta  = apply(mTheta, 2, sd)
  EPTheta   = apply(mProb, 2, mean)
  SDPTheta = apply(mProb, 2, sd)

  iCeil = ceiling(0.1*length(ETheta))
  Top10 = names(sort(EPTheta, decreasing = TRUE)[1:iCeil])

  mCoefMat = round(cbind(ETheta, SDTheta, EPTheta, SDPTheta),2)
  colnames(mCoefMat) = c("E[theta_t]", "SD[theta_t]", "E[P(theta_t)]", "SD[P(theta_t)]")

  ForcPerf = round(BacktestDMA(object, iBurnPeriod),3)

  vResStat = round(c(min(vRes), quantile(vRes, c(0.25,0.50,0.75)), max(vRes)),4);
  names(vResStat) = c("Min", "1Q", "Median", "3Q", "Max")
  mVariance = as.data.frame(object, which = "mvdec", iBurnPeriod = iBurnPeriod)
  vVariance = apply(mVariance,2,mean)
  vVarPerc = round((vVariance/vVariance[1])[-1]*100, 2)
  cat("\nCall:\n DMA(formula = ", paste(deparse(Call), sep = "\n"), ")")
  cat("\n\nResiduals:\n")
  print(vResStat)
  cat("\nCoefficients:\n")
  print(mCoefMat)
  cat("\nVariance contribution (in percentage points):\n")
  print(vVarPerc)
  cat("\nTop 10% included regressors:\t")
  cat(paste(Top10, collapse = ", "))
  cat("\n\nForecast Performance:\n")
  print(ForcPerf)

})

setMethod("as.data.frame", signature(x = "DMA"), function(x, which, iBurnPeriod = NULL) {
  if (!is.character(which) || length(which) != 1L || is.na(which)) stop("which must name one output series.")
  key <- if (which == "mtheta") "mmhat" else which
  out <- x@Est[[key]]
  if (!is.matrix(out) || nrow(out) != nrow(x@data$mF)) stop("Unsupported output series: ", which)
  rows <- .dma_rows(x, iBurnPeriod)
  out <- out[rows, , drop = FALSE]
  if (key == "mvdec") colnames(out) <- c("vtotal", "vobs", "vcoeff", "vmod", "vtvp")
  if (key %in% c("mincpmt", "mmhat", "mcovhat")) colnames(out) <- colnames(x@data$mF)
  if (key %in% c("mpmt", "mpmt_pred", "mLpdf")) colnames(out) <- x@model$vDelta
  .dma_indexed(out, if (is.null(x@data$vDates)) NULL else x@data$vDates[rows])
})

setMethod("plot", signature(x = "DMA", y = "missing"), function(x, which = NULL, iBurnPeriod = NULL, ...) {
  rows <- .dma_rows(x, iBurnPeriod)
  if (is.null(which)) {
    if (!interactive()) which <- "vyhat"
    else {
      choice <- menu(PlotMenu("DMA"))
      if (!choice) return(invisible(x))
      which <- PlotNumber2Label(choice)
    }
  }
  values <- as.matrix(as.data.frame(x, which = which, iBurnPeriod = iBurnPeriod))
  dates <- if (is.null(x@data$vDates)) rows else x@data$vDates[rows]
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  count <- ncol(values)
  titles <- colnames(values)
  if (is.null(titles)) titles <- if (count == 1L) which else paste(which, seq_len(count))
  perpage <- min(count, 10L)
  par(mfrow = c(ceiling(perpage / 2), if (perpage == 1L) 1L else 2L),
      mar = c(2.5, 3, 2, 1))
  for (j in seq_len(count)) {
    y <- values[, j]
    finite <- is.finite(y)
    lim <- if (any(finite)) range(y[finite]) else c(0, 1)
    y[!finite] <- NA_real_
    args <- utils::modifyList(list(x = dates, y = y, type = "l", xlab = "", ylab = "",
                                  main = titles[j], ylim = lim), list(...))
    do.call(graphics::plot, args)
    grid()
    if (which %in% c("vyhat", "vyhat_DMS")) points(dates, as.numeric(x@data$vY)[rows], col = "red")
  }
  invisible(x)
})

setMethod("coef", signature(object = "DMA"), function(object, iBurnPeriod = NULL) {
  as.data.frame(object, which = "mtheta", iBurnPeriod = iBurnPeriod)
})
setMethod("residuals", signature(object = "DMA"), function(object, standardize = FALSE, type = "DMA", iBurnPeriod = NULL) {
  if (!is.character(type) || length(type) != 1L || is.na(type) || !type %in% c("DMA", "DMS")) stop("type must be DMA or DMS.")
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) stop("standardize must be TRUE or FALSE.")
  out <- as.data.frame(object, which = if (type == "DMA") "veps" else "veps_DMS", iBurnPeriod = iBurnPeriod)
  if (standardize) {
    scale <- sd(as.numeric(out))
    if (!is.finite(scale) || scale <= 0) stop("Residual standard deviation is zero or undefined.")
    out <- out / scale
  }
  out
})
inclusion.prob <- function(object, ...) UseMethod("inclusion.prob")
setMethod("inclusion.prob", signature(object = "DMA"), function(object, iBurnPeriod = NULL) {
  as.data.frame(object, which = "mincpmt", iBurnPeriod = iBurnPeriod)
})
pred.like <- function(object, ...) UseMethod("pred.like")
setMethod("pred.like", signature(object = "DMA"), function(object, type = "DMA", iBurnPeriod = NULL) {
  if (!is.character(type) || length(type) != 1L || is.na(type) || !type %in% c("DMA", "DMS")) stop("type must be DMA or DMS.")
  as.data.frame(object, which = if (type == "DMA") "vLpdfhat" else "vLpdfhat_DMS", iBurnPeriod = iBurnPeriod)
})
getLastForecast <- function(object, ...) UseMethod("getLastForecast")
setMethod("getLastForecast", signature(object = "DMA"), function(object) {
  out <- object@Est$LastForecast
  if (!isTRUE(out$bForecast)) stop("No future forecast is available: supply a final NA response with known predictors to DMA().")
  out[c("PointForecast", "VarianceDecomposition", "ForecastDate")]
})
