DMA <- function(formula, data = NULL, vDelta = c(0.9, 0.95, 0.99), dAlpha = 0.99, vKeep = NULL, bZellnerPrior = FALSE, dG = 100, bParallelize = TRUE, iCores = NULL,
                dBeta = 1.0) {

  Start = Sys.time()

  if (any(class(data) %in% c("ts", "zoo", "xts"))) {
    bTimeSeries = TRUE
  } else {
    bTimeSeries = FALSE
  }

  if (!is(formula, "formula")) {
    formula = as.formula(formula)
  }

  if (!is.null(data)) {
    if (!is(data, "data.frame")) {
      data = as.data.frame(data)
    }
  }

  ModelFrame = model.frame(formula , data = data, na.action = na.omit_new)

  vY = as.numeric(ModelFrame[, 1])
  mF = model.matrix(formula, data = ModelFrame)

  if (is.na(tail(vY, 1))) {
    bForecast = TRUE
    vY[length(vY)] = mean(na.omit(vY[-length(vY)]))
  } else {
    bForecast = FALSE
  }

  if (bTimeSeries) {
    vDates = as.Date(rownames(mF))
  } else {
    vDates = NULL
  }

  mF = as.matrix(mF)

  if (!is.null(vDates)) vY = xts(vY, order.by = vDates)

  if (is.null(vKeep)) {
    FixedVar = FALSE
    vKeep = -9999
  } else if (vKeep[1] == "KS" ) {
    FixedVar = TRUE
    vKeep = 0:(ncol(mF) - 1)
  } else {
    FixedVar = TRUE
    if (is.character(vKeep)) {
      vKeep = which(colnames(mF) %in% vKeep)
    }
    vKeep = vKeep - 1
  }

  if (bParallelize) {
    if (is.null(iCores)) {
      iCores = detectCores() - 1
      ## reduce the number of cores if model complexity is low
      ## model complexity is defined low if less than 100 models
      ## are considered.
      if ((((ncol(mF) - length(vKeep))^2 - 1) * length(vDelta) <= 100) & (iCores > 2)) {
        iCores = 2
      }
    }
    Est = funcEstimate_Eff_par(vY, mF, vDelta, dAlpha, vKeep, dBeta, bZellnerPrior, dG, iCores)
  } else {
    Est = funcEstimate_Eff(vY, mF, vDelta, dAlpha, vKeep, dBeta, bZellnerPrior, dG)
  }

  if (FixedVar) {
    vKeep = vKeep + 1
  }

  if (bForecast) {

    dPointForecast = c(tail(Est$vyhat, 1))
    vVariancePrediction = c(tail(Est$mvdec, 2)[1, ])
    names(vVariancePrediction) = c("vtotal", "vobs", "vcoeff", "vmod", "vtvp")

    # remove last filtered values given that we cannot do backtest
    for (i in 1:length(Est)) {
      if (is(Est[[i]], "matrix")) {
        Est[[i]] = Est[[i]][-nrow(Est[[i]]), ,drop = FALSE]
      }
    }

    vY = vY[-length(vY)]
    mF = mF[-nrow(mF), ,drop = FALSE]

    Est[["LastForecast"]] = list(PointForecast = dPointForecast,
                                 VarianceDecomposition = vVariancePrediction,
                                 bForecast = bForecast)

  } else {

    dPointForecast = NULL
    vVariancePrediction = NULL

    Est[["LastForecast"]] = NULL

  }

  elapsedTime = Sys.time() - Start

  out <- new("DMA",
             model = list(vDelta = vDelta, dAlpha = dAlpha, bZellnerPrior = bZellnerPrior, dG = dG, dBeta = dBeta,
                          elapsedTime = elapsedTime, bParallelize = bParallelize,
                          FixedVar = FixedVar, vKeep = vKeep, Call = formula ),
             data = list(vY = vY, mF = mF, formula = formula, vDates = vDates),
             Est  = Est)

  return(out)
}


