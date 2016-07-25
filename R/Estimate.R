
DMA <- function(formula, data, vDelta = c(0.9, 0.95, 0.99), dAlpha = 0.99, vKeep = NULL, bZellnerPrior = FALSE, dG = 100, bParallelize = TRUE, iCores = NULL) {

  Start = Sys.time()

  if (any(class(data) %in% c("ts", "zoo", "xts"))) {
    bTimeSeries = TRUE
  } else {
    bTimeSeries = FALSE
  }

  if (!is(formula, "formula")){
    formula = as.formula(formula)
  }
  if (!is(data, "data.frame")){
    data = as.data.frame(data)
  }

  vY = model.frame(formula , data = data)[, 1]
  mF = model.matrix(formula, data)

  if (bTimeSeries) {
    vDates = as.Date(rownames(mF))
  } else {
    vDates = NULL
  }

  mF = as.matrix(mF)

  if(!is.null(vDates)) vY = xts(vY, order.by = vDates)

  if (is.null(vKeep)) {
    FixedVar = FALSE
    vKeep = -9999
  } else if (vKeep[1] == "KS" ) {
    FixedVar = TRUE
    vKeep = 0:(ncol(mF) - 1)
  } else {
    FixedVar = TRUE
    if (is.character(vKeep))
      vKeep = which(colnames(mF) %in% vKeep)
    vKeep = vKeep - 1
  }

  if (bParallelize) {
    if (is.null(iCores)) {
      iCores = detectCores() - 1
    }
    Est = funcEstimate_Eff_par(vY, mF, vDelta, dAlpha, vKeep, bZellnerPrior, dG, iCores)
  } else {
    Est = funcEstimate_Eff(vY, mF, vDelta, dAlpha, vKeep, bZellnerPrior, dG)
  }

  if (FixedVar){
    vKeep = vKeep + 1
  }

  elapsedTime = Sys.time() - Start

  out <- new("DMA",
             model = list(vDelta = vDelta, dAlpha = dAlpha, bZellnerPrior = bZellnerPrior, dG = dG, elapsedTime = elapsedTime, bParallelize = bParallelize,
                                 FixedVar = FixedVar, vKeep = vKeep, Call = formula ),
             data = list(vY = vY, mF = mF, formula = formula, vDates = vDates),
             Est  = Est)

  return(out)
}


