BacktestDMA <- function(object, iBurnPeriod = NULL){

  vY = object@data$vY
  iT = length(vY)

  if (is.null(iBurnPeriod)) {
    iBurnPeriod = round(iT / 3)
    warnings(paste("No burin in period has been selected. ", iBurnPeriod,
                   " periods (approximately 1/3 of the sample) has been automatically selected.", sep = ""))
  }

  vBurnIN = 1:iBurnPeriod

  vEps_DMA  = (object@Est$veps)[-vBurnIN]
  vLpdf_DMA = object@Est$vLpdfhat[-vBurnIN]

  vEps_DMS  = (object@Est$vyhat_DMS - vY)[-vBurnIN]
  vLpdf_DMS = object@Est$vLpdfhat_DMS[-vBurnIN]

  mOut = matrix(NA, 3, 2, dimnames = list(c("MSE", "MAD", "Log-predictive Likehood"), c("DMA", "DMS")))

  mOut["MSE", "DMA"] = mean(vEps_DMA^2.0)
  mOut["MSE", "DMS"] = mean(vEps_DMS^2.0)
  mOut["MAD", "DMA"] = mean(abs(vEps_DMA))
  mOut["MAD", "DMS"] = mean(abs(vEps_DMS))

  mOut["Log-predictive Likehood", "DMA"] = sum(vLpdf_DMA)
  mOut["Log-predictive Likehood", "DMS"] = sum(vLpdf_DMS)

  return(mOut)

}


