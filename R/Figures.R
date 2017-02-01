PlotMenu <- function(sClass) {
    if (sClass == "DMA") {
        vplotMenu = c("Point forecast",
                      "Predictive likelihood",
                      "Posterior weighted average of delta",
                      "Posterior inclusion probabilities of the predictors",
                      "Posterior probabilities of the forgetting factors",
                      "Filtered estimates of the regression coefficients",
                      "Variance decomposition",
                      "Observational variance",
                      "Variance due to errors in the estimation of the coefficients, theta",
                      "Variance due to model uncertainty",
                      "Variance due to uncertainty with respect to the choice of\n the degrees of time-variation in the regression coefficients",
                      "Expected number of predictors (average size)",
                      "Number of predictors (highest posterior model probability) (DMS)",
                      "Highest posterior model probability (DMS)",
                      "Point forecasts (highest posterior model probability) (DMS)",
                      "Predictive likelihood (highest posterior model probability) (DMS)")
    }
    return(vplotMenu)
}

PlotLabel2Number <- function(sWhich) {
    vplotMenu = c("vyhat", "vLpdfhat", "vdeltahat", "mincpmt", "mpmt", "mmhat", "mvdec", "vobs", "vcoeff",
                  "vmod", "vtvp", "vsize", "vsize_DMS", "vhighmp_DMS", "vyhat_DMS", "vLpdfhat_DMS")
    return(which(sWhich == vplotMenu))
}

PlotNumber2Label <- function(PlotType) {
  vplotMenu = c("vyhat", "vLpdfhat", "vdeltahat", "mincpmt", "mpmt", "mmhat", "mvdec", "vobs", "vcoeff",
                "vmod", "vtvp", "vsize", "vsize_DMS", "vhighmp_DMS", "vyhat_DMS", "vLpdfhat_DMS")
  return(vplotMenu[PlotType])
}


LegendFunDMAclass <- function(PlotType) {
    if (PlotType == 1 | PlotType == 15)
        legend("topleft", legend = c("Predicted", "Realised"), cex = 1, text.font = 2, col = c("black", "red"), lty = c(1, NA), bg = "white",
            pch = c(NA, 1))
}
LegendFun <- function(sClass, PlotType) {
    if (sClass == "DMA")
        LegendFunDMAclass(PlotType)
}

TitleFunClass <- function(PlotType) {
    vNames = PlotMenu("DMA")
    return(vNames[PlotType])
}
TitleFun <- function(sClass, PlotType) {
    if (sClass == "DMA")
        TitleFunClass(PlotType)
}
