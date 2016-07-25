setClass("DMA", representation(model = "list", data = "list", Est = "list"))

setMethod("show", "DMA", function(object) {

  iK = ncol(object@data$mF)
  iT = nrow(object@data$mF)
  iM = object@Est$iM
  vDelta = object@model$vDelta
  iD = length(vDelta)
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
  cat(paste("\n-        Dynamic Model Ageraging         -"))
  cat(paste("\n------------------------------------------"))
  cat("\n\nModel Specification\t")
  cat(paste("\nT     =", iT))
  cat(paste("\nn     =", iK))
  cat(paste("\nd     =", iD))
  cat(paste("\nAlpha =", dAlpha))
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
  object = x
  vY  = object@data$vY
  Est = object@Est

  if (is(vY, "xts")) {
    vDates = index(vY)
  } else {
    vDates = 1:length(vY)
  }

  if (!is.null(iBurnPeriod)) {
    vY = vY[-c(1:iBurnPeriod)]
    vDates = vDates[-c(1:iBurnPeriod)]
    for(v in 1:length(Est)){
      if (is.matrix(Est[[v]])) {
        Est[[v]] = Est[[v]][-c(1:iBurnPeriod), ]
      } else {
        Est[[v]] = Est[[v]][-c(1:iBurnPeriod)]
      }
    }
  }

  if (which == "mvdec") {
    Out =  Est[["mvdec"]]
    colnames(Out) = c("vtotal", "vobs", "vcoeff", "vmod", "vtvp")
  } else if (which == "mtheta") {
    Out = Est[["mmhat"]]
  } else {
    Out = Est[[which]]
  }

  if (which == "mincpmt" | which == "mtheta") {
    vNames = colnames(object@data$mF)
    if (is.null(vNames))
      vNames = paste("var", 1:ncol(Out), sep = ".")
    colnames(Out) = vNames
  }

  if (which == "mpmt") {
    vDelta = object@model$vDelta
    colnames(Out) = vDelta
  }

  if (is(vY, "xts")) {
    Out = xts(Out, vDates)
  }

  return(Out)
})

setMethod("plot", signature(x = "DMA", y = "missing"), function(x, which = NULL, iBurnPeriod = NULL, ...) {

  Est = x@Est
  vY = x@data$vY
  mF = x@data$mF
  vDates = x@data$vDates

  if (is.null(vDates)) {
    vDates = 1:length(vY)
  }
  if (!is.null(iBurnPeriod)) {
    vY = vY[-c(1:iBurnPeriod)]
    mF = mF[-c(1:iBurnPeriod), ]
    vDates = vDates[-c(1:iBurnPeriod)]

    for(v in 1:length(Est)){
      if (is.matrix(Est[[v]])) {
        Est[[v]] = Est[[v]][-c(1:iBurnPeriod), ]
      } else {
        Est[[v]] = Est[[v]][-c(1:iBurnPeriod)]
      }
    }
  }

  iT = length(vY)

  PlotType = 1
  while (PlotType > 0) {
    if (is.null(which)) {
      cat(paste("Print 1-11 or 0 to exit"))
      PlotType = menu(PlotMenu("DMA"))

      if (PlotType > 0) {
        PlotLabel   = PlotNumber2Label(PlotType)
        series2plot = Est[[PlotLabel]]
      }

      if(PlotLabel == "mvdec")
        colnames(series2plot) = c("vtotal", "vobs", "vcoeff", "vmod", "vtvp")

    } else {
      if (which == "mtheta")
        which = "mmhat"
      series2plot = as.matrix(Est[[which]])
      if (is.null(series2plot))
        stop(paste("which =", which, "is not supported."))
      PlotType = PlotLabel2Number(which)
    }

    if (PlotType > 0) {

      if (dev.cur() != 1)
        dev.off()

      if (ncol(series2plot) == 1) {
        sTitle = TitleFun("DMA", PlotType)
        plot(vDates, vY, type = "n", xaxt = "n", xlab = "", ylab = "", las = 1, ylim = c(min(series2plot[-1]), max(series2plot[-1])),
             main = sTitle)
        grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
        lines(vDates[-1], series2plot[-1], col = "black")
        if (PlotType == 1 | PlotType == 15)
          points(vDates, vY, col = "red")
        if( !is.numeric(vDates) ){
          axis.Date(1, at = seq(min(vDates), max(vDates) + 300, "year"))
          axis.Date(1, at = seq(min(vDates), max(vDates) + 300, "quarter"), labels = FALSE, tcl = -0.2)
        } else {
          foo = vDates[c(1, seq(0, iT, ceiling((iT)/20))[-1])]
          axis(1, at = foo, labels = foo)
        }
        LegendFun("DMA", PlotType)

      } else {
        iV = ncol(series2plot)
        VarNames = colnames(series2plot)
        if (PlotType == 4)
          VarNames = colnames(mF)
        if (is.null(VarNames))
          VarNames = 1:iV

        if (iV <= 5) {
          layout(matrix(1:iV, iV, 1), heights = c(rep(2, iV - 1), 2.5))
          for (i in 1:(iV)) {
            if (i == 1)
              par(mar = c(0, 4, 0.1, 2))
            if (i != 1 & i != iV)
              par(mar = c(0, 4, 0, 2))
            if (i == iV)
              par(mar = c(3, 4, 0, 2))

            if (any(PlotType == c(4, 5))) {
              vLim = c(0, 1)
            } else {
              vLim = c(min(series2plot[, i]), max(series2plot[, i]))
            }

            plot(vDates, series2plot[, i], type = "n", xaxt = "n", xlab = "", ylab = "", las = 1, ylim = vLim)
            grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
            lines(vDates[-1], series2plot[-1, i], col = "black")
            axis(4, at = mean(vLim), labels = VarNames[i], tick = F, padj = -1)
          }
          if (is(vY, "xts")) {
            axis.Date(1, at = seq(min(vDates), max(vDates), "year"))
            axis.Date(1, at = seq(min(vDates), max(vDates), "quarter"), labels = FALSE, tcl = -0.2)
          } else {
            foo = vDates[c(1, seq(0, iT, ceiling((iT)/20))[-1])]
            axis(1, at = foo, labels = foo)
          }

        } else if (iV == 6 | iV == 8 ) {
          if (iV == 6) {
            plotSeq = seq(1, iV + 1, 3)
            layout(matrix(1:6, 3, 2), heights = c(rep(2, 2), 2.5, rep(2, 2), 2.5))
          }
          if (iV == 8) {
            plotSeq = seq(1, iV + 1, 4)
            layout(matrix(1:8, 4, 2), heights = c(rep(2, 3), 2.5, rep(2, 3), 2.5))
          }

          for (i in 1:iV) {
            if (i <= iV) {
              if (any(i == plotSeq))
                par(mar = c(0, 4, 0.1, 2))
              if (all(i != plotSeq) & all(i != plotSeq - 1))
                par(mar = c(0, 4, 0, 2))
              if (any(i == plotSeq - 1))
                par(mar = c(3, 4, 0, 2))

              if (any(PlotType == c(4, 5))) {
                vLim = c(0, 1)
              } else {
                vLim = c(min(series2plot[, i]) * 1.1, max(series2plot[, i]) * 1.1)
              }

              plot(vDates, series2plot[, i], type = "n", xaxt = "n", xlab = "", ylab = "", las = 1, ylim = vLim)
              grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
              lines(vDates[-1], series2plot[-1, i], col = "black")
              axis(4, at = mean(vLim), labels = VarNames[i], tick = F, padj = -1)

              if (any(i == plotSeq - 1) | (i == iV)) {

                if (is(vY, "xts")) {
                  axis.Date(1, at = seq(min(vDates), max(vDates), "year"))
                  axis.Date(1, at = seq(min(vDates), max(vDates), "quarter"), labels = FALSE, tcl = -0.2)
                } else {
                  foo = vDates[c(1, seq(0, iT, ceiling((iT)/20))[-1])]
                  axis(1, at = foo, labels = foo)
                }

              }
            }
          }

        } else {
          nPlot = ceiling(iV/10)
          plotSeq = seq(1, iV + 1, 5)
          Start = 1
          PlotType2 = ""

          for (j in 1:nPlot) {
            if (PlotType2 != "0") {

              layout(matrix(1:10, 5, 2), heights = c(rep(2, 4), 2.5, rep(2, 4), 2.5))

              for (i in Start:(Start + 9)) {
                if (i <= iV) {
                  if (any(i == plotSeq))
                    par(mar = c(0, 4, 0.1, 2))
                  if (all(i != plotSeq) & all(i != plotSeq - 1))
                    par(mar = c(0, 4, 0, 2))
                  if (any(i == plotSeq - 1))
                    par(mar = c(3, 4, 0, 2))

                  if (any(PlotType == c(4, 5))) {
                    vLim = c(0, 1)
                  } else {
                    vLim = c(min(series2plot[, i]) * 1.1, max(series2plot[, i]) * 1.1)
                  }

                  plot(vDates, series2plot[, i], type = "n", xaxt = "n", xlab = "", ylab = "", las = 1, ylim = vLim)
                  grid(nx = 10, ny = 10, col = "gray", lty = "dotted")
                  lines(vDates[-1], series2plot[-1, i], col = "black")
                  axis(4, at = mean(vLim), labels = VarNames[i], tick = F, padj = -1)

                  if (any(i == plotSeq - 1) | (i == iV)) {

                    if (is(vY, "xts")) {
                      axis.Date(1, at = seq(min(vDates), max(vDates), "year"))
                      axis.Date(1, at = seq(min(vDates), max(vDates), "quarter"), labels = FALSE, tcl = -0.2)
                    } else {
                      foo = vDates[c(1, seq(0, iT, ceiling((iT)/20))[-1])]
                      axis(1, at = foo, labels = foo)
                    }

                  }
                }
              }
              Start = Start + 10
              if (j < nPlot)
                PlotType2 = readline("Print enter for next figures or 0 to exit\n:")

            }
          }

        }

      }
    }
    if (!is.null(which))
      PlotType = 0
  }
})

setMethod("coef", signature(object = "DMA"), function(object, iBurnPeriod = NULL ){
  mTheta = as.data.frame(object, which = "mtheta", iBurnPeriod = iBurnPeriod)
  return(mTheta)
})

setMethod("residuals",  signature(object = "DMA"), function(object, standardize = FALSE,  Type = "DMA" ,iBurnPeriod = NULL ) {
  if (Type=="DMA") vres = as.data.frame(object, which = "veps", iBurnPeriod = iBurnPeriod)
  if (Type=="DMS") vres = as.data.frame(object, which = "veps_DMS", iBurnPeriod = iBurnPeriod)

  if (standardize) vres = vres/sd(vres)

  return(vres)
})

inclusion.prob = function(object, ...) {
  UseMethod("inclusion.prob")
}

setMethod("inclusion.prob",  signature(object = "DMA"), function(object, iBurnPeriod = NULL ) {
  mProb = as.data.frame(object, which = "mincpmt", iBurnPeriod = iBurnPeriod)
  return(mProb)
})

pred.like = function(object, ...) {
  UseMethod("pred.like")
}

setMethod("pred.like",  signature(object = "DMA"), function(object, Type = "DMA", iBurnPeriod = NULL ) {

  if (Type=="DMA") vPredLike = as.data.frame(object, which = "vLpdfhat", iBurnPeriod = iBurnPeriod)
  if (Type=="DMS") vPredLike = as.data.frame(object, which = "vLpdfhat_DMS", iBurnPeriod = iBurnPeriod)

  return(vPredLike)
})
