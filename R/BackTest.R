BacktestDMA <- function(object, iBurnPeriod = NULL) {
  if (is.null(iBurnPeriod)) {
    iBurnPeriod <- min(round(nrow(object@data$mF) / 3), nrow(object@data$mF) - 1L)
    warning("No burn-in period selected; using ", iBurnPeriod,
            " periods (approximately one third of the sample).", call. = FALSE)
  }
  rows <- .dma_rows(object, iBurnPeriod)
  out <- matrix(NA_real_, 3L, 2L,
                dimnames = list(c("MSE", "MAD", "Log-predictive Likelihood"), c("DMA", "DMS")))
  for (type in colnames(out)) {
    suffix <- if (type == "DMS") "_DMS" else ""
    e <- object@Est[[paste0("veps", suffix)]][rows, 1L]
    lp <- object@Est[[paste0("vLpdfhat", suffix)]][rows, 1L]
    out[, type] <- c(mean(e^2), mean(abs(e)), sum(lp))
  }
  out
}
