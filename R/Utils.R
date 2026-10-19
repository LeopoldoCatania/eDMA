.dma_rows <- function(object, burn = NULL) {
  n <- nrow(object@data$mF)
  if (is.null(burn)) burn <- 0L
  if (!is.numeric(burn) || length(burn) != 1L || !is.finite(burn) ||
      burn != floor(burn) || burn < 0 || burn >= n)
    stop("iBurnPeriod must be an integer between 0 and the number of observations minus one.")
  seq.int(burn + 1L, n)
}

.dma_indexed <- function(x, dates) {
  if (is.null(dates)) return(x)
  if (inherits(dates, "Date") || inherits(dates, "POSIXt")) xts::xts(x, order.by = dates)
  else zoo::zoo(x, order.by = dates)
}
