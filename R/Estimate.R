DMA <- function(formula, data = NULL, vDelta = c(0.9, 0.95, 0.99), dAlpha = 0.99,
                vKeep = NULL, bZellnerPrior = FALSE, dG = 100,
                bParallelize = TRUE, iCores = NULL, dBeta = 1.0,
                dS0 = 1.0, dN0 = 2.0) {
  start <- Sys.time()
  scalar <- function(x, label, positive = FALSE, discount = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        (positive && x <= 0) || (discount && (x <= 0 || x > 1)))
      stop(label, " must be a finite numeric scalar", if (discount) " in (0, 1]." else if (positive) " greater than zero." else ".")
  }
  for (nm in c("bZellnerPrior", "bParallelize")) {
    x <- get(nm)
    if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(nm, " must be TRUE or FALSE.")
  }
  if (!is.numeric(vDelta) || !length(vDelta) || any(!is.finite(vDelta)) || any(vDelta <= 0 | vDelta > 1))
    stop("vDelta must be a nonempty numeric vector with values in (0, 1].")
  scalar(dAlpha, "dAlpha", discount = TRUE)
  scalar(dBeta, "dBeta", discount = TRUE)
  scalar(dG, "dG", positive = TRUE)
  scalar(dS0, "dS0", positive = TRUE)
  scalar(dN0, "dN0", positive = TRUE)
  if (!is.null(iCores)) {
    scalar(iCores, "iCores", positive = TRUE)
    if (iCores != floor(iCores) || iCores > 65536) stop("iCores must be an integer no greater than 65536.")
  }
  dates <- if (inherits(data, "ts")) as.numeric(stats::time(data)) else if (inherits(data, "zoo")) zoo::index(data) else NULL
  if (!inherits(formula, "formula")) formula <- as.formula(formula)
  if (!is.null(data)) {
    data <- as.data.frame(data)
    rownames(data) <- as.character(seq_len(nrow(data)))
  }
  frame <- model.frame(formula, data = data, na.action = stats::na.pass)
  if (!nrow(frame)) stop("No observations are available.")
  response <- stats::model.response(frame)
  if (!is.numeric(response) || !is.null(dim(response))) stop("The response must be a numeric vector.")
  if (!is.null(stats::model.offset(frame))) stop("Formula offsets are not supported.")
  forecast <- is.na(tail(response, 1))
  complete <- stats::complete.cases(frame)
  if (forecast) {
    final <- frame[nrow(frame), , drop = FALSE]
    final[[1L]] <- 0
    if (!all(stats::complete.cases(final))) stop("The forecast row must have complete predictors.")
    complete[nrow(frame)] <- TRUE
  }
  frame <- frame[complete, , drop = FALSE]
  if (nrow(frame) <= as.integer(forecast)) stop("At least one observed response is required.")
  vY <- as.numeric(stats::model.response(frame))
  mF <- model.matrix(attr(frame, "terms"), data = frame)
  if (!ncol(mF) || any(!is.finite(mF)) || any(!is.finite(vY[seq_len(length(vY) - forecast)])))
    stop("The design must have at least one column and all observed data must be finite.")
  if (!is.null(dates)) dates <- dates[as.integer(rownames(frame))]
  fixed <- !is.null(vKeep)
  if (!fixed) {
    keep <- numeric(0)
  } else if (identical(vKeep, "KS")) {
    keep <- seq_len(ncol(mF))
  } else {
    if (!length(vKeep) || anyNA(vKeep) || anyDuplicated(vKeep)) stop("vKeep must contain unique, nonmissing names or indices.")
    if (is.character(vKeep)) {
      keep <- match(vKeep, colnames(mF))
      if (anyNA(keep)) stop("Unknown vKeep column: ", paste(vKeep[is.na(keep)], collapse = ", "))
    } else {
      if (!is.numeric(vKeep) || any(!is.finite(vKeep)) || any(vKeep != floor(vKeep)) || any(vKeep < 1 | vKeep > ncol(mF)))
        stop("vKeep indices must be integers between 1 and the number of design columns.")
      keep <- vKeep
    }
  }
  free <- ncol(mF) - length(keep)
  if (free > 20) stop("At most 20 optional design columns are supported; use vKeep to reduce the model space.")
  models <- 2^free - as.integer(!fixed)
  if (is.null(iCores)) {
    available <- detectCores()
    if (length(available) != 1L || is.na(available) || available < 2) available <- 2L
    iCores <- max(1L, available - 1L)
    if (models * length(vDelta) <= 100) iCores <- min(2L, iCores)
    iCores <- min(iCores, models * length(vDelta), 65536)
  }
  nativekeep <- if (fixed) keep - 1 else -9999
  if (bParallelize) {
    Est <- funcEstimate_Eff_par(vY, mF, vDelta, dAlpha, nativekeep, dBeta,
                               bZellnerPrior, dG, iCores, dS0, dN0)
  } else {
    Est <- funcEstimate_Eff(vY, mF, vDelta, dAlpha, nativekeep, dBeta,
                           bZellnerPrior, dG, dS0, dN0)
  }
  if (forecast) {
    last <- length(vY)
    variance <- c(Est$mvdec[last, ])
    names(variance) <- c("vtotal", "vobs", "vcoeff", "vmod", "vtvp")
    lastForecast <- list(PointForecast = c(Est$vyhat[last, ]),
                         VarianceDecomposition = variance, bForecast = TRUE,
                         ForecastDate = if (is.null(dates)) NULL else dates[last])
    for (nm in names(Est)) if (is.matrix(Est[[nm]]) && nrow(Est[[nm]]) == last)
      Est[[nm]] <- Est[[nm]][-last, , drop = FALSE]
    vY <- vY[-last]
    mF <- mF[-last, , drop = FALSE]
    if (!is.null(dates)) dates <- dates[-last]
    Est$LastForecast <- lastForecast
  } else Est$LastForecast <- NULL
  if (!is.null(dates)) vY <- .dma_indexed(matrix(vY, ncol = 1L), dates)
  new("DMA", model = list(vDelta = vDelta, dAlpha = dAlpha,
      bZellnerPrior = bZellnerPrior, dG = dG, dBeta = dBeta, dS0 = dS0, dN0 = dN0,
      elapsedTime = Sys.time() - start, bParallelize = bParallelize,
      FixedVar = fixed, vKeep = keep, Call = formula),
      data = list(vY = vY, mF = mF, formula = formula, vDates = dates), Est = Est)
}
