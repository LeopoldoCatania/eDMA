Lag <- function(mY, iK) {
  if (!is.numeric(iK) || length(iK) != 1L || !is.finite(iK) || iK < 0 || iK != floor(iK))
    stop("iK must be a nonnegative integer.")
  if (!(is.atomic(mY) && is.null(dim(mY))) && !is.matrix(mY))
    stop("mY must be a vector or a matrix.")
  n <- if (is.matrix(mY)) nrow(mY) else length(mY)
  out <- mY
  out[] <- NA
  if (iK < n) {
    to <- seq.int(iK + 1L, n)
    from <- seq_len(n - iK)
    if (is.matrix(mY)) out[to, ] <- mY[from, , drop = FALSE]
    else out[to] <- mY[from]
  }
  out
}
