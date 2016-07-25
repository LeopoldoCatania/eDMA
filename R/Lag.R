Lag <- function(mY, iK) {

  if (is.vector(mY)){
    out = c(rep(NA, iK), mY)[1 : length(mY)]
  } else {
    out = matrix(NA, nrow(mY), ncol(mY))
    for(i in 1:ncol(mY)){
      out[, i] = c(rep(NA, iK), mY[, i])[1 : nrow(mY)]
    }
    dimnames(out) = dimnames(mY)
  }

  return(out)

}
