na.omit_new <- function(mX){
  if (is.na(mX[nrow(mX), 1])) {
    mX_new = rbind(na.omit(mX), mX[nrow(mX), ])
  } else {
    mX_new = na.omit(mX)
  }
  return(mX_new)
}
