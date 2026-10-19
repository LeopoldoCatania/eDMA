library(eDMA)
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop(label, call. = FALSE)
  checks <<- checks + 1L
  cat(sprintf("PASS %03d %s\n", checks, label))
}
close <- function(a, b, tolerance = 2e-9) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tolerance))
fails <- function(expr, pattern = NULL) {
  result <- tryCatch({force(expr); NULL}, error = identity)
  inherits(result, "error") && (is.null(pattern) || grepl(pattern, conditionMessage(result)))
}

# Independent small-model implementation: ordinary probability arithmetic,
# base R dt(), and the subtractive Kalman covariance formula.
reference <- function(y, X, deltas, alpha, beta, g = 4, S0 = 1, n0 = 2, required = integer()) {
  K <- ncol(X); T <- length(y); D <- length(deltas)
  models <- lapply(seq_len(2^K - 1L), function(mask) which(as.logical(intToBits(mask)[seq_len(K)])))
  models <- Filter(function(idx) all(required %in% idx), models)
  M <- length(models)
  states <- lapply(seq_len(M * D), function(h) {
    p <- length(models[[(h - 1L) %% M + 1L]])
    list(m = numeric(p), C = diag(g, p), S = S0)
  })
  a <- matrix(1/M, M, D); b <- rep(1/D, D); n <- n0
  output <- list(vyhat = numeric(T), vLpdfhat = numeric(T), mmhat = matrix(0,T,K),
                 mcovhat = matrix(0,T,K), mincpmt = matrix(0,T,K), mpmt = matrix(0,T,D),
                 mpmt_pred = matrix(0,T,D), mLpdf = matrix(0,T,D), mvdec = matrix(0,T,5),
                 vShat = numeric(T), vsize = numeric(T), vdeltahat = numeric(T),
                 vyhat_DMS = numeric(T), vLpdfhat_DMS = numeric(T),
                 vsize_DMS = numeric(T), vhighmp_DMS = numeric(T), vhighmpTop01_DMS = numeric(T))
  for (t in seq_len(T)) {
    ap <- a^alpha; ap <- sweep(ap, 2, colSums(ap), "/")
    bp <- b^alpha / sum(b^alpha)
    output$mpmt_pred[t, ] <- bp
    nu <- beta * n + 1
    f <- obs <- coeff <- density <- matrix(0, M, D)
    means <- variances <- array(0, c(K, M, D))
    for (d in seq_len(D)) for (j in seq_len(M)) {
      h <- (d-1L)*M+j; s <- states[[h]]; idx <- models[[j]]; x <- X[t, idx]
      R <- s$C / if (t == 1L) 1 else deltas[d]
      f[j,d] <- sum(x*s$m); coeff[j,d] <- drop(t(x)%*%R%*%x); obs[j,d] <- s$S
      q <- coeff[j,d]+s$S; e <- y[t]-f[j,d]
      density[j,d] <- dt(e/sqrt(q), df=nu)/sqrt(q)
      gain <- drop(R%*%x)/q
      s$m <- s$m+gain*e
      s$C <- R-tcrossprod(gain)*q
      s$S <- s$S+(s$S/nu)*(e^2/q-1)
      states[[h]] <- s
      means[idx,j,d] <- s$m; variances[idx,j,d] <- diag(s$C)
    }
    mu <- colSums(ap*f); overall <- sum(bp*mu)
    output$vyhat[t] <- overall
    vo <- sum(bp*colSums(ap*obs))*nu/(nu-2)
    vc <- sum(bp*colSums(ap*coeff))*nu/(nu-2)
    vm <- sum(bp*colSums(ap*sweep(f,2,mu,"-")^2))
    vt <- sum(bp*(mu-overall)^2)
    output$mvdec[t,] <- c(vo+vc+vm+vt,vo,vc,vm,vt)
    ds <- which.max(bp); js <- which.max(ap[,ds])
    output$vyhat_DMS[t] <- f[js,ds]
    output$vLpdfhat_DMS[t] <- log(density[js,ds])
    ld <- colSums(ap*density); normalizer <- sum(bp*ld)
    output$vLpdfhat[t] <- log(normalizer); output$mLpdf[t,] <- log(ld)
    a <- sweep(ap*density,2,ld,"/"); b <- bp*ld/normalizer
    joint <- sweep(a,2,b,"*"); marginal <- rowSums(joint)
    for (k in seq_len(K)) {
      output$mmhat[t,k] <- sum(joint*means[k,,])
      output$mcovhat[t,k] <- sum(joint*(variances[k,,]+(means[k,,]-output$mmhat[t,k])^2))
      output$mincpmt[t,k] <- sum(marginal[vapply(models,function(idx) k %in% idx,logical(1))])
    }
    output$mpmt[t,] <- b
    output$vShat[t] <- sum(as.vector(joint)*vapply(states,`[[`,numeric(1),"S"))
    output$vsize[t] <- sum(marginal*lengths(models))
    output$vdeltahat[t] <- sum(b*deltas)
    output$vsize_DMS[t] <- lengths(models)[which.max(marginal)]
    output$vhighmp_DMS[t] <- max(marginal)
    output$vhighmpTop01_DMS[t] <- sum(sort(marginal,decreasing=TRUE)[seq_len(ceiling(.1*M))])
    n <- nu
  }
  output
}
set.seed(181)
dat <- data.frame(y = rnorm(45), x = rnorm(45), z = rnorm(45), w = rnorm(45))
args <- list(formula = y~x+z+w, data = dat, vDelta = c(.8,.97), dAlpha = .83,
             dBeta = .96, dG = 4, bParallelize = FALSE)
fit <- do.call(DMA,args)
ref <- reference(dat$y,model.matrix(y~x+z+w,dat),args$vDelta,args$dAlpha,args$dBeta)
for (nm in names(ref)) check(close(fit@Est[[nm]],ref[[nm]]),paste("Independent numerical reference:",nm))
check(all(fit@Est$vhighmpTop01_DMS > fit@Est$vhighmp_DMS),"F11 top ten percent includes more than the largest model")
parargs <- args; parargs$bParallelize <- TRUE; parargs$iCores <- 2L
parallel <- suppressWarnings(do.call(DMA,parargs))
for (nm in names(ref)) check(close(fit@Est[[nm]],parallel@Est[[nm]],0),paste("Sequential/parallel equality:",nm))
fixedargs <- args; fixedargs$vKeep <- c(1,3)
fixed <- do.call(DMA,fixedargs)
fixedref <- reference(dat$y,model.matrix(y~x+z+w,dat),args$vDelta,args$dAlpha,args$dBeta,required=c(1,3))
check(close(fixed@Est$vyhat,fixedref$vyhat) && fixed@Est$iM==4,"Direct required-column enumeration")
check(all(abs(fixed@Est$mincpmt[,c(1,3)]-1)<1e-12),"Required columns have unit inclusion probability")

# Alter the current/future response without affecting any current prediction moments.
changed <- args; changed$data$y[25:45] <- changed$data$y[25:45]+20
altered <- do.call(DMA,changed)
check(close(fit@Est$vyhat[1:25],altered@Est$vyhat[1:25]),"F1 point forecast causality")
check(close(fit@Est$mvdec[1:25,],altered@Est$mvdec[1:25,]),"F5 predictive variance causality")
check(!close(fit@Est$mmhat[25,],altered@Est$mmhat[25,]),"F6 filtered means react to the current observation")
check(close(fit@Est$mvdec[,1],rowSums(fit@Est$mvdec[,-1,drop=FALSE])),"Predictive variance components add up")

futureargs <- args; futureargs$data <- rbind(dat,data.frame(y=NA_real_,x=6,z=-2,w=3))
future <- do.call(DMA,futureargs)
for (nm in names(ref)) check(close(fit@Est[[nm]],future@Est[[nm]],0),paste("F3 terminal prediction preserves history:",nm))
scoredargs <- futureargs; scoredargs$data$y[46] <- 8
scored <- do.call(DMA,scoredargs)
last <- getLastForecast(future)
check(close(last$PointForecast,scored@Est$vyhat[46]),"Future mean equals the pre-update forecast")
check(close(last$VarianceDecomposition,scored@Est$mvdec[46,]),"F4 variance belongs to the future predictor row")
check(!close(last$VarianceDecomposition,fit@Est$mvdec[45,]),"Future variance is not the prior row's variance")
native <- eDMA:::funcEstimate_Eff(c(dat$y,NA_real_), model.matrix(~x+z+w,futureargs$data), args$vDelta,.83,-9999,.96,FALSE,4)
check(all(is.na(native$veps[46,])) && is.na(native$vLpdfhat[46]),"Missing response is not scored")
check(close(native$mpmt[46,],native$mpmt[45,]),"Missing response does not update posterior probabilities")
check(close(native$mmhat[46,],native$mmhat[45,]),"Missing response does not update coefficients")

# Quadrature of the initial predictive mixture, independent of decomposition formulas.
X <- model.matrix(y~x+z+w,dat); xx <- X[1,]
subsets <- lapply(1:15,function(mask) which(as.logical(intToBits(mask)[1:4])))
q <- vapply(subsets,function(idx) 1+4*sum(xx[idx]^2),numeric(1))
nu <- .96*2+1
mixdensity <- function(y) vapply(y,function(value) mean(dt(value/sqrt(q),nu)/sqrt(q)),numeric(1))
check(abs(integrate(mixdensity,-Inf,Inf)$value-1)<1e-7,"F2 predictive density integrates to one")
second <- integrate(function(y) y^2*mixdensity(y),-Inf,Inf,rel.tol=1e-7)$value
check(abs(second-fit@Est$mvdec[1,1])<1e-5,"F5 variance agrees with Student-t density integration")

# Startup, dimensions, NA, time indexes, parameter boundaries and helper contracts.
single <- DMA(y~1,dat,vDelta=.95,bParallelize=FALSE)
check(nrow(coef(single,iBurnPeriod=0))==45,"Zero burn-in retains all observations")
check(identical(dim(coef(single,iBurnPeriod=44)),c(1L,1L)),"Single-row single-column extraction")
check(identical(dim(as.data.frame(single,which="mpmt",iBurnPeriod=5)),c(40L,1L)),"Single-delta extraction after burn-in")
check(close(BacktestDMA(single,0)[1,1],mean(single@Est$veps^2)),"Zero burn-in backtest")
warned <- FALSE
invisible(withCallingHandlers(BacktestDMA(single),warning=function(w) {warned <<- TRUE; invokeRestart("muffleWarning")}))
check(warned,"Default burn-in emits a warning")
check(fails(getLastForecast(single),"No future forecast"),"Useful error for absent future forecast")
check(fails(coef(single,iBurnPeriod=45)) && fails(coef(single,iBurnPeriod=-1)) && fails(coef(single,iBurnPeriod=.5)),"Burn-in boundaries validated")
check(fails(residuals(single,type="bad")) && fails(pred.like(single,type="bad")),"Residual and likelihood types validated")
check(close(sd(as.numeric(residuals(single,standardize=TRUE))),1),"Residual standardization uses global sample SD")
one <- DMA(y~1,data.frame(y=0),bParallelize=FALSE,vDelta=1)
check(is.finite(one@Est$mvdec[1,1]) && one@Est$vyhat[1]==0 && one@Est$mmhat[1]==0,"F8 zero first response and one observation are defined")
zero <- DMA(y~0+x,data.frame(y=c(0,1,0),x=c(0,1,2)),vDelta=1,bParallelize=FALSE)
check(all(is.finite(zero@Est$mvdec)) && all(is.finite(zero@Est$vLpdfhat)),"F8 zero initial design is supported")
smallnu <- DMA(y~1,dat,vDelta=1,dBeta=.1,bParallelize=FALSE)
check(is.infinite(smallnu@Est$vtotal[1]),"Infinite predictive variance for degrees of freedom at most two")
large <- dat; large$y[1] <- 1e50
largefit <- DMA(y~x+z+w,large,bParallelize=FALSE)
check(all(is.finite(largefit@Est$vLpdfhat)) && all(is.finite(largefit@Est$mpmt)),"F9 extreme observation does not underflow probability normalization")
for (keep in list(character(),"unknown",c("x","unknown"),c(1,1),0,5,.5,NA_real_,c("KS","x")))
  check(fails(DMA(y~x,dat,vKeep=keep,bParallelize=FALSE)),paste("Invalid vKeep rejected:",paste(keep,collapse=",")))
for (entry in list(list(vDelta=numeric()),list(vDelta=0),list(vDelta=NA_real_),list(dAlpha=1.1),list(dBeta=0),
                   list(dG=0),list(dS0=0),list(dN0=0),list(iCores=0),list(iCores=1.2),list(bParallelize=NA),list(bZellnerPrior=1))) {
  bad <- utils::modifyList(args,entry)
  check(fails(do.call(DMA,bad)),paste("Invalid parameter rejected:",names(entry)))
}
check(fails(DMA(y~0,dat)),"Zero-column design rejected")
check(fails(DMA(y~x,data.frame(y=c(NA,NA),x=1:2))),"All-missing response rejected")
check(fails(DMA(y~x,data.frame(y=factor(c("a","b")),x=1:2))),"Nonnumeric response rejected")
check(fails(DMA(y~x,data.frame(y=c(1,Inf),x=1:2))),"Nonfinite observed response rejected")
check(fails(DMA(y~x,data.frame(y=c(1,NA),x=c(1,NA)))),"Missing future predictor rejected")
check(fails(DMA(y~x,data.frame(y=1:5,x=rep(1,5)),bZellnerPrior=TRUE)),"Rank-deficient Zellner prior rejected before workers")
zargs <- args; zargs$bZellnerPrior <- TRUE
zfit <- do.call(DMA,zargs); zargs$data <- futureargs$data; zfuture <- do.call(DMA,zargs)
check(close(zfit@Est$vyhat,zfuture@Est$vyhat,0),"Zellner prior excludes the forecast-only row")
wide <- as.data.frame(matrix(rnorm(10*31),10,31)); names(wide) <- c("y",paste0("x",1:30))
widefit <- DMA(y~.,wide,vKeep="KS",vDelta=1,bParallelize=FALSE)
check(widefit@Est$iM==1,"Kitchen sink with 30 predictors avoids exponential enumeration")
check(fails(DMA(y~.,wide,bParallelize=FALSE)),"Oversized optional model space rejected")
check(fails(eDMA:::funcEstimate_Eff(1:3,matrix(1,2,1),1,.9,-9999,1)),"Native design dimension validation")

quarterly <- ts(as.matrix(dat[,1:2]),start=c(2000,1),frequency=4)
tfit <- DMA(y~x,quarterly,vDelta=1,bParallelize=FALSE)
check(close(tfit@data$vDates,time(quarterly)) && inherits(coef(tfit),"zoo"),"Quarterly ts retains numeric time index")
for (dates in list(as.Date("2000-01-01")+0:45,as.POSIXct("2000-01-01",tz="UTC")+(0:45)*3600,2000+(0:45)/4)) {
  indexed <- zoo::zoo(as.matrix(futureargs$data),dates)
  if (inherits(dates,"POSIXct")) indexed <- xts::as.xts(indexed)
  dated <- DMA(y~x+z+w,indexed,vDelta=1,bParallelize=FALSE)
  check(identical(dated@data$vDates,dates[1:45]) && identical(getLastForecast(dated)$ForecastDate,dates[46]),"Observed and future time indexes stay aligned")
}
lagged <- DMA(y~Lag(x,2),quarterly,vDelta=1,bParallelize=FALSE)
check(close(lagged@data$vDates,time(quarterly)[-(1:2)]),"Lag omissions preserve time alignment")
factors <- data.frame(y=rnorm(10), f=factor(rep(c("a","b"),5)))
factorfit <- DMA(y~f,factors,bParallelize=FALSE,vKeep="fb")
check(factorfit@Est$iM==2,"Factor-expanded columns define the model space")

graphicsfile <- tempfile(fileext=".pdf")
grDevices::pdf(graphicsfile)
device <- grDevices::dev.cur(); params <- par("mfrow")
plot(single,which="mpmt",iBurnPeriod=44)
plot(future,which="mvdec")
plot(tfit,which="mtheta")
plot(smallnu,which="vtotal")
plot(widefit,which="mtheta")
check(grDevices::dev.cur()==device && identical(par("mfrow"),params),"Plot preserves device and graphical parameters")
grDevices::dev.off(); unlink(graphicsfile)
check(fails(plot(single,which="unknown")),"Invalid plot output rejected")
check(close(Lag(1:3,1),c(NA,1,2)) && identical(Lag(numeric(),2),numeric()),"Lag ordinary and empty vectors")
m <- matrix(1:4,2,dimnames=list(c("a","b"),c("x","y")))
check(identical(dimnames(Lag(m,1)),dimnames(m)) && identical(Lag(m,0),m),"Lag dimensions and names preserved")
check(fails(Lag(1:3,-1)) && fails(Lag(1:3,.5)),"Invalid lags rejected")
check(length(PowerSet(4))==16 && max(unlist(PowerSet(4)))==3 && length(PowerSet(0))==1,"PowerSet endpoint and empty set")
check(fails(PowerSet(-1)) && fails(PowerSet(.5)) && fails(PowerSet(21)),"PowerSet allocation guards")
Xsim <- cbind(1,seq_len(20))
sim <- SimulateDLM(20,Xsim,c(1,2),matrix(0,2,2),0,1)
check(identical(dim(sim$mBeta),c(2L,20L)) && close(sim$vY,1+2*seq_len(20)),"Zero covariance simulation and K by T coefficient shape")
set.seed(99); s1 <- SimulateDLM(20,Xsim,c(0,0),matrix(1,2,2),1,1)
set.seed(99); s2 <- SimulateDLM(20,Xsim,c(0,0),matrix(1,2,2),1,1)
check(identical(s1,s2) && close(s1$mBeta[1,],s1$mBeta[2,]),"Rank-deficient PSD simulation and R seed reproducibility")
check(fails(SimulateDLM(.5,Xsim,c(0,0),diag(2),1,1)) && fails(SimulateDLM(20,Xsim,c(0,0),-diag(2),1,1)),"Simulation size and PSD guards")
check(fails(SimulateDLM(20,Xsim,c(0,0),diag(2),-1,1)),"Simulation observation variance guard")
# Exercise numerical failures inside the component loop in both execution modes.
for (parallelize in c(FALSE, TRUE)) {
  overflow <- dat; overflow$y[1] <- 1e308
  check(suppressWarnings(fails(DMA(y~x+z+w, overflow, bParallelize=parallelize, iCores=2), "Time 1, model")),
        paste("Numerical failure safely reported on main thread; parallel =",parallelize))
}
check(fails(residuals(one, standardize=TRUE), "undefined"),"Undefined residual SD rejected")
undiscounted <- args; undiscounted$dAlpha <- 1
ufit <- do.call(DMA,undiscounted)
uref <- reference(dat$y,model.matrix(y~x+z+w,dat),args$vDelta,1,args$dBeta)
check(close(ufit@Est$vyhat,uref$vyhat) && close(ufit@Est$vLpdfhat,uref$vLpdfhat),"Alpha one gives undiscounted Bayesian weights")
cat(sprintf("\nAll %d regression checks passed.\n",checks))
