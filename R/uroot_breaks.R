kpss <- function (y,
                  x=NULL,
                  lags = c("short", "long", "nil"),
                  use=c("nw","ba")) {
  #==Bandwidth
  # nw : kernel = c("ba", "pa", "qs")
  # and: kernel = c("ba", "pa", "qs", "th", "tr")
  y0=y
  y = as.matrix(y)
  if(ncol(y)>1) {stop("\nError: Only one series is allowed.")}
  lags <- match.arg(lags)
  const=rep(1,nrow(as.matrix(y)))

  if (is.null(x)) {datmat=data.frame(y,const)} else {
    datmat=data.frame(y,x)
  }


  t <- nrow(datmat)

  e=as.matrix(residuals(lm(datmat)))
  S = cumsum(e)



  if (!(is.null(use))) {

    if (length(use)==2L) {

      lmax <- trunc(cointReg::getBandwidth(y,bandwidth = use[1],
                                           kernel = use[2]))
      if (use[2]=="ba") {

        teststat <- 1/nrow(e)^2 * t(S)%*%S/Bartlett_uni(e,lmax)
      }  else if (use[2]=="qs") {
        teststat <- 1/nrow(e)^2 * t(S)%*%S/QS_uni(e,lmax)
      }  else if (use[2]=="pa") {
        teststat <- 1/nrow(e)^2 * t(S)%*%S/Parzen_uni(e,lmax)
      }

    } else if (is.numeric(use)) {

      lmax=as.integer(use)

      if (lmax < 0) {
        warning("\nuse.lag has to be positive and integer; lags='short' used.")
        lmax <- trunc(4 * (t/100)^(2/9))
      }
      teststat <- 1/nrow(e)^2 * t(S)%*%S/Bartlett_uni(e,lmax)
    }

  }  else if (lags == "short") {

    lmax <- trunc(4 * (t/100)^(2/9))

    teststat <- 1/nrow(e)^2 * t(S)%*%S/Bartlett_uni(e,lmax)

  } else if (lags == "long") {

    lmax <- trunc(12 * (t/100)^(2/9))

    teststat <- 1/nrow(e)^2 * t(S)%*%S/Bartlett_uni(e,lmax)
  }   else if (lags == "nil") {
    lmax <- 0
    teststat <- 1/nrow(e)^2 * t(S)%*%S/Bartlett_uni(e,lmax)

  } else {

    lmax <- trunc(4 * (t/100)^(2/9))

    teststat <- 1/nrow(e)^2 * t(S)%*%S/Bartlett_uni(e,lmax) }



  if (ncol(datmat) == 2) { #intercept
    cval <- cbind(0.347, 0.463, 0.574, 0.739)
    colnames(cval) <- c("10%", "5%", "2.5%", "1%")
    rownames(cval) <- "critical values"

  }  else if (ncol(datmat) > 2) { #intercept and trend
    cval <- cbind(0.119, 0.146, 0.176, 0.216)
    colnames(cval) <- c("10%", "5%", "2.5%", "1%")
    rownames(cval) <- "critical values"

  }

  return(list(y = y,
              lag = as.integer(lmax),
              teststat = as.numeric(teststat),
              cval = cval,
              resid = e))
}



pp <- function(y,
               type=c("none","const","trend"),
               d=NULL,
               lags=c("short","long","nill"),
               use=c("nw","ba")) {
  #==Bandwidth
  # nw : kernel = c("ba", "pa", "qs")
  # and: kernel = c("ba", "pa", "qs", "th", "tr")
  y0=y
  y = as.matrix(y)
  if(ncol(y)>1) {stop("\nError: Only one series is allowed.")}
  t=nrow(y)
  lags <- match.arg(lags)
  model <- match.arg(type)


  if(type=="trend") {

    datmat=.datmat4ADF(y,x=seq(t)/t,p=0,eq=2)

  } else {

    datmat=.datmat4ADF(y,p=0,eq=2)

  }

  if (is.null(d)) {datmat=datmat
  } else {
    d=as.matrix(d)
    datmat=data.frame(datmat,embed(d,2)[,seq(ncol(d))])}

  if (type=="none") {
    output=lm(dy~y.L1-1,data=datmat)
    COEF=coef(summary(output))
    tau=COEF[1,3] #t-ratio
    se=COEF[1,2]
    k=length(coef(output))

  } else if (type=="const") {

    output=lm(data=datmat)
    COEF=coef(summary(output))
    tau=COEF[2,3]#t-ratio
    se=COEF[2,2]
    k=length(coef(output))

  } else  {
    output=lm(data=datmat)
    COEF=coef(summary(output))
    tau=COEF[2,3]#t-ratio
    se=COEF[2,2]
    k=length(coef(output))
  }

  names(summary(output))
  summary(output)$terms


  e=as.matrix(resid(output))
  s2=summary(output)$sigma^2

  # Residual variance as in Eviews
  t = nrow(datmat)
  g0 = (t - k) * s2 * 1/t
  e1 = e;


  if (!(is.null(use))) {

    if (length(use)==2L) {

      lmax <- trunc(cointReg::getBandwidth(y,bandwidth = use[1],
                                           kernel = use[2]))
      if (use[2]=="ba") { lrv=Bartlett_uni(e,v=lmax)
      } else if (use[2]=="qs") {lrv=QS_uni(e,v=lmax)
      } else if (use[2]=="pa") {lrv=Parzen_uni(e,v=lmax)}

    } else if (is.numeric(use)) {

      lmax=as.integer(use)
      if (lmax < 0) {
        warning("\nuse.lag has to be positive and integer; lags='short' used.")
        lmax <- trunc(4 * (t/100)^(2/9))
      }
      lrv=Bartlett_uni(e,v=lmax)
    }

  } else if (lags == "short") {
    lmax <- trunc(4 * (t/100)^(2/9))
    lrv=Bartlett_uni(e,v=lmax)
  } else if (lags == "long") {
    lmax <- trunc(12 * (t/100)^(2/9))
    lrv=Bartlett_uni(e,v=lmax)
  }   else if (lags == "nil") {
    lmax <- 0
    lrv=Bartlett_uni(e,v=lmax)
  } else {
    lmax <- trunc(4 * (t/100)^0.25)
    lrv=Bartlett_uni(e,v=lmax)

  }

  #Eviews
  Zt = tau * sqrt(g0/lrv) - t * (lrv - g0) * se/(2 * sqrt(lrv) * sqrt(s2))

  #Stata
  lambda = 1/2 * (lrv - g0)
  mbaryy = (t * se[1]/sqrt(s2))^2

  if (type == "none") {
    Za = t * COEF[1,1] - lambda * mbaryy
  } else {Za = t * COEF[2,1] - lambda * mbaryy}


  #Critical values for Zt
  if (type == "none") {
    cv_Zt =matrix(c(-2.66,-1.95,-1.60,      # T = 25
                    -2.62,-1.95,-1.61,      #  T = 50
                    -2.60,-1.95,-1.61,      #  T = 100
                    -2.58,-1.95,-1.62,      #  T = 250
                    -2.58,-1.95,-1.62,      #  T = 500
                    -2.58,-1.95,-1.62,      #  T = 750
                    -2.58,-1.95,-1.62),7,3,byrow = TRUE)     #  T = inf

    cv_Za =matrix(c(-11.87, -7.32, -5.32,      # T = 25
                    -12.82, -7.69, -5.52,      # T = 50
                    -13.30, -7.88, -5.61,      # T = 100
                    -13.59, -7.99, -5.67,      # T = 250
                    -13.69, -8.03, -5.69,      # T = 500
                    -13.72, -8.04, -5.70,      # T = 750
                    -13.78, -8.07, -5.71),7,3,byrow = TRUE)     #  T = inf



  } else if (type == "const"){
    cv_Zt = matrix(c(-3.75,-2.99,-2.64,      # T = 25
                     -3.59,-2.93, -2.60,      # T = 50
                     -3.50, -2.90, -2.59,      # T = 100
                     -3.46, -2.88, -2.58,      # T = 250
                     -3.44, -2.87, -2.57,      # T = 500
                     -3.43, -2.87, -2.57,      # T = 750
                     -3.42, -2.86, -2.57),7,3,byrow = TRUE)     #  T = inf

    cv_Za = matrix(c(-17.22, -12.47, -10.23,      # T = 25
                     -18.94, -13.29,-10.75,      # T = 50
                     -19.81, -13.70, -11.00,      # T = 100
                     -20.32, -13.95, -11.16,      # T = 250
                     -20.50, -14.03, -11.21,      # T = 500
                     -20.55, -14.06, -11.23,      # T = 750
                     -20.67, -14.11, -11.26),7,3,byrow = TRUE)     #  T = inf


  } else if (type == "trend") {
    cv_Zt = matrix(c(-4.38, -3.60, -3.24,      # T = 25
                     -4.15, -3.50, -3.18,      # T = 50
                     -4.04, -3.45, -3.15,      # T = 100
                     -3.98, -3.42, -3.13,      # T = 250
                     -3.97, -3.42, -3.13,      # T = 500
                     -3.96, -3.41, -3.13,      # T = 750
                     -3.96, -3.41, -3.13),7,3,byrow = TRUE)     #  T = inf

    cv_Za = matrix(c(-22.51, -17.89, -15.55,      # T = 25
                     -22.65, -19.70, -16.84,      # T = 50
                     -27.33, -20.64, -17.50,      # T = 100
                     -28.42, -21.25, -17.93,      # T = 250
                     -28.84, -21.47, -18.08,      # T = 500
                     -29.00, -21.56, -18.13,      # T = 750
                     -29.47, -21.78, -18.28),7,3,byrow = TRUE)     #  T = inf

  }

  if (t <= 25) {
    cvZt = cv_Zt[1,,drop=FALSE]
    cvZa = cv_Za[1,,drop=FALSE]
  } else if (25 < t & t <= 50){
    cvZt = cv_Zt[2,,drop=FALSE]
    cvZa = cv_Za[2,,drop=FALSE]
  } else if (50 < t & t <= 100) {
    cvZt = cv_Zt[3,,drop=FALSE]
    cvZa = cv_Za[3,,drop=FALSE]
  } else if (100< t & t <= 250) {
    cvZt = cv_Zt[4,]
    cvZa = cv_Za[4,,drop=FALSE]
  } else if (250< t & t <= 500) {
    cvZt = cv_Zt[5,,drop=FALSE]
    cvZa = cv_Za[5,,drop=FALSE]
  } else if (500< t &  t <= 750) {
    cvZt = cv_Zt[6,,drop=FALSE]
    cvZa = cv_Za[6,,drop=FALSE]
  } else if (750 < t) {
    cvZt = cv_Zt[7,,drop=FALSE]
    cvZa = cv_Za[7,,drop=FALSE]
  }

  colnames(cvZt)=colnames(cvZa)=c("1%","5%","10%")
  rownames(cvZt)=rownames(cvZa)="critical values"

  return(list(Zt=Zt,
              cvZt=cvZt[1,],
              Za=Za,
              cvZa=cvZa[1,],
              lag = as.integer(lmax)))

}

# Zivot-Andrews unit-root test with 1 structural break
ZA_1br <- function(y,
                   model=c("intercept", "trend", "both"),
                   outlier=1,
                   pmax=8,
                   ic=c("AIC","BIC"),
                   fixed=FALSE,
                   trim=0.1,
                   eq=1,
                   season=FALSE) {

  ic <- match.arg(ic)
  model <- match.arg(model)

  if(isTRUE(season)) {

    if (is.ts(y)  | timeSeries::is.timeSeries(y)) {

      SD=forecast::seasonaldummy(as.ts(y))

    } else {stop("\nError: When season=TRUE, time series must have time series structure.")}

  }
  y0=y
  y = as.matrix(y)

  if(ncol(y)>1) {stop("\nError: Only one series is allowed.")}


t0=proc.time();if (outlier==1) {# Innovational outlier model

    t=nrow(y)
    trend=seq(t)/t
    T1 = max((3 + pmax),ceiling(trim * t))
    T2 = min((t - 3 -pmax),floor((1 - trim) * t))
    if (T1 < (pmax + 2)) {T1 = pmax + 3}
    idx = T1:T2

   if (model == "intercept") {

   roll <- function(z) {
      du1 = c(rep(0,z), rep(1, (t - z)))
      br1 = c(rep(0,z),1, rep(0,t-z-1))

      if (isTRUE(season)) {
        x=cbind(du1,br1,SD)
      } else {x=cbind(du1,br1)}

      if(isTRUE(fixed)) {
      rollmat=.datmat4ADF(y,x,p=pmax,eq=eq)
      } else {rollmat=.plag(y,x,pmax,ic,eq=eq)$data}

      COEF = coef(summary(lm(rollmat)))

      if (eq==1) {
          (COEF[2, 1] - 1)/COEF[2, 2]
      } else {
        COEF[2, 3]
        }
    }

    roll.stat <- sapply(idx, roll)
    cval = cbind(-5.34, -4.8, -4.58)
    bpoint = which.min(roll.stat)
    du1 = c(rep(0, idx[bpoint]), rep(1, (t - idx[bpoint])))
    br1 = c(rep(0,idx[bpoint]),1, rep(0,t-idx[bpoint]-1))

    if (isTRUE(season)) {
      x=cbind(du1,br1,SD)
    } else {x=cbind(du1,br1)}

    if(isTRUE(fixed)) {
    testmat= .datmat4ADF(y,x,p=pmax,eq=eq)
    p=pmax
      } else {
    out=.plag(y,x,pmax,ic,eq=eq)
    p=out$p
    testmat <- out$data
    }

    test.reg=lm(testmat)

  } else if (model == "trend") {

    roll = function(z) {
      dt1 = c(rep(0, z), 1:(t - z))/t
      br1 = c(rep(0,z),1, rep(0,t-z-1))

      if (isTRUE(season)) {
        x=cbind(trend,dt1,br1,SD)
      } else {x=cbind(trend,dt1,br1)}

      if(isTRUE(fixed)) {
        rollmat=.datmat4ADF(y,x,p=pmax,eq=eq)
      } else {rollmat=.plag(y,x,pmax,ic,eq=eq)$data}

      COEF = coef(summary(lm(rollmat)))
      if (eq==1) {
        (COEF[2, 1] - 1)/COEF[2, 2]
      } else {(COEF[2, 3])}
    }

    roll.stat <- sapply(idx, roll)
    cval = cbind(-4.93, -4.42, -4.11)
    bpoint = which.min(roll.stat)
    dt1 = c(rep(0, idx[bpoint]), 1:(t - idx[bpoint]))
    br1 = c(rep(0,idx[bpoint]),1, rep(0,t-idx[bpoint]-1))/t

    if (isTRUE(season)) {
      x=cbind(trend,dt1,br1,SD)
    } else {x=cbind(trend,dt1,br1)}

    if(isTRUE(fixed)) {
      testmat= .datmat4ADF(y,x,p=pmax,eq=eq)
      p=pmax
    } else {
      out=.plag(y,x,pmax,ic,eq=eq)
      p=out$p
      testmat <- out$data
    }
    test.reg=lm(testmat)

  }  else if (model == "both") {

    roll <- function(z) {
      du1 = c(rep(0, z), rep(1, (t - z)))
      dt1 = c(rep(0, z), 1:(t - z))/t
      br1 = c(rep(0,z),1, rep(0,t-z-1))

      if (isTRUE(season)) {
        x=cbind(trend,du1,dt1,br1,SD)
      } else {x=cbind(trend,du1,dt1,br1)}

      if(isTRUE(fixed)) {
        rollmat=.datmat4ADF(y,x,p=pmax,eq=eq)
      } else {rollmat=.plag(y,x,pmax,ic,eq=eq)$data}

      COEF <- coef(summary(lm(rollmat)))
      if (eq==1) {
        (COEF[2, 1] - 1)/COEF[2, 2]
      } else {
        COEF[2, 3] }
    }

    roll.stat <- sapply(idx, roll)
    cval <- cbind(-5.57, -5.08, -4.82)
    bpoint <- which.min(roll.stat)
    du1 = c(rep(0, idx[bpoint]), rep(1, (t - idx[bpoint])))
    dt1 = c(rep(0, idx[bpoint]), 1:(t - idx[bpoint]))/t
    br1 = c(rep(0,idx[bpoint]),1, rep(0,t-idx[bpoint]-1))

    if (isTRUE(season)) {
      x=cbind(trend,du1,dt1,br1,SD)
    } else {x=cbind(trend,du1,dt1,br1)}

    if(isTRUE(fixed)) {
      testmat= .datmat4ADF(y,x,p=pmax,eq=eq)
      p=pmax
    } else {
      out=.plag(y,x,pmax,ic,eq=eq)
      p=out$p
      testmat <- out$data
    }
    test.reg=lm(testmat)

  }
#========================================================
} else if (outlier==2) {# Additive outlier model

  t = nrow(y)
  trend=seq(t)/(t)

  T1 = max((3 + pmax),ceiling(trim * t))
  T2 = min((t - 3 -pmax),floor((1 - trim) * t))


  if (T1 < (pmax + 2)) {T1 = pmax + 3}

  idx = T1:T2

  if (model == "intercept") {

    roll <- function(z) {#z=T1

      du1 <- c(rep(0, z), rep(1, (t - z)))

      if (isTRUE(season)) {
        datmat <- data.frame(y,du1,SD)
      } else { datmat <- data.frame(y,du1) }

      e1=as.matrix(resid(lm(datmat)))
      br1 = c(rep(0,z),1, rep(0,t-z-1))

      if(isTRUE(fixed)) {
        rollmat=.datmat4ADF(y=e1,x=br1,p=pmax,eq=eq)
      } else {rollmat<-.plag(y=e1,x=br1,pmax,ic,eq)$data}


      COEF <- coef(summary(lm(rollmat)))

      if (eq==1) {
        (COEF[2, 1] - 1)/COEF[2, 2]
      } else {

        COEF[2, 3] }
    }

    roll.stat <- sapply(idx, roll)
    cval <- cbind(-5.34, -4.8, -4.58)
    bpoint <- which.min(roll.stat)
    du1 <- c(rep(0, idx[bpoint]), rep(1, (t - idx[bpoint])))

    if (isTRUE(season)) {
      datmat <- data.frame(y,du1,SD)
    } else { datmat <- data.frame(y,du1) }

    e1=as.matrix(resid(lm(datmat)))
    br1 = c(rep(0,idx[bpoint]),1, rep(0,t-idx[bpoint]-1))

    if(isTRUE(fixed)) {
      testmat= .datmat4ADF(y=e1,x=br1,p=pmax,eq=eq)
      p=pmax
    } else {
    out=.plag(y=e1,x=br1,pmax,ic,eq)
    p=out$p
    testmat <- out$data
    }

    test.reg <- lm(testmat)

  } else if (model == "trend") {

    roll <- function(z) {
      dt1 <- c(rep(0, z), 1:(t - z))/t

      if (isTRUE(season)) {
        datmat <- data.frame(y,trend,dt1,SD)
      } else { datmat <- data.frame(y,trend,dt1) }


      e1=as.matrix(resid(lm(datmat)))

      br1 = c(rep(0,z),1, rep(0,t-z-1))

      if(isTRUE(fixed)) {
        rollmat=.datmat4ADF(y=e1,x=br1,p=pmax,eq=eq)
      } else {rollmat<-.plag(y=e1,x=br1,pmax,ic,eq)$data}

      COEF <- coef(summary(lm(rollmat)))
      if (eq==1) {
        (COEF[2, 1] - 1)/COEF[2, 2]
      } else {
        COEF[2, 3] }
    }

    roll.stat <- sapply(idx, roll)
    cval <- cbind(-4.93, -4.42, -4.11)
    bpoint <- which.min(roll.stat)
    dt1 <- c(rep(0, idx[bpoint]), 1:(t - idx[bpoint]))/t

    if (isTRUE(season)) {
      datmat <- data.frame(y,trend,dt1,SD)
    } else { datmat <- data.frame(y,trend,dt1) }

    e1=as.matrix(resid(lm(datmat)))
    br1 = c(rep(0,idx[bpoint]),1, rep(0,t-idx[bpoint]-1))

    if(isTRUE(fixed)) {
      testmat= .datmat4ADF(y=e1,x=br1,p=pmax,eq=eq)
      p=pmax
    } else {
      out=.plag(y=e1,x=br1,pmax,ic,eq)
      p=out$p
      testmat <- out$data
    }

    test.reg <- lm(testmat)

  }  else if (model == "both") {

    roll <- function(z) {
      du1 <- c(rep(0, z), rep(1, (t - z)))
      dt1 <- c(rep(0, z), 1:(t - z))/t

      if (isTRUE(season)) {
        datmat <- data.frame(y,trend, du1, dt1, SD)
      } else { datmat <- data.frame(y,trend, du1, dt1) }


      e1=as.matrix(resid(lm(datmat)))

      br1 = c(rep(0,z),1, rep(0,t-z-1))

      if(isTRUE(fixed)) {
        rollmat=.datmat4ADF(y=e1,x=br1,p=pmax,eq=eq)
      } else {rollmat<-.plag(y=e1,x=br1,pmax,ic,eq)$data}

      COEF <- coef(summary(lm(rollmat)))
      if (eq==1) {
        (COEF[2, 1] - 1)/COEF[2, 2]
      } else {
        COEF[2, 3]}
    }

    roll.stat <- sapply(idx, roll)
    cval <- cbind(-5.57, -5.08, -4.82)
    bpoint <- which.min(roll.stat)
    du1 <- c(rep(0, idx[bpoint]), rep(1, (t - idx[bpoint])))
    dt1 <- c(rep(0, idx[bpoint]), 1:(t - idx[bpoint]))/t

    if (isTRUE(season)) {
      datmat <- data.frame(y,trend, du1, dt1, SD)
    } else { datmat <- data.frame(y,trend, du1, dt1) }

    e1=as.matrix(resid(lm(datmat)))
    br1 = c(rep(0,idx[bpoint]),1, rep(0,t-idx[bpoint]-1))

    if(isTRUE(fixed)) {
      testmat= .datmat4ADF(y=e1,x=br1,p=pmax,eq=eq)
      p=pmax
    } else {
      out=.plag(y=e1,x=br1,pmax,ic,eq)
      p=out$p
      testmat <- out$data
    }

    test.reg <- lm(testmat)
  }

};t1=proc.time()


  colnames(cval)=c("1%","5%","10%")
  rownames(cval) <- "critical values"
  teststat <- roll.stat[bpoint]
  return(list(teststat = teststat,
              cval = cval,
              p = p,
              bpoint = bpoint,
              tstats = roll.stat,
              testreg = test.reg,
              timeElapse=t1-t0))
}


ZA_2br <- function(y,
                   model = c("intercept","both"),
                   pmax = 8,
                   ic=c("AIC","BIC"),
                   fixed=TRUE,
                   trim=0.1,
                   eq=1,
                   trace=TRUE,
                   season=FALSE) {

  ic <- match.arg(ic)
  model <- match.arg(model)

  if(isTRUE(season)) {

    if (is.ts(y) | timeSeries::is.timeSeries(y)) {

    SD=forecast::seasonaldummy(as.ts(y))

    } else {stop("\nError: When season=TRUE, time series must have time series structure.")}

  }
  y0=y
  y = as.matrix(y)

  if(ncol(y)>1) {stop("\nError: Only one series is allowed.")}

  t = nrow(y)
  trend=seq(t)/(t)

  #T1 = round(trim * t)
  #T2 = round((1 - trim) * t)
  T1 = rbind(max((3 + pmax), ceiling(trim * t)))
  T2 = min(rbind(t - 3 -pmax),floor((1 - trim) * t))

  tb1 = T1

  if (model == "intercept") {


    if (T1 < pmax+2) {T1 = pmax + 3}


    tb2 = tb1 + 2

    idx1=tb1:T2
    idx2=tb2:T2

    roll1 <- function(z) { #

      if(isTRUE(trace)) {(T2-z+1)}

      du1 = c(rep(0,z),rep(1,t-z)) #tb1

      if(isTRUE(fixed)) {

        if(isTRUE(season)) {
          datmat=.datmat4ADF(y,x=cbind(du1,SD),p=pmax,eq=eq)
          p=pmax
        } else {
        datmat=.datmat4ADF(y,x=du1,p=pmax,eq=eq)
        p=pmax}



      } else {

        if(isTRUE(season)) {
          out=.plag(y,x=cbind(du1,SD),pmax,ic,eq)
          datmat = out$data
          p=out$p

        } else {
          out=.plag(y,x=du1,pmax,ic,eq)
          datmat = out$data
          p=out$p
        }

        }

      roll2 <- function(z) {#
        du2 = c(rep(0,z),rep(1,t-z))
        du2=embed(du2,p+2)[,1,drop=F]
        rollmat <- data.frame(datmat,du2)
        COEF <- coef(summary(lm(rollmat)))
        if (eq==1) {
          (COEF[2, 1] - 1)/COEF[2, 2]
        } else {
          COEF[2, 3]
        }

      }
      sapply(idx2, roll2)
    }

t0=proc.time();rollStat <- sapply(idx1, roll1); t1=proc.time()


  } else if (model == "both") {

    tb2 = tb1 + 3

    idx1=tb1:T2
    idx2=tb2:T2


    roll1 <- function(z) { #z=68

      if(isTRUE(trace)) {(T2-z+1)}

      du1 = c(rep(0,z),rep(1,t-z))
      dt1 = c(rep(0,z),seq(t-z))/t
      trend=seq(nrow(datmat))/nrow(datmat)

      if(isTRUE(fixed)) {

        if(isTRUE(season)) {
          datmat=.datmat4ADF(y,x=cbind(du1,dt1,trend,SD),p=pmax,eq=eq)
          p=pmax
        } else {
          datmat=.datmat4ADF(y,x=cbind(du1,dt1,trend),p=pmax,eq=eq)
          p=pmax
        }

      } else {

        if(isTRUE(season)) {
          out=.plag(y,x=cbind(du1,dt1,trend,SD),pmax,ic,eq)
          datmat = out$data
          p=out$p
        } else {
          out=.plag(y,x=cbind(du1,dt1,trend),pmax,ic,eq)
          datmat = out$data
          p=out$p}

      }

      roll2 <- function(z) {
        du2 = c(rep(0,z),rep(1,t-z))
        dt2 = c(rep(0,z),seq(t-z))
        du2=embed(du2,p+2)[,1,drop=F]
        dt2=embed(dt2,p+2)[,1,drop=F]/nrow(datmat)
        rollmat <- data.frame(datmat,du2,dt2)
        COEF <- coef(summary(lm(rollmat)))

        if (eq==1) {
          (COEF[2, 1] - 1)/COEF[2, 2]
        } else {COEF[2, 3]}

      }

      sapply(idx2, roll2)
    }

    t0=proc.time();rollStat <- sapply(idx1, roll1); t1=proc.time()

  }

  # Critical Values (see, Narayan & Popp, 2010, Table 3)
  if (model == "intercept") {
    if (t <= 50) {
      cv = cbind(-5.259, -4.514,-4.143)
    } else if (50 < t & t <= 200) {

      cv = cbind(-4.958,-4.316,-3.980)

    } else if (200 < t & t <= 400) {

      cv = cbind(-4.731,-4.136,-3.825)

    } else if (400 < t) {

      cv = cbind(-4.672,-4.081,-3.772)

    }
  } else if (model == "both") {
    if (t <= 50) {
      cv =  cbind(-5.949,-5.181,-4.789)
    } else if (50 < t & t <= 20){
      cv =  cbind(-5.576,-4.937,-4.596)
    }   else if (200 < t & t <= 400) {
      cv =  cbind(-5.318,-4.741,-4.430)
    } else if (400 < t) {
      cv =  cbind(-5.287,-4.692,-4.396)
    }
  }
  colnames(cv)=c("1%","5%","10%")
  rownames(cv) <- "critical values"
  R=which.min(apply(rollStat,1,min)) #row number, tb2
  C=which.min(rollStat[R,]) #column number, tb1
  bp1=idx1[C]
  bp2=idx2[R]

  return(list(y = y,
              teststat = min(rollStat),
              cval = cv,
              bpoint1 = min(bp1,bp2),
              bpoint2 = max(bp1,bp2),
              timeElapse=t1-t0))
}



#==================================


GHansen <- function(y,x,model=1,trim=0.1, use=c("nw","ba")) {
  # Cointegration with structural breaks
  # Gregory and Hansen(1996 Journal of Econometrics)


  #==Bandwidth and kernels
  # nw : kernel = c("ba", "pa", "qs")
  # and: kernel = c("ba", "pa", "qs", "th", "tr")
  y0=y
  y=as.matrix(y)
  x=as.matrix(x)

  k=ncol(x)
  pmax = 12
  t=nrow(y)
  trend=seq(t)/t
  T1 = round(trim*t)
  T2 = round((1-trim)*t)
  idx = T1:T2

  if (model == 1) {

    roll <- function(z) {
      du1 = c(rep(0, z), rep(1, (t - z)))
      rollmat=data.frame(y,x,du1)
      e=as.matrix(resid(lm(rollmat)))
      ADF=urca::ur.df(e,type = "none", lags=pmax,
                      selectlags = c("Fixed", "AIC", "BIC")[2])
      PP=pp(e,type="none", use=use)

      c(ADF@teststat[1,1], PP$Zt, PP$Za)
    }

    roll.stat0 <- sapply(idx, roll)
    roll.stat <-matrix(roll.stat0,length(idx),3,byrow=TRUE)
    colnames(roll.stat)=c("ADF_t","PP_Zt","PP_Za")
    which_adf=which.min(roll.stat[,1])
    which_zt=which.min(roll.stat[,2])
    which_za=which.min(roll.stat[,3])
    adf_t=roll.stat[which_adf,1]
    pp_zt=roll.stat[which_zt,2]
    pp_za=roll.stat[which_za,3]
    bpoint_adf = idx[which_adf]
    bpoint_zt = idx[which_zt]
    bpoint_za = idx[which_za]

    du1.adf = c(rep(0, bpoint_adf), rep(1, (t - bpoint_adf)))
    testmat.adf <- data.frame(y,x,du1.adf)
    test.reg.adf<-lm(testmat.adf)

    du1.za = c(rep(0, bpoint_za), rep(1, (t - bpoint_za)))
    testmat.za <- data.frame(y,x,du1.za)
    test.reg.za<-lm(testmat.za)

    du1.zt = c(rep(0, bpoint_zt), rep(1, (t - bpoint_zt)))
    testmat.zt <- data.frame(y,x,du1.zt)
    test.reg.zt<-lm(testmat.zt)


  } else if (model==2) {


    roll <- function(z) {
      du1 = c(rep(0, z), rep(1, (t - z)))
      rollmat=data.frame(y,x,du1,trend)
      e=as.matrix(resid(lm(rollmat)))
      ADF=urca::ur.df(e,type = "none", lags=pmax,
                      selectlags = c("Fixed", "AIC", "BIC")[2])
      PP=pp(e,type="none",use=use)

      c(ADF@teststat[1,1], PP$Zt, PP$Za)
    }

    roll.stat0 <- sapply(idx, roll)

    roll.stat <-matrix(roll.stat0,length(idx),3,byrow=TRUE)
    colnames(roll.stat)=c("ADF_t","PP_Zt","PP_Za")
    which_adf=which.min(roll.stat[,1])
    which_zt=which.min(roll.stat[,2])
    which_za=which.min(roll.stat[,3])
    adf_t=roll.stat[which_adf,1]
    pp_zt=roll.stat[which_zt,2]
    pp_za=roll.stat[which_za,3]
    bpoint_adf = idx[which_adf]
    bpoint_zt = idx[which_zt]
    bpoint_za = idx[which_za]

    du1.adf = c(rep(0, bpoint_adf), rep(1, (t - bpoint_adf)))
    testmat.adf <- data.frame(y, x, du1.adf, trend)
    test.reg.adf<-lm(testmat.adf)

    du1.za = c(rep(0, bpoint_za), rep(1, (t - bpoint_za)))
    testmat.za <- data.frame(y,x,du1.za,trend)
    test.reg.za<-lm(testmat.za)

    du1.zt = c(rep(0, bpoint_zt), rep(1, (t - bpoint_zt)))
    testmat.zt <- data.frame(y,x,du1.zt,trend)
    test.reg.zt<-lm(testmat.zt)

  } else if (model==3)  {

    roll <- function(z) {
      du1 = c(rep(0, z), rep(1, (t - z)))
      x1=du1*x
      rollmat=data.frame(y,x,du1,x1)
      e=as.matrix(resid(lm(rollmat)))
      ADF=urca::ur.df(e,type = "none", lags=pmax,
                      selectlags = c("Fixed", "AIC", "BIC")[2])
      PP=pp(e,type="none",use=use)

      c(ADF@teststat[1,1], PP$Zt, PP$Za)
    }

    roll.stat0 <- sapply(idx, roll)

    roll.stat <-matrix(roll.stat0,length(idx),3,byrow=TRUE)
    colnames(roll.stat)=c("ADF_t","PP_Zt","PP_Za")
    which_adf=which.min(roll.stat[,1])
    which_zt=which.min(roll.stat[,2])
    which_za=which.min(roll.stat[,3])
    adf_t=roll.stat[which_adf,1]
    pp_zt=roll.stat[which_zt,2]
    pp_za=roll.stat[which_za,3]
    bpoint_adf = idx[which_adf]
    bpoint_zt = idx[which_zt]
    bpoint_za = idx[which_za]

    du1.adf = c(rep(0, bpoint_adf), rep(1, (t - bpoint_adf)))
    x1.adf=du1.adf*x
    testmat.adf <- data.frame(y, x, du1.adf, x1.adf)
    test.reg.adf<-lm(testmat.adf)

    du1.za = c(rep(0, bpoint_za), rep(1, (t - bpoint_za)))
    x1.za=du1.za*x
    testmat.za <- data.frame(y,x,du1.za,x1.za)
    test.reg.za<-lm(testmat.za)

    du1.zt = c(rep(0, bpoint_zt), rep(1, (t - bpoint_zt)))
    x1.zt=du1.zt*x
    testmat.zt <- data.frame(y,x,du1.zt,x1.zt)
    test.reg.zt<-lm(testmat.zt)

  } else if (model==4) {

    roll <- function(z) {
      du1 = c(rep(0, z), rep(1, (t - z)))
      dt1=c(rep(0, z), seq(t - z))/(t-z)
      x1=du1*x
      rollmat=data.frame(y,x,trend,dt1,du1,x1)
      e=as.matrix(resid(lm(rollmat)))
      ADF=urca::ur.df(e,type = "none", lags=pmax,
                      selectlags = c("Fixed", "AIC", "BIC")[2])
      PP=pp(e,type="none",use=use)

      c(ADF@teststat[1,1], PP$Zt, PP$Za)
    }

    roll.stat0 <- sapply(idx, roll)

    roll.stat <-matrix(roll.stat0,length(idx),3,byrow=TRUE)
    colnames(roll.stat)=c("ADF_t","PP_Zt","PP_Za")
    which_adf=which.min(roll.stat[,1])
    which_zt=which.min(roll.stat[,2])
    which_za=which.min(roll.stat[,3])
    adf_t=roll.stat[which_adf,1]
    pp_zt=roll.stat[which_zt,2]
    pp_za=roll.stat[which_za,3]
    bpoint_adf = idx[which_adf]
    bpoint_zt = idx[which_zt]
    bpoint_za = idx[which_za]

    du1.adf = c(rep(0, bpoint_adf), rep(1, (t - bpoint_adf)))
    dt1.adf=c(rep(0, bpoint_adf), seq(t - bpoint_adf))/(t - bpoint_adf)
    x1.adf=du1.adf*x
    testmat.adf <- data.frame(y,x,trend,dt1.adf,du1.adf,x1.adf)
    test.reg.adf<-lm(testmat.adf)

    du1.za = c(rep(0, bpoint_za), rep(1, (t - bpoint_za)))
    dt1.za=c(rep(0, bpoint_za), seq(t - bpoint_za))/(t - bpoint_za)
    x1.za=du1.za*x
    testmat.za <- data.frame(y,x,trend,dt1.za,du1.za,x1.za)
    test.reg.za<-lm(testmat.za)

    du1.zt = c(rep(0, bpoint_zt), rep(1, (t - bpoint_zt)))
    dt1.zt=c(rep(0, bpoint_zt), seq(t - bpoint_zt))/(t - bpoint_zt)
    x1.zt=du1.zt*x
    testmat.zt <- data.frame(y,x,trend,dt1.zt,du1.zt,x1.zt)
    test.reg.zt<-lm(testmat.zt)

  }

  #### Critical values

  if (k==1) {
    if (model==1){
      cvADFZt=cbind(-5.13,-4.61,-4.34)
      cvZa= cbind(-50.07,-40.48,-36.19)
    } else if (model==2) {
      cvADFZt=cbind(-5.45,-4.99,-4.72)
      cvZa= cbind(-57.28,-47.96,-43.22)
    } else if (model==3) {
      cvADFZt=cbind(-5.47,-4.95,-4.68)
      cvZa= cbind(-57.17,-47.04,-41.85)
    } else if (model==4) {
      cvADFZt=cbind(-6.02,-5.50,-5.24)
      cvZa= cbind(-69.37,-58.58,-53.31)
    }

  }  else if (k==2) {
    if (model==1) {
      cvADFZt=cbind(-5.44,-4.92,-4.69)
      cvZa= cbind(-57.01,-46.98,-42.49)
    } else if (model==2) {
      cvADFZt=cbind(-5.80,-5.29,-5.03)
      cvZa= cbind(-64.77,-53.92,-48.94)
    } else if (model==3) {
      cvADFZt=cbind(-5.97,-5.50,-5.23)
      cvZa= cbind(-68.21,-58.33,-52.85)
    } else if (model==4){
      cvADFZt=cbind(-6.45,-5.96,-5.72)
      cvZa= cbind(-79.65,-68.43,-63.10)
    }

  } else if (k==3) {
    if (model==1) {
      cvADFZt=cbind(-5.77,-5.28,-5.02)
      cvZa= cbind(-63.64,-53.58,-48.65)

    } else if (model==2){
      cvADFZt=cbind(-6.05,-5.57,-5.33)
      cvZa= cbind(-70.27,-59.76,-54.94)

    } else if (model==3) {
      cvADFZt=cbind(-6.51,-6.00,-5.75)
      cvZa= cbind(-80.15,-68.94,-63.42)

    } else if (model==4) {
      cvADFZt=cbind(-6.89,-6.32,-6.16)
      cvZa= cbind(-90.84,-78.87,-72.75)
    }
  } else if (k==4) {

    if (model==1) {
      cvADFZt=cbind(-6.05,-5.56,-5.31)
      cvZa= cbind(-70.18,-59.40,-54.38)

    } else if (model==2) {
      cvADFZt=cbind(-6.36,-5.83,-5.59)
      cvZa= cbind(-76.95,-65.44,-60.12)

    } else if (model==3){
      cvADFZt=cbind(-6.92,-6.41,-6.17)
      cvZa= cbind(-90.35,-78.52,-72.56)

    } else if (model==4)
      cvADFZt=cbind(-7.31,-6.84,-6.58)
    cvZa=cbind(-100.69,-88.47,-82.30)

  }

  cvZt=cvADFZt
  if (k < 5) {
    CV=rbind(cvADFZt, cvZt, cvZa)
    rownames(CV)=c("ADF_t","PP_Zt","PP_Za")
    colnames(CV)=c("1%","5%","10%")

    breakpoint=rbind(bpoint_adf,bpoint_zt,bpoint_za)
    rownames(breakpoint)=c("ADF_t","PP_Zt","PP_Za")
    colnames(breakpoint)="break point"
    statistic=rbind(adf_t,pp_zt,pp_za)
    rownames(statistic)=c("ADF_t","PP_Zt","PP_Za")
    colnames(statistic)="statistics"
  } else {message("Gregory-Hansen test critical values unavailable for k>4")}

  rownames(roll.stat)=rownames(y)[idx]

  return(list(result=cbind(CV,statistic,breakpoint),
              teststat=timeSeries::as.timeSeries(roll.stat),
              test.reg.adf=test.reg.adf,
              test.reg.za=test.reg.za,
              test.reg.zt=test.reg.zt))

}


#=== KPSS with 1 break
kpss_1br <-function (y,
                      lags=c("short", "long", "nil"),
                      model=c("intercept","both"),
                      use=c("nw","ba"),
                      trim=0.1){
  #==Bandwidth
  # nw : kernel = c("ba", "pa", "qs")
  # and: kernel = c("ba", "pa", "qs", "th", "tr")
  lags <- match.arg(lags)
  model <- match.arg(model)
  y0=y
  y = as.matrix(y)
  if(ncol(y)>1) {stop("\nError: Only one series is allowed.")}

  t=nrow(y)
  trend=seq(t)/t
  T1 = round(trim * t);
  T2 = round((1 - trim) * t);
  idx = T1:T2

  if (model == "intercept") {

    roll <- function(z) {
      du1 = c(rep(0, z), rep(1, (t - z)))
      x=cbind(du1)
      kpss(y,x,lags=lags,use=use)$teststat
    }

    roll.stat <- sapply(idx, roll)
    bpoint = idx[which.min(roll.stat)]
    teststat=min(roll.stat)
    cval = matrix(c(0.28299, 0.37538, 0.60388,
                    0.22915,0.30212, 0.48265,
                    0.18678,0.24247, 0.38052,
                    0.16007,0.20106, 0.30162,
                    0.15176,0.18688, 0.26842),5,3,byrow=TRUE)


  } else if (model=="both") {

    roll <- function(z) {
      du1 = c(rep(0, z), rep(1, (t - z)))
      dt1 = c(rep(0, z), 1:(t - z))/t
      x=cbind(trend,dt1,du1)
      kpss(y,x,lags,use=use)$teststat
    }

    roll.stat <- sapply(idx, roll)
    bpoint = idx[which.min(roll.stat)]
    teststat=min(roll.stat)

    cval = matrix(c(0.09724, 0.12046, 0.17704,
                    0.09724, 0.12046, 0.14208,
                    0.06485, 0.07889, 0.11308,
                    0.05570, 0.06615, 0.09122,
                    0.05267, 0.06163, 0.08216),5,3,byrow=TRUE)


  }

  TBi = round(10 * (bpoint/t))
  if (TBi > 5) { TBi = 10 - TBi}
  cv=t(as.matrix(rev(cval[TBi,,drop=F])))

  rownames(cv)="critical values"
  colnames(cv)=c("1%","5%","10%")

  return(list(teststat=teststat,
              cval=cv,
              bpoint=bpoint,
              tstats=roll.stat))
}




#=== KPSS with two breaks
kpss_2br <-function (y,
                      lags=c("short", "long", "nil"),
                      model=1,
                      use=c("nw","ba"),
                      trace=TRUE){
  #==Bandwidth
  # nw : kernel = c("ba", "pa", "qs")
  # and: kernel = c("ba", "pa", "qs", "th", "tr")
  lags <- match.arg(lags)
  y=as.matrix(y)

  if(ncol(y)>1) {stop("\nError: Only one series is allowed.")}

  t = nrow(y)
  trend=seq(t)/(t)

  tb1 = 2
  tb2 = tb1 + 2
  idx1=tb1:(t-4)
  idx2=tb2:(t-2)

  if (model == 1) {#cbind(du1,du2)


    roll1 <- function(z) { #z=680

      if(isTRUE(trace)) {(t-z-3)}

      du1 = c(rep(0,z),rep(1,t-z)) #tb1

      roll2 <- function(z) {#z=70

        du2 = c(rep(0, z), rep(1, (t - z)))
        x=cbind(du1,du2)
        kpss(y,x,lags=lags, use=use)$teststat
      }

      sapply(idx2, roll2)

    }


    t0=proc.time();rollStat <- sapply(idx1, roll1); t1=proc.time()


  } else if (model == 2) {#cbind(trend,du1,du2)


    roll1 <- function(z) { #tb1 loop

      if(isTRUE(trace)) {t-z-3}

      du1 = c(rep(0,z),rep(1,t-z))


      roll2 <- function(z) { #tb2 loop
        du2 = c(rep(0,z),rep(1,t-z))

        x=cbind(trend,du1,du2)
        kpss(y,x,lags=lags, use=use)$teststat
      }

      sapply(idx2, roll2)
    }

    t0=proc.time(); rollStat <- sapply(idx1, roll1); t1=proc.time()

  } else if (model == 3) {#cbind(trend,dt1,dt2)


    roll1 <- function(z) { #tb1 loop

      if(isTRUE(trace)) {t-z-3}

      dt1 = c(rep(0,z),seq(t-z))/t

      roll2 <- function(z) { #tb2 loop

        dt2 = c(rep(0,z),seq(t-z))/t
        x=cbind(trend,dt1,dt2)
        kpss(y,x,lags=lags, use=use)$teststat
      }

      sapply(idx2, roll2)
    }

    t0=proc.time(); rollStat <- sapply(idx1, roll1); t1=proc.time()

  } else if (model == 4) {#cbind(trend,du1,du2,dt1,dt2)


    roll1 <- function(z) { #tb1 loop

      if(isTRUE(trace)) {(t-z-3)}

      du1 = c(rep(0,z),rep(1,t-z))
      dt1 = c(rep(0,z),seq(t-z))/t

      roll2 <- function(z) { #tb2 loop
        du2 = c(rep(0,z),rep(1,t-z))
        dt2 = c(rep(0,z),seq(t-z))/t
        x=cbind(trend,du1,du2,dt1,dt2)
        kpss(y,x,lags=lags, use=use)$teststat
      }

      sapply(idx2, roll2)
    }

    t0=proc.time(); rollStat <- sapply(idx1, roll1); t1=proc.time()

  }




  # Critical Values
  if (model == 1) {
    mat_cv1 = matrix(c(0.4758,0.3659,0.2802,0.2299,0.2275,0.2883,0.3664,0.4758,
                       0,0.3682,0.2832,0.2109,0.1835,0.2075,0.2897,0.3612,
                       0,0,0.2874,0.2077,0.1588,0.1620,0.2178,0.2937,
                       0,0,0,0.2330,0.1733,0.1648,0.1811,0.2292,
                       0,0,0,0,0.2271,0.2109,0.2027,0.2239,
                       0,0,0,0,0,0.2919,0.2846,0.2862,
                       0,0,0,0,0,0,0.3666,0.3692,
                       0,0,0,0,0,0,0,0.4853),8,8,byrow=TRUE)

    mat_cv5 = matrix(c(0.2992,0.2344,0.1890,0.1560,0.1581,0.1883,0.2339,0.3001,
                       0,0.2339,0.1802,0.1423,0.1289,0.1390,0.1821,0.2366,
                       0,0,0.1846,0.1401,0.1153,0.1165,0.1421,0.1859,
                       0,0,0,0.1585,0.1266,0.1165,0.1275,0.1571,
                       0,0,0,0,0.1564,0.1443,0.1388,0.1547,
                       0,0,0,0,0,0.1834,0.1830,0.1847,
                       0,0,0,0,0,0,0.2328,0.2332,
                       0,0,0,0,0,0,0,0.3009),8,8,byrow=TRUE)

    mat_cv10 = matrix(c(0.2262,	0.1789,	0.1441,	0.1258,	0.1286,	0.1455,	0.177,	0.226,
                        0,0.1775,	0.1402,	0.1148,	0.1053,	0.1116,	0.1394,	0.1791,
                        0,0,0.1432,0.1128,	0.0959,	0.0963,	0.1136,	0.1441,
                        0,0,0,0.1276,0.1049,0.0965,0.1043,0.1279,
                        0,0,0,0,0.1564,0.1151,0.1120,	0.1264,
                        0,0,0,0,0,0.1440,0.1409,0.1438,
                        0,0,0,0,0,0,0.1774,0.1785,
                        0,0,0,0,0,0,0,0.2289),8,8,byrow=TRUE)


  } else if (model == 2) {
    mat_cv1 = matrix(c(0.1456,0.1169,0.1247,0.1601,0.1616,0.127,0.1157,0.1444,
                       0,0.1192,0.1002,0.1219,0.14,0.1175,0.1028,0.1214,
                       0,0,0.1252,0.1187,0.1253,0.1248,0.1147,0.1268,
                       0,0,0,0.1604,0.1405,0.1272,0.1399,0.1574,
                       0,0,0,0,0.1679,0.1191,0.1159,0.1590,
                       0,0,0,0,0,0.1223,0.1020,0.1260,
                       0,0,0,0,0,0,0.1205,0.1180,
                       0,0,0,0,0,0,0,0.1421),8,8,byrow=TRUE)

    mat_cv5 = matrix(c(0.0988,0.0834,0.0899,0.1073,0.108,0.091,0.0825,0.0992,
                       0,0.0843,0.075,0.086,0.0948,0.0849,0.0749,0.0845,
                       0,0,0.0895,0.0848,0.0898,0.0886,	0.0831,0.0898,
                       0,0,0,0.1069,0.0958,0.0896,0.0934,0.1051,
                       0,0,0,0,0.1087,0.0833,0.0836,0.1061,
                       0,0,0,0,0,0.0894,0.0746,0.0900,
                       0,0,0,0,0,0,0.0856,0.0843,
                       0,0,0,0,0,0,0,0.0999),8,8,byrow=TRUE)

    mat_cv10 = matrix(c(0.0797,0.0699,0.0748,0.0858,0.0855,0.075,0.0693,0.0801,
                        0,0.0701,0.0641,0.071,0.0771,0.0702,0.0643,0.0696,
                        0,0,0.0739,0.0706,0.074,0.0736,0.0691,0.0741,
                        0,0,0,0.0860,0.0771,0.0742,0.0747,0.085,
                        0,0,0,0,0.0859,0.0699,0.0698,0.0867,
                        0,0,0,0,0,0.0741,0.0633,0.0743,
                        0,0,0,0,0,0,0.0714,0.0699,
                        0,0,0,0,0,0,0,0.0815),8,8,byrow=TRUE)

  } else if (model == 3) {
    mat_cv1 = matrix(c(0.1524,0.1298,0.1098,0.1040,0.1003,0.1175,0.1393,0.1565,
                       0,0.1205,0.1021,0.0892,0.0879,0.0963,0.1134,0.1357,
                       0,0,0.0966,0.0820,0.0780,0.0827,0.0974,0.1153,
                       0,0,0,0.0797,0.0751,0.0766,0.0888,0.1040,
                       0,0,0,0,0.0803,0.0829,0.0910,0.1032,
                       0,0,0,0,0,0.0968,0.0993,0.1089,
                       0,0,0,0,0,0,0.1200,0.1282,
                       0,0,0,0,0,0,0,0.1518),8,8,byrow=TRUE)

    mat_cv5 = matrix(c(0.1060,0.0897,0.0766,0.0723,0.0738,0.0803,0.0937,0.1068,
                       0,0.0822,0.0715,0.0650,0.0634,0.0681,0.0785,0.0929,
                       0,0,0.0679,0.0601,0.0573,0.0599,0.0683,0.0802,
                       0,0,0,0.0591,0.0556,0.0565,0.0633,0.0737,
                       0,0,0,0,0.0590,0.0590,0.0646,0.0731,
                       0,0,0,0,0,0.0684,0.0709,0.0772,
                       0,0,0,0,0,0,0.0818,0.0886,
                       0,0,0,0,0,0,0,0.1021),8,8,byrow=TRUE)

    mat_cv10 = matrix(c(0.0848,0.0729,0.0630,0.0603,0.0613,0.0664,0.0755,0.0864,
                        0,0.0669,0.0593,0.0540,0.0529,0.0562,0.0644,0.0754,
                        0,0,0.0561,0.0504,0.0484,0.0502,0.0568,0.0659,
                        0,0,0,0.0500,0.0473,0.0483,0.0528,0.0609,
                        0,0,0,0,0.0502,0.0500,0.0539,0.0606,
                        0,0,0,0,0,0.0570,0.0584,0.0634,
                        0,0,0,0,0,0,0.0674,0.0722,
                        0,0,0,0,0,0,0,0.0831),8,8,byrow=TRUE)

  } else if (model == 4) {

    mat_cv1 = matrix(c(0.1439,0.1116,0.0852,0.0691,0.0704,0.0856,0.1098,0.1430,
                       0,0.1100,0.0835,0.0644,0.0556,0.0634,0.0849,0.1113,
                       0,0,0.0855,0.0652,0.0506,0.0503,0.0637,0.0845,
                       0,0,0,0.0699,0.0560,0.0501,0.0550,0.0695,
                       0,0,0,0,0.0704,0.0629,0.0637,0.0707,
                       0,0,0,0,0,0.0874,0.0840,0.0858,
                       0,0,0,0,0,0,0.1122,0.1124,
                       0,0,0,0,0,0,0,0.1425),8,8,byrow=TRUE)

    mat_cv5 = matrix(c(0.0972,0.0772,0.0605,0.0518,0.052,0.0606,0.0765,0.0965,
                       0,0.0763,0.0591,0.047,0.0424,0.0466,0.0586,0.0757,
                       0,0,0.0601,0.0474,0.039,0.0389,0.0467,0.0605,
                       0,0,0,0.0518,0.0425,0.0389,0.0423,0.0518,
                       0,0,0,0,0.0521,0.0466,0.047,0.0523,
                       0,0,0,0,0,0.0615,0.0586,0.0608,
                       0,0,0,0,0,0,0.0768,0.0766,
                       0,0,0,0,0,0,0,0.0966),8,8,byrow=TRUE)

    mat_cv10 = matrix(c(0.0788,0.0626,0.0508,0.0441,0.0444,0.0504,0.0625,0.0775,
                        0,0.0621,0.0487,0.0399,0.036,0.0397,0.0485,0.0619,
                        0,0,0.0501,0.0398,0.0339,0.0340,0.0397,0.0500,
                        0,0,0,0.0442,0.0367,0.0340,0.0365,0.0442,
                        0,0,0,0,0.0444,0.0397,0.0399,0.0443,
                        0,0,0,0,0,0.0507,0.0486,0.0503,
                        0,0,0,0,0,0,0.0620,0.0625,
                        0,0,0,0,0,0,0,0.0780),8,8,byrow=TRUE)


  }

  lam1 = tb1/t
  lam2 = tb2/t

  if  (lam1 <= 0.15) {lam1row = 1
  } else if (0.15< lam1 & lam1 <= 0.85) {
    lam1row = round(10*lam1)
  } else if (lam1 > 0.85) {lam1row = 8}


  if  (lam2 <= 0.25) {lam2col = 1
  } else if (0.25< lam2 & lam2 <= 0.95) {lam2col = round(10*lam2)
  } else if (lam2 > 0.95) {lam2col = 8}

  cv1 = mat_cv1[lam1row,lam2col]
  cv5 = mat_cv5[lam1row,lam2col]
  cv10= mat_cv10[lam1row,lam2col]



  cv=cbind(cv1,cv5,cv10)
  colnames(cv)=c("1%","5%","10%")
  rownames(cv) <- "critical values"

  R=which.min(apply(rollStat,1,min)) #row number, tb2
  C=which.min(rollStat[R,]) #column number, tb1
  bp1=idx1[C]
  bp2=idx2[R]

  return(list(teststat = min(rollStat),
              cval = cv,
              bpoint1 = min(bp1,bp2),
              bpoint2 = max(bp1,bp2),
              timeElapse=t1-t0))
}





#==Kernels used for PP and KPSS tests
Bartlett_uni <-function (e, v) {

  e=as.matrix(e)

  t = nrow(e)

  lrv = t(e)%*%e/t

  for (i in 1:v) {
    w   = 1- i/(v+1)
    lrv = lrv + 2 * t(e[seq(t-i)]) %*% e[(1+i):t] %*% w/t
  }

  return(as.numeric(lrv))

}



QS_uni <- function (e, v) {

  e=as.matrix(e)

  t = nrow(e)
  lrv = t(e)%*%e/t

  for (i in seq(v)) {
    x1 = i/v
    x2 = 6*pi*x1/5
    w  = (25/(12*(pi*x1)^2)) * (sin(x2)/x2 - cos(x2));
    lrv= lrv + 2 * t(e[1:(t-i)]) %*% e[(1+i):t] * w/t
  }

  return(as.numeric(lrv))
}


Parzen_uni <-function(e,v){
  e=as.matrix(e)

  if (v > nrow(e)) {v = nrow(e)-1}

  if (v >= 1) { weights = rep(0,v) }

  A = 0
  for ( i  in seq(v)) {#i=1
    if (i<=max(seq(v/2))) {
      m  =  1  -  6*((i/(v+1))^2) + 6*((abs(i)/(v+1))^3)
      t1  =  apply(e[-seq(i),,drop=F],2, function(x) x-mean(x))
      t2  =  apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x))
      A  =  A  +  m * (t(t1)%*%t2)
      weights[i] = m}

    else {
      m  =  2*(1 - (abs(i)/(v+1)))^3
      t1  =  apply(e[-seq(i),,drop=F],2, function(x) x-mean(x))
      t2  =  apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x))
      A  =  A  +  m*(t(t1)%*%t2)
      weights[i] = m

    }
  }
  colnames(A)=rownames(A)=colnames(e)
  return(as.numeric(A/nrow(e)))
}

SPC_Bartlett <- function (e, v) {


  e=as.matrix(e)
  t = nrow(e)
  if (v <= 0) {v = 1}
  min_BIC = 1e25

  rho=NULL;for (i in 1:v) { #i=1


    temp = embed(e, i+1)

    # Autoregressive coefficients
    REG1=lm(V1~.-1,data=as.data.frame(temp))
    rho_temp = coef(REG1)

    # OLS residuals
    res_temp = as.matrix(resid(REG1))

    BIC = log(t(res_temp) %*% res_temp/(t-v)) + (i * log(t - v)/(t - v))

    if (BIC < min_BIC) {
      min_BIC = BIC
      k = i
      rho = as.matrix(rho_temp)

      #Prewithening
      res = res_temp
    }

  }

  res=as.matrix(res)
  temp = embed(res,2)
  tail(temp)
  # We use an AR(1) approximation as in Andrews and Monahan (1992, p. 958)
  REG2=lm(V1~V2-1,data=as.data.frame(temp))
  a = coef(REG2)

  # Obtaining the bandwidth for the spectral window
  l = as.numeric(1.1447 * (4 * a^2*t/((1 + a)^2 * (1 - a)^2))^(1/3))

  # Truncate the estimated bandwidth
  l = trunc(l)
  res=as.matrix(res)
  #Short-run variance
  lrv = t(res)%*% res/t

  if (l==0) {
    i=1
    w = (1 - i/(l + 1))
    lrv = lrv + 2 * t(res[1:(nrow(res)-i)])%*% res[(1+i):nrow(res)] * w/t

  } else {
    for (i in 1:l) {

      # Bartlett kernel
      w = (1 - i/(l + 1))
      lrv = lrv + 2 * t(res[1:(nrow(res)-i)])%*% res[(1+i):nrow(res)] * w/t
    }
  }


  #Recoloring
  lrv_recolored = lrv/(1 - sum(rho))^2

  # Sul, Phillips and Choi (2003) boundary rule
  lrv = min(lrv_recolored,t * lrv)
  return(as.numeric(lrv))
}



SPC_QS <- function (e, v) {


  e=as.matrix(e)
  t = nrow(e)
  if (v <= 0) {v = 1}
  min_BIC = 1e25

  rho=NULL;for (i in 1:v) { #i=1


    temp = embed(e, i+1)

    # Autoregressive coefficients
    REG1=lm(V1~.-1,data=as.data.frame(temp))
    rho_temp = coef(REG1)

    # OLS residuals
    res_temp = as.matrix(resid(REG1))

    BIC = log(t(res_temp) %*% res_temp/(t-v)) + (i * log(t - v)/(t - v))

    if (BIC < min_BIC) {
      min_BIC = BIC
      k = i
      rho = as.matrix(rho_temp)

      #Prewithening
      res = res_temp
    }

  }

  res=as.matrix(res)
  temp = embed(res,2)
  tail(temp)
  # We use an AR(1) approximation as in Andrews and Monahan (1992, p. 958)
  REG2=lm(V1~V2-1,data=as.data.frame(temp))
  a = coef(REG2)

  # Obtaining the bandwidth for the spectral window
  l = as.numeric(1.3221 * (4 * a^2 * t/((1 + a)^2 * (1 - a)^2))^(1/5))

  # Truncate the estimated bandwidth
  l = trunc(l)
  res=as.matrix(res)
  #Short-run variance
  lrv = t(res)%*% res/t

  if (l==0) {
    i=1
    w = 25/(12 * pi^2 * (i/l)^2) * (sin(6 * pi * i/(l * 5))/(6 * pi * i/(l * 5)) - cos(6 * pi * i/(l * 5)))
    lrv = lrv + 2 * t(res[1:(nrow(res)-i)])%*% res[(1+i):nrow(res)] * w/t

  } else {
    for (i in 1:l) {

      # Bartlett kernel
      w = 25/(12 * pi^2 * (i/l)^2) * (sin(6 * pi * i/(l * 5))/(6 * pi * i/(l * 5)) - cos(6 * pi * i/(l * 5)))
      lrv = lrv + 2 * t(res[1:(nrow(res)-i)])%*% res[(1+i):nrow(res)] * w/t
    }
  }


  #Recoloring
  lrv_recolored = lrv/(1 - sum(rho))^2

  # Sul, Phillips and Choi (2003) boundary rule
  lrv = min(lrv_recolored,t * lrv)
  return(as.numeric(lrv))
}


Kurozumi_Bartlett<- function(e) {

  e=as.matrix(e)
  t = nrow(e)

  # AR(1) estimate
  temp = embed(e, 2)
  a = coef(lm(V1~V2-1,data=as.data.frame(temp)))

  # Defines the upper bound
  k=0.7

  # Bandwidth
  l = min(1.1447 * (4 * a^2 * t/((1 + a)^2 * (1 - a)^2))^(1/3)
          , 1.1447 * (4 * k^2 * t/((1 + k)^2 * (1 - k)^2))^(1/3))

  # Truncate the estimated bandwidth
  l = trunc(l);

  # Conventional variance estimation
  lrv = t(e)%*%e/t;

  if (l==0) {
    i=1
    w = (1 - i/(l + 1))
    lrv = lrv + 2 * t(e[1:(nrow(e)-i)])%*% e[(1+i):nrow(e)] * w/t

  } else {
    for (i in 1:l) {

      # Bartlett kernel
      w = (1 - i/(l + 1))
      lrv = lrv + 2 * t(e[1:(nrow(e)-i)])%*% e[(1+i):nrow(e)] * w/t
    }
  }
  as.numeric(lrv)
}



Kurozumi_QS <-function(e) {

  e=as.matrix(e)
  t = nrow(e)

  # AR(1) estimate
  temp = embed(e, 2)
  a = coef(lm(V1~V2-1,data=as.data.frame(temp)))

  # Defines the upper bound
  k=0.7

  # Bandwidth
  l = min(1.3221 * (4 * a^2*t/((1 + a)^2 * (1 - a)^2))^(1/5)
          ,1.3221 * (4 * k^2*t/((1 + k)^2 * (1 - k)^2))^(1/5))

  # Truncate the estimated bandwidth
  l = trunc(l)

  # Short-run variance
  lrv = t(e)%*%e/t;
  if (l==0) {
    i=1
    w = 25/(12 * pi^2 * (i/l)^2) * (sin(6 * pi * i/(l * 5))/(6 * pi * i/(l * 5)) - cos(6 * pi * i/(l * 5)))
    lrv = lrv + 2 * t(e[1:(nrow(e)-i)])%*% e[(1+i):nrow(e)] * w/t

  } else {
    for (i in 1:l) {

      # Bartlett kernel
      w = 25/(12 * pi^2 * (i/l)^2) * (sin(6 * pi * i/(l * 5))/(6 * pi * i/(l * 5)) - cos(6 * pi * i/(l * 5)))
      lrv = lrv + 2 * t(e[1:(nrow(e)-i)])%*% e[(1+i):nrow(e)] * w/t
    }
  }
  as.numeric(lrv)
}


































.ZA_2br_Devlp <-function(y,
                            ic=c("AIC","BIC")[1],
                            model = c("intercept","both")[1],
                            pmax = 8,
                            trim=0.1,
                            eq=1,
                        trace=TRUE) {
  y=as.matrix(y)
  if(ncol(y)>1) {stop("\nError: Only one series is allowed.")}

  t = nrow(y)
  trend=seq(t)/(t)

  #T1 = round(trim * t)
  #T2 = round((1 - trim) * t)
  T1 = rbind(max((3 + pmax), ceiling(trim * t)))
  T2 = min(rbind(t - 3 -pmax),floor((1 - trim) * t))

  tb1 = T1

  #if (model == "intercept") {



  if (T1 < pmax+2) {T1 = pmax + 3}


  tb2 = tb1 + 2

  idx1=tb1:T2
  idx2=tb2:T2


  t0=proc.time();rollStat=matrix(rep(0,length(idx1)*length(idx2)),length(idx2),length(idx1));
  for (tb1 in idx1) { #tb1=68

    if(isTRUE(trace)) {(T2-tb1+1)}

    du1 = c(rep(0,tb1),rep(1,t-tb1)) #tb1

    out=.plag(y,x=NULL,pmax,ic,eq)
    datmat = out$data
    p=out$p

    roll2 <- function(z) {#tb2
      du2 = c(rep(0,z),rep(1,t-z))
      du1=embed(du1,p+2)[,1,drop=F]
      du2=embed(du2,p+2)[,1,drop=F]
      rollmat <- data.frame(datmat,du1, du2)

      COEF <- coef(summary(lm(rollmat)))
      if (eq==1) {
        (COEF[2, 1] - 1)/COEF[2, 2]
      } else {
        COEF[2, 3]
      }
    }
    rollStat[,T1-67] <- sapply(idx2, roll2)


  };t1=proc.time()



  # Critical Values (see, Narayan & Popp, 2010, Table 3)
  if (model == "intercept") {
    if (t <= 50) {
      cv = cbind(-5.259, -4.514,-4.143)
    } else if (50 < t & t <= 200) {

      cv = cbind(-4.958,-4.316,-3.980)

    } else if (200 < t & t <= 400) {

      cv = cbind(-4.731,-4.136,-3.825)

    } else if (400 < t) {

      cv = cbind(-4.672,-4.081,-3.772)

    }
  } else if (model == "both") {
    if (t <= 50) {
      cv =  cbind(-5.949,-5.181,-4.789)
    } else if (50 < t & t <= 20){
      cv =  cbind(-5.576,-4.937,-4.596)
    }   else if (200 < t & t <= 400) {
      cv =  cbind(-5.318,-4.741,-4.430)
    } else if (400 < t) {
      cv =  cbind(-5.287,-4.692,-4.396)
    }
  }

  colnames(cv)=c("1%","5%","10%")
  rownames(cv) <- "critical values"
  R=which.min(apply(rollStat,1,min)) #row number, tb2
  C=which.min(rollStat[R,]) #column number, tb1
  bp1=idx1[C]
  bp2=idx2[R]


  return(list(y = y,
              teststat = min(rollStat),
              cval = cv,
              bpoint1 = min(bp1,bp2),
              bpoint2 = max(bp1,bp2),
              timeElapse=t1-t0,
              test.name = "Narayan-Popp"))


}










.ADF.1br <-function(y,ic=c("AIC","BIC"),
                   model=c("intercept", "trend","both"),
                   outlier=1,
                   pmax=8,
                   trim=0.1) {
  y=as.matrix(y)

  if(ncol(y)>1) {stop("\nError: Only one series is allowed.")}

  ic <- match.arg(ic)
  model <- match.arg(model)

  t = nrow(y)
  tb1_min = 0
  ADF_min = 1000

  T1 = round(trim * t)
  T2 = round((1 - trim) * t)

  if (T1 < pmax+2) {T1 = pmax + 3}

  tb1 = T1

  #== Searching for breaks
  t0=proc.time(); Stat=NULL;for (tb1 in T1:T2) { #tb1=T1

    dt = seq(t)/t # trend dummy

    du1 = c(rep(0,tb1),rep(1,t-tb1))
    dt1 =  c(rep(0,tb1),seq(t-tb1))/t
    br1 = c(rep(0,tb1),1, rep(0,t-tb1-1))

    if (outlier == 1)   {# Innovational outlier model
      taup = aicp = sicp = tstatp = NULL
      for (p in c(0,seq(pmax))) { #p=3

        if (model == "intercept") { #dummy matrix
          z = cbind(du1,br1)
        } else if (model == "trend") {
          z = cbind(dt,dt1,br1)
        } else if (model == "both") {
          z = cbind(du1,dt,dt1,br1)
        }

        datmat <- .datmat4ADF(y,x=z,p)

        REG=lm(datmat)
        e1=as.matrix(resid(REG))
        COEF=coef(summary(REG))
        taup = c(taup,(COEF[2,1]-1)/COEF[2,2])
        aicp = c(aicp,log((t(e1)%*%e1)/nrow(e1)) + 2 * (ncol(datmat) -1 + 2)/nrow(e1))
        sicp=c(sicp,log((t(e1)%*%e1)/nrow(e1)) +  (ncol(datmat) -1 + 2) * log(nrow(e1))/nrow(e1))
        tstatp = c(tstatp,abs(COEF[nrow(COEF),3]))
      }

      if (ic=="AIC") {p = which.min(aicp)} else if (ic=="BIC") {p = which.min(sicp)}

      stat = taup[p]

      if (stat < ADF_min) {
        tb1_min = tb1 # First break date


        # T-statistic with breaks
        ADF_min = stat

        # Optimal lag
        opt_lag = p-1
      }

    Stat=c(Stat,stat)

    } else if (outlier == 2) { #Additive outlier model

      if (model == "intercept") {

        x = du1

      } else if (model=="trend") {

        x = cbind(dt, dt1)

      } else if (model=="both") {

        x = cbind(du1, dt, dt1)

      }

      resid =as.matrix(resid(lm(y~x)))
      taup = aicp = sicp = tstatp = NULL
      for (p in c(0,seq(pmax))) {

        datmat <- .datmat4ADF(resid,x=br1,p)

        REG=lm(datmat)
        e1=as.matrix(resid(REG))
        COEF=coef(summary(REG))
        taup = c(taup,(COEF[2,1]-1)/COEF[2,2])
        aicp = c(aicp,log((t(e1)%*%e1)/nrow(e1)) + 2*(ncol(datmat) -1 + 2)/nrow(e1))
        sicp=c(sicp,log((t(e1)%*%e1)/nrow(e1)) +  (ncol(datmat) -1 + 2)*log(nrow(e1))/nrow(e1))
        tstatp = c(tstatp,abs(COEF[nrow(COEF),3]))

}

      if (ic=="AIC") {p  =which.min(aicp)} else if (ic=="BIC") {p  =which.min(sicp)}
      stat = taup[p]

      if (stat < ADF_min) {
        tb1_min = tb1 # First break date


        # T-statistic with breaks
        ADF_min = stat

        # Optimal lag
        opt_lag = p-1
      }
      Stat=c(Stat,stat)
    }   # end of outlier==2


    if (model=="intercept"){cval <- cbind(-5.34, -4.80, -4.58)
    } else if (model=="trend") {    cval <-cbind(-4.93, -4.42, -4.11)
    } else {cval <- cbind(-5.57, -5.08, -4.82)}

  }; t1=proc.time()  # end of tb1 loop


  timeElapse=t1-t0
  colnames(cval)=c("1%","5%","10%")
  rownames(cval) <- "critical values"
  return(list(teststat=ADF_min,
              bpoint=tb1_min,
              p=opt_lag,
              cval=cval,
              tstats=Stat,
              timeElapse=timeElapse,
              test.name = "Zivot-Andrews"))
}



.ADF.2br <-function(y,
                   ic=c("AIC","BIC"),
                   model=c("intercept", "both"),
                   pmax=8,
                   trim=0.1) {
  y=as.matrix(y)
  if(ncol(y)>1) {stop("\nError: Only one series is allowed.")}

  ic <- match.arg(ic)
  model <- match.arg(model)

  t = nrow(y)
  tb1_min = 0
  tb2_min = 0
  ADF_min = 1000

  T1 = round(trim * t)
  T2 = round((1 - trim) * t)
  t1 = max((3 + pmax), ceiling(trim * t))
  t2 = min((t - 3 -pmax),floor((1 - trim) * t))

  if (T1 < pmax+2) {T1 = pmax + 3}

  tb1 = T1

  #== Loop: Searching for breaks

  t0=proc.time();Stat=NULL;  for (tb1 in T1:T2) { #tb1=T1

    # Bounds as in LS
    if (model == "intercept") {tb2 = tb1 + 2
    } else  {tb2 = tb1 + 3}

    id=tb2

    for (tb2 in id:T2) { #tb2=id

      dt = seq(t)/t # Deterministic trend

      if (model == "intercept"){
        du1 = c(rep(0,tb1),rep(1,t-tb1))
        du2 = c(rep(0,tb2),rep(1,t-tb2))
        z = cbind(du1,du2)

      } else if (model == "trend") {

        dt1 =  c(rep(0,tb1),seq(t-tb1))/t
        dt2 =  c(rep(0,tb2),seq(t-tb2))/t
        z = cbind(dt,du2,dt1,dt2)

      } else if (model == "both") {
        du1 =  c(rep(0,tb1),rep(1,t-tb1))
        du2 =  c(rep(0,tb2),rep(1,t-tb2))
        dt1 =  c(rep(0,tb1),seq(t-tb1))/t
        dt2 =  c(rep(0,tb2),seq(t-tb2))/t
        z = cbind(dt,du1,du2,dt1,dt2)
      }

      t = nrow(y)

      Y=embed(y,2)
      dy = as.matrix(Y[,1]-Y[,2])
      y1 = Y[,2,drop=F]
      z1 = embed(z,2)[,-seq(ncol(z)),drop=F]

      dep = dy
      ly = y1
      temp=cbind(dy,ly)
      taup = aicp = sicp = tstatp = NULL
      for (p in c(0,seq(pmax))) { #p=3
        dep=embed(temp,p+1)[,1,drop=F]
        ly=embed(temp,p+1)[,2,drop=F]
        ldy = embed(dy,p+1)[,-1,drop=F]
        if(p==0) {lz = z1} else {
          lz = as.matrix(z1)[-seq(p),,drop=F] }
        rownames(dep)=rownames(ldy)=rownames(ly)=rownames(lz) =NULL


        if (p == 0) {
          x = cbind(ly,lz)
        } else if (p > 0) {
          x = cbind(ly,lz,ldy[, 1:p])
        }
        REG=lm(dep~x)
        e1=as.matrix(resid(REG))
        COEF=coef(summary(lm(REG)))
        taup = c(taup,COEF[2,1]/COEF[2,2])
        aicp = c(aicp,log((t(e1)%*%e1)/nrow(e1)) + 2 * (ncol(x) -1 + 2)/nrow(e1))
        sicp=c(sicp,log((t(e1)%*%e1)/nrow(e1)) +  (ncol(x) -1 + 2) * log(nrow(e1))/nrow(e1))
        tstatp = c(tstatp,abs(COEF[nrow(COEF),3]))
      }

      if (ic=="AIC") {p  =which.min(aicp)} else if (ic=="BIC") {p  =which.min(sicp)}

      stat = taup[p]
      Stat=c(Stat,stat)
      if (stat < ADF_min) {
        tb1_min = tb1 # First break date


        tb2_min=tb2 # Second break date

        # T-statistic with breaks
        ADF_min = stat

        # Optimal lag
        opt_lag = p-1
      }

      message(paste0("tb1=",tb1,"/",T2,"; tb2=",tb2))
    } # end of tb2 loop


  };t1=proc.time()  # end of tb1 loop


  timeElapse=t1-t0

  # Critical Values (see, Narayan & Popp, 2010, Table 3)
  if (model == "intercept") {
    if (t <= 50) {
      cv = cbind(-5.259, -4.514,-4.143)
    } else if (50 < t & t <= 200) {

      cv = cbind(-4.958,-4.316,-3.980)

    } else if (200 < t & t <= 400) {

      cv = cbind(-4.731,-4.136,-3.825)

    } else if (400 < t) {

      cv = cbind(-4.672,-4.081,-3.772)

    }
  } else {
    if (t <= 50) {
      cv =  cbind(-5.949,-5.181,-4.789)
    } else if (50 < t & t <= 20){
      cv =  cbind(-5.576,-4.937,-4.596)
    }   else if (200 < t & t <= 400) {
      cv =  cbind(-5.318,-4.741,-4.430)
    } else if (400 < t) {
      cv =  cbind(-5.287,-4.692,-4.396)
    }
  }
  colnames(cv)=c("1%","5%","10%")
  rownames(cv) <- "critical values"
  return(list(teststat=ADF_min,
              tb1_min=tb1_min,
              tb2_min=tb2_min,
              p=opt_lag,
              cval=cv,
              tstats=Stat,
              timeElapse=timeElapse,
              test.name = "Narayan-Popp"))
}








.datmat4ADF <- function (y,x=NULL,p=8,eq=1) {

  Y=na.omit(as.matrix(y))
  if(ncol(Y)>1) {stop("\nError: Only one series is allowed.")}
  temp=embed(Y,2)


  if(eq==1) { #AR(1) y(t)=y(t-1)+....

    if (p==0)  {
      y=temp[,1,drop=F]
      y.L1=temp[,2,drop=F]
      datmat=cbind(y=y, y.L1=y.L1)
      colnames(datmat)=c("y","y.L1")
    } else {

      y=temp[,1,drop=F]
      y.L1=temp[,2,drop=F]
      Dy=y-y.L1
      dy=embed(Dy,p+1)[,-1,drop=F]
      colnames(dy)=paste("dy.L", 1:p, sep = "")

      y=embed(y,p+1)[,1,drop=F]
      y.L1=embed(y.L1,p+1)[,1,drop=F]
      datmat=cbind(y=y, y.L1=y.L1)
      colnames(datmat)=c("y","y.L1")
      datmat=cbind(datmat,dy)
    }

  } else if (eq==2) { # Dy(t)=y(t-1)+....

    if (p==0)  {
      y=temp[,1,drop=F]
      y.L1=temp[,2,drop=F]
      Dy=y-y.L1
      datmat=cbind(y=Dy, y.L1=y.L1)
      colnames(datmat)=c("dy","y.L1")
    } else {

      y=temp[,1,drop=F]
      y.L1=temp[,2,drop=F]
      Dy=y-y.L1
      dy=embed(Dy,p+1)[,-1,drop=F]
      colnames(dy)=paste("dy.L", 1:p, sep = "")

      y=embed(Dy,p+1)[,1,drop=F]
      y.L1=embed(y.L1,p+1)[,1,drop=F]
      datmat=cbind(y=y, y.L1=y.L1)
      colnames(datmat)=c("dy","y.L1")
      datmat=cbind(datmat,dy) }
  }

  if (is.null(x)) { datmat=datmat

  } else {
    x=as.matrix(x)
    X=embed(x,p+2)[,-seq(ncol(x)*(p+1)),drop=F]
    datmat=cbind(datmat,X)
  }

  as.data.frame(datmat)
}

.plag <-function(y,x=NULL,pmax,ic=c("AIC","BIC"),eq=1){
  y=as.matrix(y)

  if(ncol(y)>1) {stop("\nError: Only one series is allowed.")}

  aicp=bicp=rep(1,(pmax+1))
  datap=list()
  for (i in c(0,seq(pmax))) {

    if (is.null(x)) {datmat <- .datmat4ADF(y,x=NULL,p=i,eq)

    } else {datmat <- .datmat4ADF(y,x,p=i,eq)

    }
    datap[[i+1]]=datmat
    e1=as.matrix(resid(lm(datmat)))
    aicp[i+1] = log((t(e1)%*%e1)/nrow(e1)) + 2 * (ncol(datmat)-1 + 2)/nrow(e1)
    bicp[i+1]=log((t(e1)%*%e1)/nrow(e1)) +  (ncol(datmat)-1 + 2) * log(nrow(e1))/nrow(e1)

  }

  if(ic=="AIC") {p=which.min(aicp)-1} else {p=which.min(bicp)-1}

  return(list(p=p,data=datap[[p+1]]))

}




if (FALSE) {

  if (lam[1] <= 0.15) {mat_cv1  =mat_cv1[1,]
  } else if (0.15 < lam[1]  & lam[1] <= 0.25) {
    mat_cv1 = mat_cv1[2,]
  } else if (0.25 < lam[1] & lam[1] <= 0.35) {
    mat_cv1 = mat_cv1[3,]
  } else if (0.35 < lam[1] & lam[1] <= 0.45) {
    mat_cv1 = mat_cv1[4,]
  } else if (0.45 < lam[1] & lam[1] <= 0.55) {
    mat_cv1 = mat_cv1[5,]
  } else if (0.55 < lam[1] & lam[1] <= 0.65) {
    mat_cv1 = mat_cv1[6,]
  } else if (0.65 < lam[1] & lam[1] <= 0.75) {
    mat_cv1 = mat_cv1[7,]
  } else if (0.75 < lam[1] & lam[1] < 1) {
    mat_cv1 = mat_cv1[8,]}



  if (lam[2] <= 0.25) {
    mat_cv1 = mat_cv1[1]
  } else if (0.25 < lam[2] & lam[2] <= 0.35){
    mat_cv1 = mat_cv1[2]
  } else if (0.35 < lam[2] & lam[2] <= 0.45){
    mat_cv1 = mat_cv1[3]
  } else if (0.45 < lam[2] &  lam[2] <= 0.55){
    mat_cv1 = mat_cv1[4]
  } else if (0.55 < lam[2] & lam[2] <= 0.65){
    mat_cv1 = mat_cv1[5]
  } else if (0.65 < lam[2] & lam[2] <= 0.75){
    mat_cv1 = mat_cv1[6]
  } else if (0.75 < lam[2] & lam[2] <= 0.85){
    mat_cv1 = mat_cv1[7]
  } else if (0.85 < lam[2] & lam[2] < 1){
    mat_cv1 = mat_cv1[8]}


  if (lam[1] <= 0.15) {
    mat_cv5 = mat_cv5[1,]
  } else if (0.15 < lam[1] & lam[1] <= 0.25) {
    mat_cv5 = mat_cv5[2,]
  } else if  (0.25 < lam[1] & lam[1] <= 0.35) {
    mat_cv5 = mat_cv5[3,]
  } else if  (0.35 < lam[1] & lam[1] <= 0.45) {
    mat_cv5 = mat_cv5[4,]
  } else if  (0.45 < lam[1] & lam[1] <= 0.55) {
    mat_cv5 = mat_cv5[5,]
  } else if  (0.55 < lam[1] & lam[1] <= 0.65) {
    mat_cv5 = mat_cv5[6,]
  } else if  (0.65 < lam[1] & lam[1] <= 0.75) {
    mat_cv5 = mat_cv5[7,]
  } else if  (0.75 < lam[1] & lam[1] < 1) {
    mat_cv5 = mat_cv5[8,]}



  if (lam[2] <= 0.25) {
    mat_cv5 = mat_cv5[1]
  } else if  (0.25 < lam[2] & lam[2] <= 0.35){
    mat_cv5 = mat_cv5[2];
  } else if  (0.35 < lam[2] & lam[2] <= 0.45){
    mat_cv5 = mat_cv5[3];
  } else if  (0.45 < lam[2] & lam[2] <= 0.55){
    mat_cv5 = mat_cv5[4];
  } else if  (0.55 < lam[2] & lam[2] <= 0.65){
    mat_cv5 = mat_cv5[5];
  } else if  (0.65 < lam[2] & lam[2] <= 0.75){
    mat_cv5 = mat_cv5[6];
  } else if  (0.75 < lam[2] & lam[2] <= 0.85){
    mat_cv5 = mat_cv5[7];
  } else if  (0.85 < lam[2] & lam[2] < 1){
    mat_cv5 = mat_cv5[8]}


  if (lam[1] <= 0.15) {
    mat_cv10 = mat_cv10[1,]
  } else if   (0.15 < lam[1] & lam[1] <= 0.25) {
    mat_cv10 = mat_cv10[2,]
  } else if   (0.25 < lam[1] & lam[1] <= 0.35) {
    mat_cv10 = mat_cv10[3,]
  } else if   (0.35 < lam[1] & lam[1] <= 0.45) {
    mat_cv10 = mat_cv10[4,]
  } else if   (0.45 < lam[1] & lam[1] <= 0.55) {
    mat_cv10 = mat_cv10[5,]
  } else if   (0.55 < lam[1] & lam[1] <= 0.65) {
    mat_cv10 = mat_cv10[6,]
  } else if   (0.65 < lam[1] & lam[1] <= 0.75) {
    mat_cv10 = mat_cv10[7,]
  } else if   (0.75 < lam[1] & lam[1] < 1) {
    mat_cv10 = mat_cv10[8, ]}


  if (lam[2] <= 0.25) {
    mat_cv10 = mat_cv10[1]
  } else if   (0.25 < lam[2] & lam[2] <= 0.35) {
    mat_cv10 = mat_cv10[2]
  } else if   (0.35 < lam[2] & lam[2] <= 0.45) {
    mat_cv10 = mat_cv10[3]
  } else if   (0.45 < lam[2] & lam[2] <= 0.55) {
    mat_cv10 = mat_cv10[4]
  } else if   (0.55 < lam[2] & lam[2] <= 0.65) {
    mat_cv10 = mat_cv10[5]
  } else if   (0.65 < lam[2] & lam[2] <= 0.75) {
    mat_cv10 = mat_cv10[6]
  } else if   (0.75 < lam[2] & lam[2] <= 0.85) {
    mat_cv10 = mat_cv10[7]
  } else if   (0.85 < lam[2] & lam[2] < 1) {
    mat_cv10 = mat_cv10[8]}

}



