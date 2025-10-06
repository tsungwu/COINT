.lrvar <- function(data,v,ker_fun,aband=0,filter=1) {
  e=as.matrix(data)
  if (filter==1) {
  a=.xahat(e)$a
  new_e=.xahat(e)$new_e
  tmp = solve(diag(ncol(e)) - a)
  s  = (t(new_e)%*%new_e)/nrow(new_e)

  io=.kacf(new_e,v,ker_fun,aband)
  lr = t(tmp) %*% (s + io + t(io)) %*% tmp
  } else {

    s  = (t(e)%*%e)/nrow(e)
    io=.kacf(e,v,ker_fun,aband)
    lr =(s + io + t(io))
  }

  rownames(lr)=colnames(lr)=colnames(e)
  return(lr)
}


.Delta <- function(data,v,ker_fun,aband=0,filter=1) {
  e=as.matrix(data)
  if (filter==1) {
  a=.xahat(e)$a
  new_e=as.matrix(.xahat(e)$new_e)
  tmp = as.matrix(solve(diag(ncol(e)) - a))
  s  = (t(new_e)%*%new_e)/nrow(new_e)

  io=.kacf(data=new_e,v,ker_fun,aband)

  su = (t(e)%*%e)/(nrow(e))
  lr = su + t(tmp) %*%io%*%tmp + (t(tmp) %*% t(a)%*% su)
} else {
  io = .kacf(data=e,v,ker_fun,aband)
  lr = (t(e)%*%e)/(nrow(e)) + io
}

  rownames(lr)=colnames(lr)=colnames(e)
  return(t(lr))

}



.xahat<-function(data) {
  e=as.matrix(data)
  t1 = e[-1,,drop=F]
  t2 = embed(e,2)[,-seq(ncol(e)),drop=F]
  rownames(t1)=rownames(t2)=NULL
  tail(t1)
  a  = solve(t(t2)%*%t2)%*%(t(t2)%*%t1)
  new_e  = t1 - (t2%*%a)
  rownames(new_e)=NULL
  return(list(a=a,new_e=new_e))

}

.kacf <- function(data, v, ker_fun="parzen", aband=0) {
  e=as.matrix(data)

  if (aband==1) { # Automatic bandwidth enabled
    eb = embed(e,2)[,-seq(ncol(e)),drop=F]
    ef = e[-1,,drop=F]
    ae = as.matrix(colSums(eb*ef)/colSums(eb^2))
    if(ncol(e)>1) {
      ee = ef - t(apply(eb,1,function(x) x*t(ae)))
    } else {
      ee = ef - eb%*%t(ae)
    }
    se = colMeans(ee^2)
    ad = colSums((se/((1-ae)^2))^2)
    a1 = 4*colSums((ae*se/(((1-ae)^3)*(1+ae)))^2)/ad
    a2 = 4*colSums((ae*se/((1-ae)^4))^2)/ad
    nobs = nrow(e)
    if (ker_fun == "qs") {
      v = 1.3221*((a2*nobs)^0.2)-1
    }  else if (ker_fun == "parzen") {
      v = 2.6614*((a2*nobs)^0.2)-1
    } else if (ker_fun == "bartlett") {
      v = 1.1447*((a1*nobs)^0.333)-1
    } else if (ker_fun == "tukham") {

      v = 1.7462*((a2*nobs)^0.2)-1}

  }
  return(match.fun(ker_fun)(data=e,v)$amat)
}

.detrend <- function (data,p){
  #p=order of the time polynomial in the fitted regression
  #p=-1, type="none"
  #p=0, type="const"
  #p>=1, type="trend
data=as.matrix(data)
  nobs     = nrow(data)
  const    = rep(1,nobs)

  if (p < -1 ) {stop("Error: p<-1 is not allowed.")}

  if (p == -1) {return(data)

    } else {

    if (p ==0) {

      xmat = const

    } else {

    td = seq(nobs)/nobs
    trend=NULL
    for (i in seq(p)) {
      trend=cbind(trend,td^i)
    }
    xmat= cbind(const,trend)
  }

  invx     = solve(t(xmat)%*%xmat)
  beta     = invx %*% (t(xmat)%*%data)
  resid    = data - xmat%*%beta
  return(resid)
    }

}


.covarf<-function(e,v,ker_fun="parzen",aband=0,filter=1) {

  Value=.Delta(e,v,ker_fun,aband,filter)-((t(e)%*%e)/nrow(e))

  return(Value)
}




