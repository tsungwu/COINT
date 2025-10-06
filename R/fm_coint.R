###======= CCR
ccr<-function(y,
              x,
              type=c("const","trend","season","all"),
              v=15,
              ker_fun="parzen",
              aband=0,
              filter=0){

  type=match.arg(type)
  y0=timeSeries::as.timeSeries(y)
  x0=timeSeries::as.timeSeries(x)

  y=as.matrix(y)
  x=as.matrix(x)

  nobs  = nrow(y)
  m     = ncol(x)
  n     = ncol(y)

  ## Deterministic parts

  D=.Dummies4FM(data=x0,type)
  d=D$d
  xd=D$xd
  z=D$z

  a     = solve(t(z)%*%z)%*%t(z)%*%y
  yres = y-(z%*%a)
  what  = cbind(yres[-1,],diff(xd))
  colnames(what)=NULL

  var   = t(what)%*%what/nrow(what)
  sig   = .lrvar(what,v, ker_fun,aband,filter)
  del   = t(.Delta(what,v,ker_fun,aband,filter))

  d12      = del[,(n+1):(n+m)]
  true_vec = solve(sig[(n+1):(n+m),(n+1):(n+m)])%*%sig[(n+1):(n+m),1:n]
  adj      = solve(var)%*%d12%*%a[1:m,] + rbind(matrix(rep(0,n*1),1,n),true_vec)

  yn     =   y[-1,,drop=F] - (what %*% adj)
  xn     =   x[-1,,drop=F] - (what %*% solve(var) %*% d12)

  if (type %in% c("season","all")) {

  xn=ts(xn,end=end(as.ts(x0)),frequency=frequency(as.ts(x0)))
  xn=timeSeries::as.timeSeries(xn)

  }
  ## Deterministic parts for yn and xn

  data=.Dummies4FM(data=xn,type)$z

  beta   =   solve(t(data)%*%data)%*%t(data)%*%yn

  fit=data%*%beta
  resid  =   y[-1,,drop=F] - data%*%beta
  colnames(resid)=colnames(fit)=colnames(y)

  rownames(resid)=rownames(y)[-1]
  rownames(fit)=rownames(y)[-1]


  rss    =   t(resid)%*%resid

  sigma2  =   sig[1:n,1:n] - t(sig[(n+1):(m+n),1:n])%*%true_vec
  if(ncol(sigma2)==1 & nrow(sigma2)==1) {sigma2=as.numeric(sigma2)}
  vcov   =   sigma2 * solve(t(data)%*%data)
  stderr=sqrt(diag(vcov))
  colnames(resid)=paste0("u_",colnames(y))
  tstat=as.matrix(beta/stderr)
  beta=as.matrix(beta)
  rownames(beta)=c(colnames(x),rownames(beta)[-seq(ncol(x))])

  colnames(tstat)="t value"
  colnames(beta)="Estimate"
  stderr=as.matrix(stderr)
  colnames(stderr)="Std. Error"


  return(list(coefTable=cbind(beta,stderr,tstat),
              vcov=vcov,
              sigma=sqrt(sigma2),
              rss=rss,
              fit=fit,
              resid=resid))

}


###======= CCRQ
ccrQ<-function(y,
               x,
               type=c("trend","all"),
               v=15,
               q=2,#degree of time polynomial
               ker_fun="parzen",
               aband=0,
               filter=0){
  type=match.arg(type)

  y0=timeSeries::as.timeSeries(y)
  x0=timeSeries::as.timeSeries(x)
  y=as.matrix(y)
  x=as.matrix(x)
  nobs  = nrow(y)
  m     = ncol(x)
  n     = ncol(y)

  ## Deterministic parts
  D=.Dummies4FMQ(data=x0,type,q)
  d=D$d
  xd=D$xd
  z=D$z

  a     = solve(t(z)%*%z)%*%t(z)%*%y
  yres = y-(z %*%a )
  what  = cbind(yres[-1,],diff(xd))
  colnames(what)=NULL

  var   = t(what) %*% what/nrow(what)
  sig   = .lrvar(what,v, ker_fun,aband,filter)
  del   = t(.Delta(what,v,ker_fun,aband,filter))

  d12      = del[,(n+1):(n+m)]
  true_vec = solve(sig[(n+1):(n+m),(n+1):(n+m)])%*%sig[(n+1):(n+m),1:n]
  adj      = solve(var)%*%d12%*%a[1:m,] + rbind(matrix(rep(0,n*1),1,n),true_vec)

  yn     =   y[-1,,drop=F] - (what %*% adj)
  xn     =   x[-1,,drop=F] - (what %*% solve(var) %*% d12)

  if (type == "all") {

    xn=ts(xn,end=end(as.ts(x0)),frequency=frequency(as.ts(x0)))
    xn=timeSeries::as.timeSeries(xn)

  }


  ## Deterministic parts for xn and yn
  data=.Dummies4FMQ(data=xn,type,q)$z

  beta   =   solve(t(data)%*%data)%*%t(data)%*%yn

  fit=data%*%beta
  resid  =   y[-1,,drop=F] - data%*%beta
  colnames(resid)=colnames(fit)=colnames(y)

    rownames(resid)=rownames(y)[-1]
    rownames(fit)=rownames(y)[-1]

  rss    =   t(resid)%*%resid

  sigma2  =   sig[1:n,1:n] - t(sig[(n+1):(m+n),1:n])%*%true_vec
  if(ncol(sigma2)==1 & nrow(sigma2)==1) {sigma2=as.numeric(sigma2)}
  vcov   =   sigma2 * solve(t(data)%*%data)

  stderr=as.matrix(sqrt(diag(vcov)))
  tstat=beta/stderr
  colnames(resid)=paste0("u_",colnames(y))
  colnames(tstat)="t value"
  colnames(beta)="Estimate"
  colnames(stderr)="Std. Error"


  return(list(coefTable=cbind(beta,stderr,tstat),
              vcov=vcov,
              sigma=sqrt(sigma2),
              rss=rss,
              fit=fit,
              resid=resid))

}



###======= Phillips and Hansen (1992) FMOLS
fm <- function(y,
               x,
               type=c("const","trend","season","all"),
               v=15,
               ker_fun="parzen",
               aband=0,
               filter=0,
               sb_start=0.15) {
  type=match.arg(type)
  y0=timeSeries::as.timeSeries(y)
  x0=timeSeries::as.timeSeries(x)
  y=as.matrix(y)
  x=as.matrix(x)
  t  = nrow(y)
  m     = ncol(x)
  n     = ncol(y)

  ## Deterministic parts
  D=.Dummies4FM(data=x0,type)
  d=D$d
  xd=D$xd
  z=D$z

  xxi   = solve(t(z)%*%z,tol=.Machine$double.eps^2)
  xy    = t(z)%*%y
  coef  = xxi%*%xy #OLS coefficients

  u = as.matrix(y - ( z%*% coef)) # OLS Residual


  ##=== Fully-modifying begins

  e = cbind(u[-1,],diff(xd))

  colnames(e)=c("u_ols",colnames(xd))

  del   = .Delta(e,v,ker_fun,aband,filter)
  sig   = .lrvar(e,v,ker_fun,aband,filter)

  sig11 = sig[seq(n),seq(n)]
  A=chol(sig[(n+1):(m+n),(n+1):(m+n)])
  s21   = chol2inv(A)%*%sig[(n+1):(m+n),seq(n)]
  rho   = (sig[seq(n),(n+1):(m+n)]%*%s21)%*%solve(sig11)
  s112  = sig11 - (sig[seq(n),(n+1):(m+n)])%*%s21

  ys    = y[-1,] - diff(xd)%*%s21
  del21 = nrow(e)*(del[(n+1):(m+n),seq(n)]-(del[(n+1):(m+n),(n+1):(m+n)]%*%s21))

  if (type != "none") {del21 = rbind(del21,matrix(rep(0,ncol(d)*n),ncol(d),n))}

  xk    = z[-1,]
  ixx   = solve(t(xk)%*%xk)
  beta  = ixx%*%((t(xk) %*% ys)-del21) #FM coefficients

  if(max(dim(s112))==1) {s112=as.numeric(s112)}
  vcov   = ixx*s112 # FM vcov
  stderr    =as.matrix(sqrt(diag(vcov)))

  fit=xk%*%beta
  colnames(fit)=colnames(y)

  rsd     = y[-1,] - fit #FM Residuals
  colnames(rsd)="resid"

  ## Hansen (1992) stability tests
  B=ys-(xk%*%beta)
  rownames(B)=NULL
  C0=rep(c(t(del21)/nrow(ys)),nrow(xk))

  C=matrix(C0,dim(xk)[1],dim(xk)[2],byrow = T)

  sc    = apply(xk,2,function(x) x*B) - C


  sc    = apply(sc, 2,cumsum)
  lc    = sum(diag((t(sc)%*%sc) %*% ixx)) / (s112*nrow(ys))

  t1    = round(t*sb_start)
  t2    = round(t*(1-sb_start))

  f     = rep(0,t2-t1+1)

  for (j in t1:t2){
    sj  = t(sc[j,])
    vj  = t(xk[1:j,])%*% xk[1:j,] #moment(xk[1:j,.],0)
    mj  = vj - vj %*%ixx %*%vj
    f[j-t1+1] = sj %*% chol2inv(chol(mj)) %*% t(sj)
  }


  fstat = as.matrix(f/s112)
  stests = cbind(Lc=lc,MeanF=colMeans(fstat),SupF=max(fstat))
  rownames(stests)="Hansen(1992)"


    rownames(rsd)=rownames(y)[-1]
    rownames(fit)=rownames(y)[-1]


  tstat=beta/stderr
  colnames(tstat)="t value"
  colnames(beta)="Estimate"
  colnames(stderr)="Std. Error"


  return(list(coefTable=cbind(beta,stderr,tstat),
              vcov=vcov,
              sigma=sqrt(s112),
              rss=t(rsd)%*%rsd,
              fit=fit,
              stests=stests,
              resid=rsd))
}


###======= Phillips and Hansen (1992) FM
fmQ <- function(y,
               x,
               type=c("trend","all"),
               v=15,
               q=2,#degree of time polynomial
               ker_fun="parzen",
               aband=0,
               filter=0,
               sb_start=0.15) {
  type=match.arg(type)

  y0=timeSeries::as.timeSeries(y)
  x0=timeSeries::as.timeSeries(x)

  y=as.matrix(y)
  x=as.matrix(x)


  t     =nrow(y) - 1
  m     = ncol(x)
  n     = ncol(y)

  ## Deterministic parts
  D=.Dummies4FMQ(data=x0,type,q)
  d=D$d
  xd=D$xd
  z=D$z

  xxi   = solve(t(z)%*%z,tol=.Machine$double.eps^2)
  xy    = t(z)%*%y
  coef  = xxi%*%xy #OLS coefficients

  u = as.matrix(y - ( z%*% coef)) # OLS Residual


  ##=== Fully-modifying begins here

  e = cbind(u[-1,],diff(xd))

  colnames(e)=c("u_ols",colnames(xd))

  del   = .Delta(e,v,ker_fun,aband,filter)
  sig   = .lrvar(e,v,ker_fun,aband,filter)

  sig11 = sig[seq(n),seq(n)]
  A=chol(sig[(n+1):(m+n),(n+1):(m+n)])
  s21   = chol2inv(A)%*%sig[(n+1):(m+n),seq(n)]
  rho   = (sig[seq(n),(n+1):(m+n)]%*%s21)%*%solve(sig11)
  s112  = sig11 - (sig[seq(n),(n+1):(m+n)])%*%s21

  ys    = y[-1,] - diff(xd)%*%s21
  del21 = nrow(e)*(del[(n+1):(m+n),seq(n)]-(del[(n+1):(m+n),(n+1):(m+n)]%*%s21))

  if (type != "none") {del21 = rbind(del21,matrix(rep(0,ncol(d)*n),ncol(d),n))}

  xk    = z[-1,]
  ixx   = solve(t(xk)%*%xk)
  beta  = ixx%*%((t(xk) %*% ys)-del21) #FM coefficients

  if(max(dim(s112))==1) {s112=as.numeric(s112)}
  vcov   = ixx*s112 # FM vcov
  stderr    =as.matrix(sqrt(diag(vcov)))

  fit=xk%*%beta
  colnames(fit)=colnames(y)

  rsd     = y[-1,] - fit #FM Residuals
  colnames(rsd)="resid"

  ## Hansen (1992) stability tests
  B=ys-(xk%*%beta)
  rownames(B)=NULL
  C0=rep(c(t(del21)/nrow(ys)),nrow(xk))

  C=matrix(C0,dim(xk)[1],dim(xk)[2],byrow = T)

  sc    = apply(xk,2,function(x) x*B) - C

  sc    = apply(sc, 2,cumsum)
  lc    = sum(diag((t(sc)%*%sc) %*% ixx)) / (s112*nrow(ys))

  t1    = round(t*sb_start)
  t2    = round(t*(1-sb_start))

  f     = rep(0,t2-t1+1)

  for (j in t1:t2){
    sj  = t(sc[j,])
    vj  = t(xk[1:j,])%*% xk[1:j,] #moment(xk[1:j,.],0)
    mj  = vj - vj %*%ixx %*%vj
    f[j-t1+1] = sj %*% chol2inv(chol(mj)) %*% t(sj)
  }


  fstat = as.matrix(f/s112)
  stests = cbind(Lc=lc,MeanF=colMeans(fstat),SupF=max(fstat))
  rownames(stests)="Hansen(1992)"


    rownames(rsd)=rownames(y)[-1]
    rownames(fit)=rownames(y)[-1]


  tstat=beta/stderr
  colnames(tstat)="t value"
  colnames(beta)="Estimate"
  colnames(stderr)="Std. Error"


  return(list(coefTable=cbind(beta,stderr,tstat),
              vcov=vcov,
              sigma=sqrt(s112),
              rss=t(rsd)%*%rsd,
              fit=fit,
              stests=stests,
              resid=rsd))
}



###======= Multivariate Fully-Modified OLS
fmols <-function(y,
                  x,
                  type=c("const","trend","season","all"),
                  v=15,
                  ker_fun="parzen",
                  aband=0,
                  filter=1){
  y0=y
  x0=x
  type=match.arg(type)
  y=as.matrix(y)
  x=as.matrix(x)
  m     = ncol(x)
  n     = ncol(y)

  ## Deterministic parts
  D=.Dummies4FM(data=x0,type)
  d=D$d
  xd=D$xd
  z=D$z

  ixx   = solve(t(z)%*%z)
  xy    = t(z)%*%y

  beta  = solve(t(z)%*%z) %*% t(z)%*%y

  u     = as.matrix(y - (z%*%beta))

  ##=== Begin Fully-Modifying

  xd=na.omit(diff(xd))
  rownames(xd)=seq(nrow(xd))

  e = cbind(u[-1,,drop=F],xd)
  e=as.matrix(e)

  del=.Delta(e, v, ker_fun, aband,filter) #long-run covariance matrix
  del=t(t(del)[-seq(n),])
  del=del[-c((nrow(del)-m+1):nrow(del)),,drop=F]

  sig = .lrvar(e, v, ker_fun, aband,filter)
  sig=t(sig[-seq(n),,drop=F])
  sigxx = sig[-seq(n),,drop=F]


  delxx=.Delta(xd, v, ker_fun, aband,filter)

  true  = del%*%solve(sigxx)

  ys    = as.matrix(y[-1,,drop=F] - xd%*%t(true))
  dels  = as.matrix(t(del)-delxx%*%t(true))

  if(type != "none") {
    dels=rbind(dels,matrix(rep(0,ncol(d)*n),ncol(d),n))

  }

  xk    = z[-1,,drop=F]
  ixx   = solve(t(xk)%*%xk)

  beta  = ixx%*%((t(xk)%*%ys)-(nrow(xk)*dels)) # Multivariate FM

  fit=xk%*%beta


  # Compute the co-variance matrix

  temp0=NULL
  for (i in seq(ncol(u))) {

    temp0=cbind(temp0,u[-1,i] * xk)

  }

  bige=.lrvar(temp0,v,ker_fun,aband,filter)

  temp=diag(1,ncol(u),ncol(u))  %x%  (sqrt(nrow(xk))*ixx)
  vcov   = temp %*% bige %*% temp

  resid=y[-1,,drop=F]-fit
  if(is.ts(y0)) {
    resid=ts(resid,end=end(y0),frequency=frequency(y0))
    fit=ts(fit,end=end(y0),frequency=frequency(y0))
  } else {
    rownames(resid)=rownames(fit)=rownames(y)[-1]
  }
  colnames(resid)=paste0("u_",colnames(y))
  stderr=matrix(sqrt(diag(vcov)),nrow(beta), ncol(beta))
  colnames(stderr)=colnames(beta)
  rownames(stderr)=rownames(beta)

  return(list(beta=beta,
              stderr=stderr,
              tstat=beta/stderr,
              vcov=vcov,
              fit=fit,
              resid=resid))
}



###======= Estimating Fully-Modified GMM
fmgmm<-function(y,
                x,
                z,
                v=15,
                ker_fun="parzen",
                times=5){
  y0=y
  x0=x
  y=as.matrix(y)
  x=as.matrix(x)
  z=as.matrix(z)
  n = ncol(y)
  m = ncol(x)
  nz = ncol(z)
  colnames(z)=paste0("IV",seq(nz))

  if (times < 0) {times = 1}

  #   Naive IV...  /

  ahatx = (t(y)%*%z)%*%solve(t(z)%*%z)%*%(t(z)%*%x)%*%solve((t(x)%*%z)%*%solve(t(z)%*%z)%*%(t(z)%*%x))

  #   Naive Residuals....

  uhat0 = y - x%*%t(ahatx)
  #   Compute sztl matrix using uhat0 and z ....
  #   Construct observation matrix of (u times z)

  tmp = (uhat0 %x% t(as.matrix(rep(1,ncol(z))))) * (t(as.matrix(rep(1,ncol(uhat0)))) %x% z)
  sigmal = match.fun(ker_fun)(tmp,v)$amat;
  sztl   = ((t(tmp)%*%tmp)/nrow(tmp)) + sigmal + t(sigmal)

  #   Compute GMM estimator  (Step 2)

  tmp     = (diag(1,n,n) %x% (t(x)%*%z)) %*% solve(sztl)


  vecagmm = solve(tmp %*% (diag(1, n, n) %x% (t(z) %*% x))) %*% tmp%*%(as.vector(t(t(y)%*%z)))

  #   Reshape the estimator matrix back to an (n x m) matrix */

  ahatgmm   = matrix(vecagmm, n, m,byrow=TRUE)

  #   Compute the GMM residuals

  uhat0gmm = y - x%*%t(ahatgmm)

  for (i in seq(times)) {#i=1

    #Recompute sztl using GMM residuals

    tmp = (uhat0gmm %x% t(as.matrix(rep(1,ncol(z))))) * (t(as.matrix(rep(1,ncol(uhat0gmm)))) %x% z)
    sigmal = match.fun(ker_fun)(tmp,v)$amat;
    sztl   = ((t(tmp)%*%tmp)/nrow(tmp)) + sigmal + t(sigmal)

    # Compute Omega and .Delta matrices
    rsd=uhat0gmm[-1,,drop=F]
    tmp1 = cbind(rsd,diff(x),diff(z))
    tmp2 = match.fun(ker_fun)(tmp1,v)$amat

    del   = (t(tmp1)%*%tmp1)/nrow(tmp1) + tmp2
    omg = del + t(tmp2)

    # Extract the required matrices
    if (n==1) {
    A=t(as.matrix((omg[seq(n),])))
    omega0a = t(as.matrix(A[,-seq(n)]))
    omegaaa = t(t(omg[-seq(n),])[-seq(n),])
    colnames(omegaaa)=rownames(omegaaa)=NULL
    } else {

      A=t(as.matrix((omg[seq(n),])))
      omega0a = t(A[-seq(n),])
      omegaaa = t(t(omg[-seq(n),])[-seq(n),])
      colnames(omegaaa)=rownames(omegaaa)=NULL

    }

    del0z   = del[1:n,(1+n+m):(m+n+nz)]
    deluz   =  t(t(del[-seq(n),])[-seq(n+m),])
    delp0z  = del0z - (omega0a %*% solve(omegaaa) %*%deluz)

    yplus = y[-1,,drop=F] - (diff(cbind(x,z),1) %*% t(omega0a %*% solve(omegaaa)))
    yt   = y[-1,,drop=F]
    zt   = z[-1,,drop=F]
    xt   = x[-1,,drop=F]

    #  Compute the estimator ...

    tmp    = (diag(n) %x% (t(xt)%*%zt))%*% solve(sztl)
    G=t(yplus)%*%zt-(nrow(y)*delp0z)
    vecagmm = solve(tmp%*%(diag(1,n,n)%x%(t(zt)%*%xt)))%*% tmp %*% as.matrix(as.numeric(t(G)))

    ahatgmm = matrix(vecagmm,n, m, byrow=T)
    vcov = solve(tmp%*%(diag(n) %x% (t(zt) %*% xt))) # Variance covariance matrix ; 2005/04/12 */

    # Recompute the residuals

    uhat0gmm = y - x %*% t(ahatgmm)

  }  # End of iteration loop --- Times

  #===== End of FM-GMM estimation


  #===== IV validity test

  #  Step 4:

  # Recompute Sztl

  # Construct observation matrix of (u times z)  */
  uhat0gmm=as.matrix(uhat0gmm)
  tmp = (uhat0gmm %x% t(as.matrix(rep(1,ncol(z)))) ) * (t(as.matrix(rep(1,ncol(uhat0gmm)))) %x% z)
  sigmal = match.fun(ker_fun)(tmp,v)$amat;
  sztl    = ((t(tmp)%*%tmp)/nrow(tmp)) + sigmal + t(sigmal)


  # Recompute Omega and .Delta matrices */ /*tmp*/
  rsd=uhat0gmm[-1,,drop=F]
  tmp1 = cbind(rsd,diff(x),diff(z))
  tmp2 = match.fun(ker_fun)(tmp1,v)$amat

  del = (t(tmp1) %*% tmp1)/nrow(tmp1) + tmp2
  omg = del + t(tmp2)

  if (n==1) {
    #A=t(omg[-((nrow(omg)-(m+nz)+1):nrow(omg)),])
    A=t(as.matrix((omg[seq(n),])))
    omega0a = t(as.matrix(A[,-seq(n)]))
    omegaaa = t(t(omg[-seq(n),])[-seq(n),])
    colnames(omegaaa)=rownames(omegaaa)=NULL
  } else {

    A=t(as.matrix((omg[seq(n),])))
    omega0a = t(A[-seq(n),])
    omegaaa = t(t(omg[-seq(n),])[-seq(n),])
    colnames(omegaaa)=rownames(omegaaa)=NULL

  }

  # Note: omegaaa will be singular if z = x (i.e. instruments are the x)
  # So don't compute next part if this is just happens to be the case.

  st = 0
  s1 = 0
  s2 = 0
  lromega = 0
  dof = n*(nz-m)

  if (det(omegaaa) != 0 & dof > 0) {

    uhat0gmm = uhat0gmm[-1,,drop=F]

    tempz = t(omega0a %*% solve(omegaaa))

    uhatgmmp = uhat0gmm - diff(cbind(x,z))%*% tempz

    #Compute long-run covariance matrix of the GMM residuals */

    tmp2     = match.fun(ker_fun)(uhatgmmp,v)$amat;
    lromega  = t(uhatgmmp)%*%uhatgmmp/nrow(uhatgmmp) + tmp2 + t(tmp2)

    # Step 5: Long-run scores... */

    del0z   = del[1:n,(1+n+m):(m+n+nz)]
    deluz   =  t(t(del[-seq(n),])[-seq(n+m),])
    delp0z  = del0z - (omega0a %*% solve(omegaaa) %*% deluz)

    # Construct obs matrix of (u times z)  :by Y.K. */

    tmp = (uhat0gmm %x% t(as.matrix(rep(1,ncol(zt)))) ) * (t(as.matrix(rep(1,ncol(uhat0gmm)))) %x% zt)
    sigmal = match.fun(ker_fun)(tmp,v)$amat;
    sztl    = ((t(tmp)%*%tmp)/nrow(tmp)) + sigmal + t(sigmal)

    lrscores = t(uhatgmmp)%*% zt - (nrow(y)*delp0z)

    s1 = nrow(yt)* t(as.matrix(c(t(del0z)))) %*% solve(sztl) %*% as.matrix(c(t(del0z)))

    s2 = t(c(t(lrscores))) %*% solve(lromega %x% (t(zt)%*%zt)) %*% c(t(lrscores))

    st = s1 + s2

    pvalue = dchisq(st,df=dof) *2

  }
  ahatgmm=t(ahatgmm)
  rownames(ahatgmm)=colnames(x)
  colnames(ahatgmm)=colnames(y)

  fitted=x %*% ahatgmm
  colnames(fitted)=colnames(y)

  if(is.ts(y0)) {fitted=ts(fitted,end=end(y0),frequency=frequency(y0))

  } else {rownames(fitted)=rownames(y)}

  uhat0gmm = y - fitted
  uhat0gmm=as.matrix(uhat0gmm)
  colnames(uhat0gmm)=paste0("u_",colnames(y))
  stderr=matrix(sqrt(diag(vcov)),nrow(ahatgmm), ncol(ahatgmm))
  tstats=ahatgmm/stderr
  rownames(stderr)=rownames(ahatgmm)=rownames(tstats)=colnames(x)
  colnames(stderr)=colnames(ahatgmm)=colnames(tstats)=colnames(y)
  colnames(vcov)=rownames(vcov)=c(t(outer(colnames(tstats),rownames(tstats),FUN = paste,sep="_")))

  return(list(beta=ahatgmm,
              stderr=stderr,
              tstats=tstats,
              vcov=vcov,
              lromega=lromega,
              s1=as.numeric(s1),
              s2=as.numeric(s2),
              pvalue=as.numeric(pvalue),
              fit=fitted,
              resid=uhat0gmm))

}


fmgive <- function(y,x,z,v,t,ker_fun="parzen",times=5) {
  y0=y
  x0=x
  y=as.matrix(y)
  x=as.matrix(x)
  z=as.matrix(z)
  colnames(z)=paste0("IV_",colnames(z))
  n = ncol(y)
  m = ncol(x)
  nz = ncol(z)

  if (times <= 0) { times = 1 }

  # naive instrumental variable estimation */

  ahatx = (t(y)%*%z)%*%solve(t(z)%*%z)%*%(t(z)%*%x)%*%solve(t(x)%*%z%*%solve(t(z)%*%z)%*%(t(z)%*%x))
  u=y-x%*%t(ahatx)
  what_t = .wt(u, t)

  # now transform the X, Y, Z matrices using the weight matrix, what_t */


  xstar = what_t%*%(diag(n) %x% x)
  zstar = what_t%*%(diag(n) %x% z)


  vecystr = what_t%*% as.matrix(c(y))

  # compute GIVE */

  vecagive = solve((t(xstar)%*%zstar)%*%solve(t(zstar)%*%zstar)%*%(t(zstar)%*%xstar))%*%(t(xstar)%*%zstar)%*%(solve(t(zstar)%*%zstar)%*%(t(zstar)%*%vecystr))

  agive = matrix(vecagive,n,m,byrow = TRUE)

  uhatgive = y - x%*%t(agive)

  # uhatgivn is the UM(unmodified)-GIVE residual (may be used in tests) */
  uhatgivn = uhatgive

  # set up a loop to compute the FM-GIVE estimator and then iterate */

  # Iteration begins
  for (iters in seq(times)) {

    # re-calculate GLS-transformation matrix  */

    what_t = .wt(uhatgive,t)

    # now transform the X, Y, Z matrices using the weight matrix, what_t*/

      xstar = what_t%*%(diag(n) %x% x)
      zstar = what_t%*%(diag(n) %x% z)


    vecystr = what_t%*%as.matrix(c(y))

    #        iters = iters + 1

    # compute Omega and Delta matrices */
    rsd=uhatgive[-1,]
    tmp1 = cbind(rsd,na.omit(diff(cbind(x,z))))
    tmp2 = match.fun(ker_fun)(tmp1,v)$amat

    del   = t(tmp1)%*%tmp1/nrow(tmp1) + tmp2
    omega = del + t(tmp2)

#    trimr(trimr(omegaX, 0, m + nz)',n,nz)'

    # extract the required matrices ... */
    omega1a=t(omega[-c((nrow(omega)-m-nz+1):nrow(omega)),,drop=F])
    omega0ux=t(omega1a[-c(seq(n),(nrow(omega)-nz+1):nrow(omega)),,drop=F])

    omega2a=t(omega[-c(seq(n),(nrow(omega)-nz+1):nrow(omega)),,drop=F])
    omega0xx=t(omega2a[-c(seq(n),(nrow(omega2a)-nz+1):nrow(omega2a)),,drop=F])



    del1=t(del[-c(seq(n),(nrow(del)-nz+1):nrow(del)),,drop=F])
    deluxz   = t(del1[-seq(n+m),,drop=F])

    #  Compute the FM-GIVE estimator ... */

    zt = z[-1,]
    tmp = solve(omega0xx)%*%((t(na.omit(diff(x)))%*%zt)-(nrow(zt)*deluxz))%*%solve(t(zt)%*%zt)
    tmp = (omega0ux %*% tmp)
    A=solve(t(xstar)%*%zstar %*% solve(t(zstar)%*%zstar)%*% t(zstar)%*%xstar)
    B=t(xstar)%*%zstar
    C=solve(t(zstar)%*%zstar) %*% t(zstar)%*% vecystr-as.matrix(c(t(tmp)))

    vecagive = A  %*% B %*% C

    agive = matrix(vecagive,n,m,byrow=TRUE)

    uhatgive = y - x%*%t(agive)

  }   # End of iteration loop --- _Times */

  #***************************************************************/
  #          End of FM-GIVE estimation   &                       */
  #          Beginning of IV validity test procedure             */
  #***************************************************************/

  what_t = .wt(uhatgive,t)

  # now transform the X, Y, Z matrices using the weight matrix, what_t */


    xstar = what_t%*%(diag(n) %x% x)
    zstar = what_t%*%(diag(n) %x% z)

  vecystr = what_t%*%as.matrix(c(y))
  # Recompute Omega and Delta matrices */
  rsd=uhatgive[-1,]
  tmp1 = cbind(rsd,na.omit(diff(cbind(x,z))))
  tmp2 = match.fun(ker_fun)(tmp1,v)$amat

  del   = t(tmp1)%*%tmp1/nrow(tmp1) + tmp2
  omega = del + t(tmp2)

  # extract the required matrices ... */

  omega1a=t(omega[-c((nrow(omega)-m-nz+1):nrow(omega)),,drop=F])
  omega0ux=t(omega1a[-c(seq(n),(nrow(omega)-nz+1):nrow(omega)),,drop=F])

  omega2a=t(omega[-c(seq(n),(nrow(omega)-nz+1):nrow(omega)),,drop=F])
  omega0xx=t(omega2a[-c(seq(n),(nrow(omega2a)-nz+1):nrow(omega2a)),,drop=F])


  # Compute the delta matrix */
  # Note: "weights[iu]" are the kernel weights. They are re-computed in */
  # KERNEL.SRC every time match.fun(ker_fun) is called. */

  del1=t(del[-c(seq(n),(nrow(del)-nz+1):nrow(del)),,drop=F])
  deluxz   = t(del1[-seq(n+m),,drop=F])

  tmp1 = cbind(uhatgive[-1,],na.omit(diff(cbind(z))))
  tmp2 = match.fun(ker_fun)(tmp1,v)$amat

  del   = t(tmp1)%*%tmp1/nrow(tmp1) + tmp2
  del2a=t(head(del, nrow(del)-nz))
  del0z    = t(del2a[-seq(n),])

  tmp = match.fun(ker_fun)(uhatgive,v)$amat
  nobs = nrow(z)

  for (s in c(0:(n-1))) { #s=0
    tmp = diff(zstar[(1+(s*nobs)):((s+1)*nobs),])
    if (s == 0) {
      tmp1 = rbind(rep(0,ncol(tmp)),tmp)
      # zeros compensate for the differencing */
    } else {
      tmp1 = rbind(tmp1,rep(0,ncol(tmp)),tmp)
    }

  }

  # calculate delta_u0-zdot star  */
  # uhatgive = naivegive residual */
  tmp  = t(tmp1)%*% (what_t %*% as.matrix(c((uhatgivn))))/nobs

  tmp2 = what_t %*% as.matrix(c(uhatgivn))
  tmp2 = t(matrix(tmp2,n,ncol(t(uhatgivn)),byrow = TRUE))


  weights=match.fun(ker_fun)(uhatgive,v)$weights

  for (iu in seq(v)) {#iu=0
    tmp2 = rbind(tmp2[-1,,drop=F],rep(0,ncol(uhatgivn)))
    tmp = tmp + weights[iu]* t(tmp1)%*% as.matrix(c(tmp2))/nobs
  }

  # calculate the corrected GIVE residual */

  tmp1 = rbind(matrix(rep(0,ncol(x)),1,ncol(x)),na.omit(diff(x)))

  uhatgivs= as.matrix(uhatgive-(tmp1%*%solve(omega0xx)%*%t(omega0ux)))

  colnames(uhatgivs)=paste0("u_",colnames(y))

  # long-run covariance matrix */

  tmp1 = match.fun(ker_fun)(uhatgivs,v)$amat

  lromega = (t(uhatgivs)%*%uhatgivs)/nrow(uhatgivs) + tmp1 + t(tmp1)

  # Now compute the variance-covariance matrix...  */

  pz = zstar%*%solve(t(zstar)%*%zstar)%*%t(zstar)
  pp = solve(t(xstar)%*% pz %*% xstar)%*%t(xstar)%*%pz

  vcov = pp %*% (lromega %x% diag(nrow(uhatgivs))) %*%t(pp)

  # Fully modified score matrix */

  t2 = t(uhatgivs)%*%z
  del0zpls = del0z - (omega0ux)%*%solve(omega0xx)%*%deluxz
  t2 = t2 - nobs*del0zpls
  s1 = (nobs^2) * t(tmp)%*%solve(t(zstar)%*%zstar)%*%tmp


  s2 = t(as.matrix(c(t(t2)))) %*% solve(lromega %x% (t(z)%*%z))%*%as.matrix(c(t(t2)))

  st  = s1 + s2

  if (n*(nz-m) > 0 ) {pvalue = dchisq(st, n*(nz-m))*2 } else {pvalue = 0 }

stderr=matrix(sqrt(diag(vcov)),nrow(agive), ncol(agive))
tstats=agive/stderr
agive=t(agive)
stderr=t(stderr)
tstats=t(tstats)
rownames(stderr)=rownames(agive)=rownames(tstats)=colnames(x)
colnames(stderr)=colnames(agive)=colnames(tstats)=colnames(y)
colnames(vcov)=rownames(vcov)=c(t(outer(colnames(tstats),rownames(tstats),FUN = paste,sep="_")))


fitted=x%*%agive
colnames(fitted)=colnames(y)

if(is.ts(y0)) {fitted=ts(fitted,end=end(y0),frequency=frequency(y0))

} else {rownames(fitted)=rownames(y)}

return(list(beta=agive,
            stderr=stderr,
            tstats=tstats,
            vcov=vcov,
            lromega=lromega,
            s1=as.numeric(s1),
            s2=as.numeric(s2),
            pvalue=as.numeric(pvalue),
            fit=fitted,
            resid=uhatgivs))
}



.wt <- function(u,t) {
  u=as.matrix(u)
  nobs = nrow(u)
  tmp = embed(u,t+1)[,-seq(ncol(u))]
  tmp  = rbind(matrix(rep(0,t*ncol(tmp)),t,ncol(tmp)),tmp)

  all_c = solve(t(tmp)%*%tmp)%*%(t(tmp)%*%u)
  res_e = u-(tmp%*%all_c)
  all_c = t(all_c)
  s_e   = chol(solve((t(res_e)%*% res_e)/nrow(res_e)))

  for (r in c(0,seq(t))) { #r=0

    z = matrix(rep(0,nobs*nobs),nobs,nobs)

    for (iz in c(0,seq(nobs-r-1))){
      z[iz+1,iz+1+r] = 1L
    }

    # instead of setting t initial values to zero, as in the KP paper, the value
    # at time = 1 is used */

    z[1,1:(r+1)] = rep(1,r+1)
    z = t(z)
    if (r == 0) {
      w = s_e %x% z

    } else {
      tw = s_e %*% all_c[1:ncol(u),(1+((r-1)*ncol(u))):(r*ncol(u))]
      w = w - (tw %x% z)
    }
  }
  return(w)
}





#Lag Selection for fmvar

fmvar_plag <-function(data,
                      p=1,
                      lag.max=12,
                      v=15,
                      type=c("const","trend","season","all"),
                      ker_fun="parzen",
                      aband=0,
                      filter=0){
  type=match.arg(type)
  if (type %in% c("season","all")) {

    if (!(timeSeries::is.timeSeries(data))) {

      stop("\nError: For season, time series data must be in timeSeries() format")}
  }

   data=timeSeries::as.timeSeries(data)


   criteria=NULL;for (i in seq(lag.max)) {#i=q



    out=fmvar(data=data,p=p,q=i,v=v,type=type,ker_fun=ker_fun,aband=aband,filter=filter)

    if(out$type=="const") {detint=1
    } else if(out$type=="trend") {detint=2
    } else if (out$type=="season" & timeSeries::is.timeSeries(data)) {detintK=1+frequency(as.ts(data))
    } else if (out$type=="all" & timeSeries::is.timeSeries(data)) {detint=2+frequency(as.ts(data))
    } else if (out$type=="none") {detint=0}

    K=ncol(data)

    resids=out$resid
    sample <- nrow(resids)
    nstar <- ncol(as.matrix(out$beta))
    sigma.det <- det(crossprod(resids)/sample)

    aic <- log(sigma.det) + (2/sample) * (i * K^2 + K * detint)
    hq <- log(sigma.det) + (2 * log(log(sample))/sample) * (i * K^2 + K * detint)
    sc <- log(sigma.det) + (log(sample)/sample) * (i * K^2 + K * detint)
    fpe <- ((sample + nstar)/(sample - nstar))^K * sigma.det

    criteria=rbind(criteria,cbind(aic,hq,sc,fpe))

  }

  rownames(criteria) <- paste(seq(1:lag.max))
  colnames(criteria) <- c("AIC(q)", "HQ(q)", "SC(q)", "FPE(q)")

  order <- apply(criteria, 2, which.min);order

  return(list(selection = order, criteria = criteria))
}

###======= Estimating Fully-Modified VAR
fmvar <- function(data,
                  p=1,
                  q=5,
                  v=15,
                  type=c("const","trend","season","all"),
                  ker_fun="parzen",
                  aband=0,
                  filter=0) {
  type=match.arg(type)
  if (q<0) {stop("q must not be less than 0.")}

  n=ncol(data)

  D=.data4FMVAR(data,p,q,type=type)
  z=D$z
  zp=D$zp
  y=D$y
  yl=D$yl
  xd=D$xd
  x=D$x

  if(type!="none") {d =D$d}

  ixx   = solve(t(z)%*%z, tol=.Machine$double.eps^2)
  xy    = t(z)%*%y

  ##=== OLS results
  beta  = ixx %*% xy

  u     = y - (z%*%beta)
  colnames(u)=paste0("u_", colnames(y))


  ##=== Begin Fully-Modifying
  m= n*p #Nu,ber of AR terms

  if (q >=1) {

    e=cbind(u,xd[,seq(n*p),drop=F])

  } else {

    xd    = rbind(rep(0,n)|diff(xd[,seq(n*p),drop=F],1))
    e     = cbind(u,xd)
  }


  del=.lrvar(e, v, ker_fun, aband,filter) #long-run covariance matrix
  del=del[,-seq(n)][seq(n),,drop=F] #n by n*p

  sig = .lrvar(e, v, ker_fun, aband,filter)
  sigxx = sig[-seq(n),-seq(n),drop=F] #remove residual, vcov of innovation terms_(t-1)
  #n*p by n*p matrix, dZ_(t-1)

  delxx=t(.Delta(xd[,seq(n*p),drop=F], v, ker_fun, aband,filter))# innovation terms_(t-1)

  true  = del%*%solve(sigxx, tol=.Machine$double.eps^2) #n by n*p matrix

  ys    = y - xd[,seq(n*p),drop=F]%*%t(true)
  dels  = -true%*%delxx #n by n*p


  if (q > 0) {

    a = cbind(t(y)%*%x[,seq(ncol(zp)),drop=F],t(ys)%*%yl - nrow(ys)*dels)

  } else {

    a = t(ys)%*%yl - (nrow(ys)*dels)
  }

  if (type=="none") { a = a

  } else {
    a = cbind(a,t(y)%*%d)

  }
  colnames(a)=colnames(ixx)

  beta = t(a%*%ixx)

  fit=z%*%beta
  u = y - (z%*%beta)
  colnames(u)=paste0("u_",colnames(y))



  #   Okay, compute the co-variance matrix....

  vcov   = ((t(u)%*%u)/nrow(u)) %x% ixx

  stderr=matrix(sqrt(diag(vcov)),nrow(beta), ncol(beta))
  colnames(stderr)=colnames(beta)
  rownames(stderr)=rownames(beta)

    rownames(u)=rownames(data)[-seq(max(p,q)+1)]
    rownames(y)=rownames(data)[-seq(max(p,q)+1)]
    rownames(yl)=rownames(data)[-seq(max(p,q)+1)]
    rownames(fit)=rownames(data)[-seq(max(p,q)+1)]
    rownames(z)=rownames(data)[-seq(max(p,q)+1)]
    rownames(zp)=rownames(data)[-seq(max(p,q)+1)]


  return(list(beta=beta,
              stderr=stderr,
              tstat=beta/stderr,
              vcov=vcov,
              fit=fit,
              resid=u,
              data=data,
              type=type,
              p=p,
              q=q))
}



## Forecast FMVAR
fmvar_forecast<-function(output,
                         n.ahead=6) {
  newdata=timeSeries::as.timeSeries(output$data)
  p=output$p
  q=output$q
  type=output$type

  fcst=NULL
  for (i in seq(n.ahead)) {

    z= .data4FMVAR(data=newdata,p,q,type)$z

    ypred=z[nrow(z),] %*% output$beta
    newdata=rbind(newdata,ypred)
    if (output$type %in% c("season","all")) {

      newdata=ts(newdata,end=end(as.ts(output$data)),
                 frequency=frequency(as.ts(output$data)))
      newdata=timeSeries::as.timeSeries(newdata) }

    fcst=rbind(fcst,ypred)
  }
  return(fcst)
}



## FM utilities

#==============================================

.data4FMVAR <- function(data,
                        p,# VAR(p)
                        q,# number of innovation terms
                        type=c("const","trend","season","all")) {
  type=match.arg(type)
  if (type %in% c("season", "all")) {

    if (!(timeSeries::is.timeSeries(data))) {

      stop("\nError: For season, time series data must be in timeSeries() format")
    }

  }

  data=timeSeries::as.timeSeries(data)

  n=ncol(data)

  zp=embed(na.omit(diff(data)),max(p,q)+1) # differenced terms

  colnames(zp)=c(outer(paste0(colnames(data),"_dL"),0:q,FUN=paste0))
  zp=zp[,-seq(n)]
  if (q>=1) {
    y=embed(data[-1,],max(p,q)+1)[,seq(n),drop=F]
    colnames(y)=colnames(data)
    yL=embed(data[-1,],max(p,q)+1)[,,drop=F]
    colnames(yL)=c(outer(paste0(colnames(data),"_L"),0:q,FUN=paste0))
    yl=yL[,-seq(n)][,seq(n*p)]
    x=cbind(zp,yl )#All lagged innovation terms and AR(p)

  } else if (q==0){

    y=embed(data[-1,],max(p,q)+1)[,seq(n),drop=F]
    colnames(y)=colnames(data)
    yL=embed(data[-1,],max(p,q)+1)[,,drop=F]
    colnames(yL)=c(outer(paste0(colnames(data),"_L"),0:q,FUN=paste0))
    yl=yL[,-seq(n)][,seq(n*p)]
    x=yl #All lagged innovation terms and AR(p)
  } #end of q condition

  #x is all RHS variables

  ## Deterministic parts

  if (type=="const") {
    d =rep(1,nrow(x))
    xd=apply(x,2,function(x) resid(lm(x~d))) #de-mean
    z     = data.frame(x,const=d)

  } else if (type=="trend") {

    d = cbind(trend=seq(nrow(x))/nrow(x),const=rep(1,nrow(x)))
    xd=apply(x,2,function(x) resid(lm(x~d))) #de-trend
    z     = data.frame(x,d)

  }  else if (type=="season") {

    sd=ts(x[,1],end=end(as.ts(data)),frequency=frequency(as.ts(data)))
    d=cbind(forecast::seasonaldummy(sd),const=rep(1,nrow(x)))

    xd=apply(x,2,function(x) resid(lm(x~d))) #de-season
    z     = data.frame(x,d)

  } else if (type=="all") {

    sd=ts(x[,1],end=end(as.ts(data)),frequency=frequency(as.ts(data)))

    d=cbind(forecast::seasonaldummy(sd),
            trend=seq(nrow(x))/nrow(x),
            const=rep(1,nrow(x)))
    xd=apply(x,2,function(x) resid(lm(x~d))) #de-ALL3
    z     = data.frame(x,d)

  }  else if (type=="none"){

    xd    = x
    z     = x

  }

  # x is all RHS variables
  # xd is de-meaned x
  # z is xd + d

  z=as.matrix(z)
  if (type=="none") {d=NULL}
  return(list(z=z,zp=zp,y=y,yl=yl,xd=xd,x=x,d=d))
}

#==============================================
.Dummies4FM <-function(data,
                       type=c("const","trend","season","all")) {

  type=match.arg(type)

 if (type %in% c("season","all")) {

    if (!(timeSeries::is.timeSeries(data))) {

      stop("\nError: For season, time series data must be in timeSeries() format")}
 }

  x=timeSeries::as.timeSeries(data)

  ## Deterministic parts
  if (type=="const") {
    d =rep(1,nrow(x))
    xd=apply(x,2,function(x) resid(lm(x~d))) #de-mean
    z     = data.frame(x,const=d)

  } else if (type=="trend") {

    d = cbind(trend=seq(nrow(x))/nrow(x),const=rep(1,nrow(x)))
    xd=apply(x,2,function(x) resid(lm(x~d))) #de-trend
    z     = data.frame(x,d)

  }  else if (type=="season") {

    sd=as.ts(data[,1])
    d=cbind(forecast::seasonaldummy(sd),const=rep(1,nrow(x)))

    xd=apply(x,2,function(x) resid(lm(x~d))) #de-season
    z     = data.frame(x,d)

  } else if (type=="all") {

    sd=as.ts(data[,1])

    d=cbind(forecast::seasonaldummy(sd),
            trend=seq(nrow(x))/nrow(x),
            const=rep(1,nrow(x)))
    xd=apply(x,2,function(x) resid(lm(x~d))) #de-ALL3
    z     = data.frame(x,d)
  }  else if (type=="none"){

    xd    = x
    z     = x
  }

  z=as.matrix(z)
  d=as.matrix(d)
  xd=as.matrix(xd)


  rownames(z)=rownames(d)=rownames(xd)=rownames(data)


  return(list(d=d,xd=xd,z=z))

}


.Dummies4FMQ <-function(data,type=c("trend","all"),q=2) {

  if (type == "all") {
    if (!(timeSeries::is.timeSeries(data))) {

      stop("\nError: For season, time series data must be in timeSeries() format")}
  }


  x=timeSeries::as.timeSeries(data)

  type=match.arg(type)

  if (type=="trend") {

    td=seq(nrow(x))/nrow(x)
    trend=NULL
    for (i in 1:q){
      trend=cbind(trend,td^i)
    }
    colnames(trend)=paste0("trend",seq(q))
    d = cbind(const=rep(1,nrow(x)),trend=trend)
    xd=apply(x,2,function(x) resid(lm(x~d))) #de-trend
    z     = data.frame(x,d)

  }  else if (type=="all") {

    td=seq(nrow(x))/nrow(x)
    trend=NULL
    for (i in 1:q){

      trend=cbind(trend,td^i)
    }
    colnames(trend)=paste0("trend",seq(q))

    sd=as.ts(data[,1])
    d=cbind(forecast::seasonaldummy(sd),
            const=rep(1,nrow(x)),
            trend=trend)
    xd=apply(x,2,function(x) resid(lm(x~d))) #de-ALL3
    z     = data.frame(x,d)
  }

  z=as.matrix(z)
  d=as.matrix(d)
  xd=as.matrix(xd)

  return(list(d=d,xd=xd,z=z))

}




