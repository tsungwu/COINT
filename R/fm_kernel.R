##=== kernel functions
# Brillinger, David R. (1981) Time Series: Data Analysis and Theory.
# San Francisco, CA: Holden-Day.


# Computes the Fejer Bartlett window, Brillinger (1981, p.55)
bartlett <- function(data,v) {
  e=as.matrix(data)
  if (v > nrow(e)) {v = nrow(e)-1}

  if (v >= 1) { weights = rep(0,v) }

  A = 0
  for (i  in seq(v)) {
    f = abs(i)/(v+1)
    m = 1.0 - f
    t1  =  apply(e[-seq(i),,drop=F],2, function(x) x-mean(x))
    t2  =  apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x))
    A  =  A  +  m*(t(t1)%*%t2)
    weights[i] = m

  }
  colnames(A)=rownames(A)=colnames(e)
  return(A/nrow(e))
}

#Computes the Parzen window, Brillinger (1981, p.55)
parzen <-function(data,v){
  e=as.matrix(data)

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
  return(A/nrow(e))
}

#Procedure to compute the Tukey-Hanning window, Brillinger (1981, p.55)
tukhan <- function(data,v) {
  e=as.matrix(data)
  if (v > nrow(e)) {v = nrow(e)-1}

  if (v >= 1) { weights = rep(0,v) }

  A = 0

  for (i  in seq(v)) {
    m = 0.5*(1.0 + cos((pi*i)/(v+1)))
    t1  =  apply(e[-seq(i),,drop=F],2, function(x) x-mean(x))
    t2  =  apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x))
    A  =  A  +  m*(t(t1)%*%t2)
    weights[i] = m

  }
  colnames(A)=rownames(A)=colnames(e)
  return(A/nrow(e))
}


#Procedure to compute the Tukey-Hamming window, Brillinger (1981, p.55)
tukham <- function(data,v) {
  e=as.matrix(data)
  if (v > nrow(e)) {v = nrow(e)-1}
  if (v >= 1) { weights = rep(0,v) }

  q=0.23
  A = 0
  for (i  in seq(v)) {
    m = 1.0 - (2*q) + 2*q*cos((pi*i)/(v+1))
    t1  =  apply(e[-seq(i),,drop=F],2, function(x) x-mean(x))
    t2  =  apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x))
    A  =  A  +  m*(t(t1)%*%t2)
    weights[i] = m

  }
  colnames(A)=rownames(A)=colnames(e)
  return(A/nrow(e))
}


#Computes the chauchy window, Brillinger (1981, p.55)
cauchy <- function(data,v) {
  e=as.matrix(data)
  if (v > nrow(e)) {v = nrow(e)-1}

  if (v >= 1) { weights = rep(0,v) }

  A = 0
  for (i  in seq(v)) {
    f = i/(v+1)
    m = 1/(1+f^2)
    t1  =  apply(e[-seq(i),,drop=F],2, function(x) x-mean(x))
    t2  =  apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x))
    A  =  A  +  m*(t(t1)%*%t2)
    weights[i] = m

  }
  colnames(A)=rownames(A)=colnames(e)
  return(A/nrow(e))
}

#Computes the Andrews (1991) Quadratic-Spectral window
#Andrews, D. W. K. (1991) Heteroskedasticity and Autocorrelation Consistent
#Covariance Matrix Estimation," Econometrica, 59: 817-858.
qs <- function(data,v) {
  e=as.matrix(data)

  if (v > nrow(e)) {v = nrow(e)-1}

  if (v >= 1) { weights = rep(0,v) }

  jb = seq(nrow(e)-1)/(v+1)
  jband = jb*1.2*pi;
  kn = ((sin(jband)/jband - cos(jband))/(jband^2))*3

  A = 0
  for (i  in seq(nrow(e)-1)) {
    m = kn[i]
    t1  =  as.matrix(apply(e[-seq(i),,drop=F],2, function(x) x-mean(x)))
    t2  =  as.matrix(apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x)))

  if (ncol(e)>1 & min(dim(t1))==1 & min(dim(t2))==1) {
      A  =  A  +  m*(t1)%*%t(t2)}
  else {A  =  A  +  m*(t(t1)%*%t2)}
    weights[i] = m

  }
  colnames(A)=rownames(A)=colnames(e)
  return(A/nrow(e))
}

#Computes the Gauss-Weierstrass window, Brillinger (1981, p.55)
gw <- function(data,v) {
  e=as.matrix(data)
  if (v > nrow(e)) {v = nrow(e)-1}

  if (v >= 1) { weights = rep(0,v) }

  A = 0
  for (i  in seq(v)) {
    f = i/(v+1)
    m = exp(-0.5*(f^2))
    t1  =  apply(e[-seq(i),,drop=F],2, function(x) x-mean(x))
    t2  =  apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x))
    A  =  A  +  m*(t(t1)%*%t2)
    weights[i] = m

  }
  colnames(A)=rownames(A)=colnames(e)
  return(A/nrow(e))
}


#Computes the Dirichlet window, Brillinger (1981, p.55)
dchlet <- function(data,v) {
  e=as.matrix(data)
  if (v > nrow(e)) {v = nrow(e)-1}

  if (v >= 1) { weights = rep(0,v) }

  A = 0
  for (i  in seq(v)) {
    m = 1
    t1  =  apply(e[-seq(i),,drop=F],2, function(x) x-mean(x))
    t2  =  apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x))
    A  =  A  +  m*(t(t1)%*%t2)
    weights[i] = m

  }
  colnames(A)=rownames(A)=colnames(e)
  return(A/nrow(e))
}

#Computes the modified Dirichlet window, Brillinger (1981, p.55)
mdchlet <- function(data,v) {
  e=as.matrix(data)
  if (v > nrow(e)) {v = nrow(e)-1L}

  if (v >= 1L) { weights = rep(0,v) }

  A = 0
  for (i  in seq(v)) {
    if (i<(v-1L)) { m = 1L } else {m=0.5}
    t1  =  apply(e[-seq(i),,drop=F],2, function(x) x-mean(x))
    t2  =  apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x))
    A  =  A  +  m*(t(t1)%*%t2)
    weights[i] = m

  }
  colnames(A)=rownames(A)=colnames(e)
  return(A/nrow(e))
}

#Computes the Reisz window, Brillinger (1981, p.55)
reisz <- function(data,v) {
  e=as.matrix(data)
  if (v > nrow(e)) {v = nrow(e)-1L}

  if (v >= 1L) { weights = rep(0,v) }

  A = 0
  for (i  in seq(v)) {
    f = i/(v+1)
    m = 1L-f^2
    t1  =  apply(e[-seq(i),,drop=F],2, function(x) x-mean(x))
    t2  =  apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x))
    A  =  A  +  m*(t(t1)%*%t2)
    weights[i] = m

  }
  colnames(A)=rownames(A)=colnames(e)
  return(A/nrow(e))
}

#Computes the Bohman window, Brillinger (1981, p.55)
bohman <- function(data,v) {
  e=as.matrix(data)
  if (v > nrow(e)) {v = nrow(e)-1}

  if (v >= 1) { weights = rep(0,v) }

  A = 0
  for (i  in seq(v)) {
    f = i/(v+1)
    m = (1 - abs(f))*cos(pi*f) + sin(pi*abs(f))/pi
    t1  =  apply(e[-seq(i),,drop=F],2, function(x) x-mean(x))
    t2  =  apply(embed(e,i+1)[,-seq(ncol(e)*i),drop=F],2, function(x) x-mean(x))
    A  =  A  +  m*(t(t1)%*%t2)
    weights[i] = m

  }
  colnames(A)=rownames(A)=colnames(e)
  return(A/nrow(e))
}










