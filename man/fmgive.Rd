\name{fmgive}
\alias{fmgive}

\encoding{latin1}

\title{Fully-Modified GIVE Estimator}

\description{

Computes the Kitamura-Phillips (1997) Fully-Modified GIVE estimator for univariate and multivariate cointegrated regression models.
  }

\usage{fmgive(y,x,z,v,t,ker_fun="parzen",times=5) }

\arguments{
  \item{y}{The data of dependent variable(s) in a regression.}
  \item{x}{The data of independent variables in a regression.}

  \item{z}{Instruments}

  \item{v}{Number of autocovariance terms to compute the spectrum at frequency zero, default=15.}

  \item{t}{number of terms to use in the computation of the WT matrix.}


  \item{ker_fun}{Set kernel function to one of the available kernels, default="parzen". See section \code{details} below.}

  \item{times}{Number of iteration to compute GIVE residuals, default =5.}

}

\details{
  1. Like FMOLS, fmgmm allows both single equation and multivariate system of equations. The multvariate case is a system that many dependent variables to common Xs.\cr
  2. Available kernels. Technical details are referred to Brillinger (1981,P.55)\cr
    "parzen"=Parzen kernel\cr
    "bartlett"=Bartlett kernel\cr
    "dchlet"= Dirichlet kernel\cr
    "mdchlet"= Modified Dirichlet kernel\cr
    "tukham"=Tukey-Hamming kernel\cr
    "tukhan"=Tukey-Hanning kernel\cr
    "cauchy"=Cauchy kernel\cr
    "bohman"=Bohman kernel\cr
    "reisz"=Riesz,Bochner kernel\cr
    "gw"= Gauss-Weierstrass kernel\cr
    "qs"= Andrews (1991) Quadratic-Spectral\cr

  These kernels are written for FM procedures, technically different from those used in \code{pp} and \code{kpss} tests.\cr

  3. Andrews (1991) has developed data based (or automatic) bandwidth procedures for computing the spectrum. \code{COINT} implements these procedures for the Parzen, Bartlette, Tukey-Hamming, and the Quadratic-Spectral kernels. When aband is active, COINT ignores the value you specify for the band-width parameter and automatically substitutes the data-based value.\cr
  4. The aim of the AR(1) filter is to flatten the spectrum of residual around the zero frequency, thereby making it easier to estimate the true spectrum by simple averaging of the periodogram.\cr
}

\value{

  \item{beta}{Coefficient estimates.}
  \item{stderr}{Standard error of the residuals.}
  \item{tstat}{t-statistics of parameter estimates.}
  \item{vcov}{Variance-covariance matrix for the parameter estimates.}
  \item{lromega}{long-run variance-covariance matrix of residuals.}
  \item{s1}{The first statistic for testing validity of overidentifying restrictions.}
  \item{s2}{The second statistic for testing validity of overidentifying restrictions.}
  \item{pvalue}{The p-value for s1+s2.}
  \item{fit}{The fitted values, or conditional mean, of the regression.}
  \item{resid}{GMM residuals.}

}


\references{
Kitamura, Y. and P. C. B. Phillips (1997) Fully-Modified IV, GIVE and GMM Estimation with Possibly Nonstationary Regressors and Instruments. \emph{Journal of Econometrics}, 80, 85-123.
}


\examples{

data(macro)

y=macro[-1,c(1)]
x=macro[-1,c(2,4)]
z=as.matrix(na.omit(exp(diff(macro)))) #IV

out=fmgive(y,x,z,v=15,t=5,ker_fun="parzen")
out$beta
out$vcov
out$stderr
out$tstat #t-ratio
tail(out$fit)
tail(out$resid)

y=macro[-1,c(1,3)]
x=macro[-1,c(2,4)]
z=as.matrix(na.omit(exp(diff(macro)))) #IV

out=fmgive(y,x,z,v=15,t=5,ker_fun="parzen")
out$beta
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{cointegration}
