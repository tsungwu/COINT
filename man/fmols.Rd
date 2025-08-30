\name{fmols}

\alias{fmols}

\encoding{latin1}

\title{Multivariate Fully-Modified OLS Estimator}

\description{

  Phillips' (1995) Fully-Modified OLS estimator for single equation and multivariate cointegrated regression models.

  }

\usage{
fmols(y,
      x,
      type=c("const","trend","season","all"),
      v=15,
      ker_fun="parzen",
      aband=0,
      filter=1)
}
\arguments{
  \item{y}{The data of dependent variable(s) in a regression.}
  \item{x}{The data of independent variables in a regression.}

  \item{type}{The deterministic parts in the regression.\cr
  const = intercept;\cr
  trend = lnear trend;\cr
  season = seasonal dummies, this option requires \code{ts} format;\cr
  all = all three options, which requires \code{ts} format.
  }

  \item{v}{Number of autocovariance terms to compute the spectrum at frequency zero, default=15.}

  \item{ker_fun}{Set kernel function to one of the available kernels, default="parzen". See section \code{details} below.}

  \item{aband}{Whether to activate the automatic bandwidth selection.\cr
  aband=1. To activate.\cr
  aband=0. Do not activate.The default.
  }

  \item{filter}{Whether to use an AR(1) filter to compute the spectrum at frequency zero.\cr
  filter=1. To use.The default.\cr
  filter=0. Do not use.

  }


}

\details{
1. This fmols allows both single equation and multivariate system of equations. The multvariate case is a system that many dependent variables to common Xs.
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

These kernels are written for fmols procedures, technically different from those used in \code{pp} and \code{kpss} tests.\cr

3. Andrews (1991) has developed data based (or automatic) bandwidth procedures for computing the spectrum. \code{COINT} implements these procedures for the Parzen, Bartlette, Tukey-Hamming, and the Quadratic-Spectral kernels. When aband is active, COINT ignores the value you specify for the band-width parameter and automatically substitutes the data-based value.\cr
4. The aim of the AR(1) filter is to flatten the spectrum of residual around the zero frequency, thereby making it easier to estimate the true spectrum by simple averaging of the periodogram.\cr
}

\value{

  \item{beta}{Coefficient estimates.}
  \item{stderr}{Standard error of the residuals.}
  \item{tstat}{t-statistics of parameter estimates.}
  \item{vcov}{Variance-covariance matrix for the parameter estimates.}
  \item{fit}{The fitted values, or conditional mean, of the regression.}
  \item{resid}{Estimated residuals.}

}


\references{
Phillips, P. C. B (1995) Fully Modified Least Squares and Vector Autoregression, \emph{Econometrica}, 63, 1023-1078.
}


\examples{

data(macro)
y=macro[,1:2]
x=macro[,3:4]
out=fmols(y,x,type=c("const","trend","season","all")[2],v=15,ker_fun="bartlett")
out$beta
out$stderr
out$tstat #t-ratio
out$vcov
tail(out$fit)
tail(out$resid)
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{cointegration}
