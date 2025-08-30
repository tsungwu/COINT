\name{ccrQ}

\alias{ccrQ}

\encoding{latin1}

\title{Canonical Cointegrating Regression with Time Polynomial}

\description{
  Computes Park's (1992) canonical cointegrating regression estimator for cointegrated regressions with time polynomial, using OLS for the first stage regression.
  }


\usage{
ccrQ(y,
      x,
      type=c("trend","all"),
      v=15,
      q=2,
      ker_fun="parzen",
      aband=0,
      filter=0)
}

\arguments{
  \item{y}{The data of dependent variable in a regression.}
  \item{x}{The data of independent variables in a regression.}

  \item{type}{The deterministic parts in the regression.\cr
  const = intercept;\cr
  trend = lnear trend;\cr
  season = seasonal dummies, this option requires \code{ts} format;\cr
  all = all three options, which requires \code{ts} format.
  }

  \item{v}{Number of autocovariance terms to compute the spectrum at frequency zero, default=15.}
  \item{q}{degree of time polynomial, default=2.}
  \item{ker_fun}{Set kernel function to one of the available kernels, default="parzen". See section \code{details} below.}
  \item{aband}{Whether to activate the automatic bandwidth selection.\cr
  aband=1. To activate.\cr
  aband=0. Do not activate.The default.
  }

  \item{filter}{Whether to use an AR(1) filter to compute the spectrum at frequency zero.\cr
  filter=1. To use.\cr
  filter=0. Do not use.The default.
  }

}

\details{
1. Available kernels. Technical details are referred to Brillinger (1981,P.55)\cr
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

These kernels are written for ccrQ procedures, technically different from those used in \code{pp} and \code{kpss} tests.

2. Andrews (1991) has developed data based (or automatic) bandwidth procedures for computing the spectrum. \code{COINT} implements these procedures for the Parzen, Bartlette, Tukey-Hamming, and the Quadratic-Spectral kernels. When aband is active, COINT ignores the value you specify for the band-width parameter and automatically substitutes the data-based value.\cr
3. The aim of the AR(1) filter is to flatten the spectrum of residual around the zero frequency, thereby making it easier to estimate the true spectrum by simple averaging of the periodogram.\cr
}

\value{

  \item{coefTable}{Coefficients table.}
  \item{vcov}{Variance-covariance matrix for the parameter estimates.}
  \item{sigma}{Standard error of the residuals.}
  \item{rss}{Residual sum of squares.}
  \item{fit}{The fitted values, or conditional mean, of the regression.}
  \item{resid}{Estimated residuals.}
}

\references{
Park, J. Y. (1992) Canonical Cointegrating Regressions. \emph{Econometrica}, 60, 119-144.
}


\examples{

data(macro)
y=macro[,1]
x=macro[,-1]
out=ccrQ(y,x,q=3,type=c("trend","all")[1],v=15,ker_fun="parzen")
out$coefTable
out$vcov
tail(out$fit)
tail(out$resid)
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{cointegration}
