\name{fmQ}

\alias{fmQ}

\encoding{latin1}

\title{Fully-Modified OLS Estimator with Time Polynomial}

\description{
  Computes the Phillips-Hansen (1990) Fully-Modified estimator for cointegrated regressions with Time Polynomial, using OLS for the first stage regression.
  }

\usage{
fmQ(y,x,type=c("trend","all"),v=15,
    q=2, ker_fun="parzen",aband=0,filter=0,sb_start=0.15)
}

\arguments{
  \item{y}{The data of dependent variable in a regression.}
  \item{x}{The data of independent variables in a regression.}

  \item{type}{The deterministic parts in the regression.\cr
  trend = lnear trend;\cr
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
  filter=0. Do not use. The default.
  }

  \item{sb_start}{The percentage specifies the beginning of sub-sample for stability test, and the end sample is \emph{(1-sb_start)}.}

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

These kernels are written for fmQ procedures, technically different from those used in \code{pp} and \code{kpss} tests.

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
Andrews, D. W. K. (1991) Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.\emph{Econometrica}, 59: 817-858. \cr
Brillinger, David R. (1981) \emph{Time Series Data Analysis and Theory}. San Francisco, CA: Holden-Day.\cr
Hansen, B. E. (1992) Tests for Parameter Instability in Regressions with I(1) Processes. \emph{Journal of Business and Economic Statistics}, 10, 321-335.\cr
Phillips, P. C. B. and Hansen B. E.(1990) Statistical Inference in Instrumental Variables Regression with I(1) Processes. \emph{Review of Economic Studies}, 57, 99-125.\cr
}


\examples{

data(macro)
y=macro[,1]
x=macro[,-1]
out=fmQ(y,x,type=c("trend","all")[1],v=15,q=3,ker_fun="parzen")
out$coefTable
out$vcov
out$stests
tail(out$fit)
tail(out$resid)


}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{cointegration}
