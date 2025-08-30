\name{sw}

\alias{sw}

\encoding{latin1}

\title{Stock-Watson Common Trends Statistic}

\description{

Computes Stock and Watson (1988) common trends statistic for the null hypothesis that the data is a noncointegrated system (after allowing for a p-th
order polynomial time trend).
  }

\usage{
sw(data,p,v=15,ker_fun="parzen",aband=0,filter=0)
}
\arguments{
  \item{data}{Matrix of k-time series variables.}
  \item{p}{Order of the time polynomial in the null hypothesis.}
  \item{v}{Number of autocovariance terms to compute the spectrum at frequency zero, default=15.}

  \item{ker_fun}{Set kernel function to one of the available kernels, default="parzen". See section \code{details} below.}

  \item{aband}{Whether to activate the automatic bandwidth selection.\cr
  aband=1. To activate.\cr
  aband=0. Do not activate. The default.
  }

  \item{filter}{Whether to use an AR(1) filter to compute the spectrum at frequency zero.\cr
  filter=1. To use.\cr
  filter=0. Do not use. The default.}

}


\details{
1. Available kernels. Technical details are referred to Brillinger (1981,P.55)\cr
  "parzen"=Parzen kernel\cr
  "fejer"=Bartlett kernel\cr
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

2. Andrews (1991) has developed data based (or automatic) bandwidth procedures for computing the spectrum. \code{COINT} implements these procedures for the Parzen, Bartlette, Tukey-Hamming, and the Quadratic-Spectral kernels. When aband is active, \code{COINT} ignores the value you specify for the band-width parameter and automatically substitutes the data-based value.\cr
3. The aim of the AR(1) filter is to flatten the spectrum of residual around the zero frequency, thereby making it easier to estimate the true spectrum by simple averaging of the periodogram.\cr
}

\value{

  \item{sw_stat}{Test statistic. Reject the null of a unit root if the SW statistic < critical value}
  \item{sw_cv}{Critical values.}

}


\references{
Stock, J. & M. K. Watson (1988) Testing for Common Trends. \emph{Journal of the American Statistical Association}, 83, 1097-1107.
}


\examples{

data(macro)
sw(macro,p=1,v=15)


}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{cointegration}
