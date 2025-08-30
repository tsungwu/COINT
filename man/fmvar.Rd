\name{fmvar}

\alias{fmvar}

\encoding{latin1}

\title{Fully-Modified VAR Estimator}

\description{

Computes the Phillips' (1995) Fully-Modified" VAR estimator for cointegrated regressions, using OLS for the first stage regression.

  }

\usage{
fmvar(data,p=1,q=5,v=15,type=c("const","trend","season","all"),
                   ker_fun="parzen",aband=0,filter=0)
}
\arguments{
  \item{data}{The dependent variables for a VAR system.}
  \item{p}{The number of lags for dependent variables, as in a VAR(p).}
  \item{q}{The number of lagged innovation terms to include in the fitted FMVAR(p,q).}
  \item{v}{Number of autocovariance terms to compute the spectrum at frequency zero, default=15.}

  \item{type}{The deterministic parts in the regression. Please note that fmvar will "de-" before inclusion, for example, \code{trend} is to pre-detrend the data.\cr
  const = intercept;\cr
  trend = lnear trend;\cr
  season = seasonal dummies, this option requires \code{ts} format;\cr
  all = all three options, which requires \code{ts} format.
  }


  \item{ker_fun}{Set kernel function to one of the available kernels, default="parzen". See section \code{details} below.}

  \item{aband}{Whether to activate the automatic bandwidth selection.\cr
  aband=1. To activate.\cr
  aband=0. Do not activate.The default.
  }

  \item{filter}{Whether to use an AR(1) filter to compute the spectrum at frequency zero.\cr
  filter=1. To use.\cr
  filter=0. Do not use. The default.

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

These kernels are written for FM procedures, technically different from those used in \code{pp} and \code{kpss} tests.\cr

2. Andrews (1991) has developed data based (or automatic) bandwidth procedures for computing the spectrum. \code{COINT} implements these procedures for the Parzen, Bartlette, Tukey-Hamming, and the Quadratic-Spectral kernels. When aband is active, COINT ignores the value you specify for the band-width parameter and automatically substitutes the data-based value.\cr
3. The aim of the AR(1) filter is to flatten the spectrum of residual around the zero frequency, thereby making it easier to estimate the true spectrum by simple averaging of the periodogram.\cr
}

\value{

  \item{beta}{Coefficient estimates.}
  \item{stderr}{Standard error of the residuals.}
  \item{tstat}{t-statistics of parameter estimates.}
  \item{vcov}{Variance-covariance matrix for the parameter estimates.}
  \item{fit}{The fitted values, or conditional mean, of the regression.}
  \item{resid}{Estimated residuals.}
  \item{data}{The data used in \code{fmvar}}
  \item{type}{The type used in \code{fmvar}}
  \item{p}{The p argument used in \code{fmvar}}
  \item{q}{The q argument used in \code{fmvar}}

}


\references{
Phillips, P. C. B (1995) Fully Modified Least Squares and Vector Autoregression. \emph{Econometrica}, 63, 1023-1078.
}


\examples{

data(macro)
out=fmvar(macro,p=1,q=6,v=15,type="trend",ker_fun="parzen",aband=0,filter=0)
out$beta
out$stderr
out$tstat
out$vcov
tail(out$data)
tail(out$resid)

ID1=grep(rownames(out$beta),pattern="_dL")
ID2=grep(rownames(out$beta),pattern="_L")
ID3=rownames(out$beta)[-c(ID1,ID2)]

out$beta[ID1,]; #innovation terms
out$beta[ID2,]; #VAR(1)
out$beta[ID3,]; #deterministic parts


}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{cointegration}
