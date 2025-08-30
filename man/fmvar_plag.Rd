\name{fmvar_plag}

\alias{fmvar_plag}

\encoding{latin1}

\title{Select the q in a FMVAR(p,q) by Specific Criterion}

\description{

Select the q for a FMVAR(p,q) for cointegrated regressions.

  }

\usage{
fmvar_plag(data,
          p=1,
          lag.max=12,
          v=15,
          type=c("const","trend","season","all"),
          ker_fun="parzen",
          aband=0,
          filter=0)
}
\arguments{
  \item{data}{The dependent variables for a VAR system.}
  \item{p}{The number of lags for dependent variables, as in a standard VAR(p).}
  \item{lag.max}{The maxum number of lags used to search for optimal q in fmvar(p,q).}
  \item{v}{Number of autocovariance terms to compute the spectrum at frequency zero, default=15.}

  \item{type}{The deterministic parts in the regression.\cr
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

These kernels are written for fm.ols procedures, technically different from those used in \code{pp} and \code{kpss} tests.\cr

2. Andrews (1991) has developed data based (or automatic) bandwidth procedures for computing the spectrum. \code{COINT} implements these procedures for the Parzen, Bartlette, Tukey-Hamming, and the Quadratic-Spectral kernels. When aband is active, COINT ignores the value you specify for the band-width parameter and automatically substitutes the data-based value.\cr
3. The aim of the AR(1) filter is to flatten the spectrum of residual around the zero frequency, thereby making it easier to estimate the true spectrum by simple averaging of the periodogram.\cr
}

\value{

  \item{selection}{The selected order of q: the selected lagged innovation terms to include in the fitted FMVAR(p,q)}

  \item{criteria}{The matrix of all lags and the values of four criteria: "AIC(q)", "HQ(q)", "SC(q)", "FPE(q)". }

}


\references{
Phillips, P. C. B (1995) Fully Modified Least Squares and Vector Autoregression. \emph{Econometrica}, 63, 1023-1078.
}


\examples{

data(macro)

Q=fmvar_plag(macro,
            p=1,
            v=15,
            lag.max=16,type="trend",
            ker_fun="parzen")$selection[1]

out=fmvar(macro,p=1,q=Q,v=15,type="trend", ker_fun="parzen")

}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{cointegration}
