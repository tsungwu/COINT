\name{Za}

\alias{Za}

\encoding{latin1}

\title{Phillips' (1987) Za and Zt Test for Unit Root}

\description{

Compute Phillips' (1987) Za and Zt statistics for the null hypothesis that y has a unit root.
  }

\usage{Za(y,p=1,v=15,ker_fun="parzen",aband=0,filter=0)}

\arguments{
  \item{y}{The data of dependent variable in a regression.}
  \item{p}{Order of the time polynomial in the cointegrating regressio.\cr
  p =-1: No deterministic term in the cointegrating regression.\cr
  p = 0: For a constant term.\cr
  p = 1: For a constant term and trend. Default.\cr
  P > 2: For time polynomial.\cr
Critical values are available for p within [1,5].}

  \item{v}{Number of autocovariance terms to compute the spectrum at frequency zero, default=15.}

  \item{ker_fun}{Set kernel function to one of the available kernels, default="parzen". See section \code{details} below.}

  \item{aband}{Whether to activate the automatic bandwidth selection.\cr
  aband=1. To activate.\cr
  aband=0. Do not activate.The default.
  }

  \item{filter}{Whether to use an AR(1) filter to compute the spectrum at frequency zero.\cr
  filter=1. To use.\cr
  filter=0. Do not use. The default.}


}


\value{

  \item{alpha}{Estimate of the AR(1) coefficient.}
  \item{za}{Za statistic for the series under the null has a unit root. Reject the null hypothesis of a unit root if the test statistic < critical value.}
  \item{za_cv}{Critical values of Za.}
  \item{zt}{Zt statistic for the series under the null has a unit root. Reject the null hypothesis of a unit root if the test statistic < critical value.}
  \item{zt_cv}{Critical values of Zt.}

}


\references{
Ouliaris, S., J. Y. Park, and P. C. B. Phillips (1989) Testing for a Unit Root in the Presence of a Maintained Trend. Ch. 1 in Baldev Raj (ed.), \emph{Advances in Econometrics and Modelling}. Netherlands: Kluwer Academic Publishers. \cr
Phillips, P. C. B. (1987) Time Series Regression with a Unit Root. \emph{Econometrica}, 55, 277-301.
}


\examples{

data(macro)
y=macro[,1]
Za(y,p=1,v=10)

}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{unit root test}
