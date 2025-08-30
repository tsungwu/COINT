\name{CZa}

\alias{CZa}

\encoding{latin1}

\title{Phillips' (1987) Za and Zt test for cointegration}

\description{

Test the null hypothesis of no cointegration between y and x using Phillips' (1987) Za and Zt statistics and Phillips and Ouliaris (1990) limit theory.
  }

\usage{CZa(y,x,p=1,v=15)}

\arguments{
  \item{y}{The data of dependent variable in a regression.}
  \item{x}{The data of independent variables in a regression.}
  \item{p}{Order of the time polynomial in the cointegrating regressio. Critical values are available for p within [1,5].\cr
  p =-1: No deterministic term in the cointegrating regression.\cr
  p = 0: For a constant term.\cr
  p = 1: For a constant term and trend. Default.\cr
  P > 2: For time polynomial.
  }

  \item{v}{Number of autocovariance terms to compute the spectrum at frequency zero, default=15.}

}


\value{

  \item{alpha}{Estimate of the AR(1) coefficient.}
  \item{cza}{Za statistic for non-cointegration.Reject the null hypothesis of no cointegration if the Z statistic < critical value.}
  \item{cza_cv}{Critical values of cza.}
  \item{czt}{Zt statistic for non-cointegration.Reject the null hypothesis of no cointegration if the Z statistic < critical value.}
  \item{czt_cv}{Critical values of czt.}

}


\references{
Phillips, P. C. B. (1987) Time Series Regression with a Unit Root. \emph{Econometrica}, 55, 277-301.\cr
Phillips, P. C. B. and Ouliaris S. (1990) Asymptotic Properties of Residual Based Tests for Cointegration. \emph{Econometrica}, 58, 165-193.

}


\examples{

data(macro)
y=macro[,1]
x=macro[,-1]
CZa(y,x,p=1,v=10)


}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{cointegration}
