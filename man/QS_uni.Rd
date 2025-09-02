\name{QS_uni}

\alias{QS_uni}

\encoding{latin1}

\title{Quadratic Spectral Kernel for Consistent Estimate of Long-run Variance}
\description{
  Compute the QS kernel to obtain consistent estimate of long-run variance.
  }

\usage{QS_uni(e,v)}

\arguments{
  \item{e}{A univariate time series for computing consistent long-run variance, normally, regression residuals.}
  \item{v}{Number of lag terms used to compute the long-run variance.}

}


\value{

  Return the consistent estimate of long-run variance, that PP and KPSS tests require. This procedure handles single time series only.

}

\references{
Brillinger, David R. (1981) \emph{Time Series Data Analysis and Theory}. San Francisco, CA: Holden-Day.

}


\examples{
data(macro)
y=macro[,"INF"]
e=y-mean(y)
QS_uni(e,v=15)
}


\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{kernel}
