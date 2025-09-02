\name{gw}

\alias{gw}

\encoding{latin1}

\title{Gauss-Weierstrass Kernel for Consistent Estimate of Long-run Variance}
\description{
  Computes the Gauss-Weierstrass window to obtain consistent estimate of long-run variance of multivariate time series.
  }
\usage{gw(data,v)}

\arguments{
  \item{data}{Data matrix for computing consistent long-run variance, normally, multivariate regression residuals.}
  \item{v}{Number of autocovariance terms in the kernel.}
}


\value{
Return the consistent estimate of long-run variance. This procedure handles both multivariate and single time series, which is basically designed for "fmvar","fmgmm" and "fmols".
}

\references{
Brillinger, David R. (1981) \emph{Time Series Data Analysis and Theory}. San Francisco, CA: Holden-Day.
}

\examples{
data(macro)
e=apply(macro, 2, function(x) x-mean(x))
gw(e,v=15)
}


\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{kernel}
