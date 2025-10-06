\name{bartlett}
\alias{bartlett}

\encoding{latin1}

\title{Bartlett Kernel for Consistent Estimate of Long-run Variance}
\description{
  Compute the Bartlett kernel to obtain consistent estimate of long-run variance of multivariate time series.
  }

\usage{bartlett(data,v)}

\arguments{
  \item{data}{Data matrix for computing consistent long-run variance, normally, multivariate regression residuals.}
  \item{v}{Number of autocovariance terms in the kernel.}

}


\value{

  \item{amat}{Return the consistent estimate of long-run variance. This procedure handles both multivariate and single time series, which is basically designed for "fmvar","fmgmm" and "fmols".}
  \item{weights}{The weights vector, used by function \code{fmgive}}
}

\references{
Brillinger, David R. (1981) Time Series Data Analysis and Theory. San Francisco, CA: Holden-Day.
}

\examples{
data(macro)
e=apply(macro, 2, function(x) x-mean(x))
bartlett(e,v=15)
}


\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{kernel}
