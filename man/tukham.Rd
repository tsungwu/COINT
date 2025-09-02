\name{tukham}

\alias{tukham}

\encoding{latin1}

\title{Tukey-Hamming Kernel for Consistent Estimate of Long-run Variance}

\description{
Compute the Tukey-Hamming window to obtain consistent estimate of long-run variance of multivariate time series.
}

\usage{
tukham(data,v)
}

\arguments{
  \item{data}{Data matrix for computing consistent long-run variance, normally, multivariate regression residuals.}
  \item{v}{Number of autocovariance terms in the kernel.}

}


\value{
Brillinger, David R. (1981) \emph{Time Series Data Analysis and Theory}. San Francisco, CA: Holden-Day.
}

\references{
Brillinger, David R. (1981) Time Series Data Analysis and Theory. San Francisco, CA: Holden-Day.
}

\examples{
data(macro)
e=apply(macro, 2, function(x) x-mean(x))
tukham(e,v=15)
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{kernel}
