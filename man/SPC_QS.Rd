\name{SPC_QS}

\alias{SPC_QS}

\encoding{latin1}

\title{Quadratic Spectral Kernel for Consistent Estimate of Long-run Variance}

\description{
  Compute the Quadratic Spectral kernel of Sul, Phillips and Choi (2003) to obtains consistent estimate of long-run variance.
  }
\usage{

SPC_QS(e,v)

}
\arguments{
  \item{e}{data that needs to compute consistent long-run variance, normally, regression residuals.}
  \item{v}{Number of lag terms to use when computing the long-run variance.}

}


\value{
Return the consistent estimate of long-run variance.
}

\references{
Carrion-i-Silvestre, J. L. and Sanso, A. (2006) A guide to the computation of stationarity tests.\emph{Empirical Economics}, 31(2), 433-448.\cr
Carrion-i-Silvestre, J. L. and Sanso, A. (2007) The KPSS test with two structural breaks. \emph{Spanish Economic Review}, 9(2), 105-127.\cr
Sul, D., Phillips, P.C.B., and Choi, C.Y.(2005)  Prewhitening Bias in HAC Estimation. \emph{Oxford Bulletin of Economics and Statistics}, 67(4), 517-546.
}

\examples{
data(macro)
y=macro[,"INF"]
e=y-mean(y)
SPC_QS(e,v=15)
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{kernel}
