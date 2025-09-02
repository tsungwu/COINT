\name{Kurozumi_Bartlett}

\alias{Kurozumi_Bartlett}

\encoding{latin1}

\title{Bartlett Kernel for Consistent Estimate of Long-run Variance}
\description{
  Compute the Bartlett kernel proposed by Kurozumi (2002) to obtains consistent estimate of long-run variance.
  }

\usage{Kurozumi_Bartlett(e)}

\arguments{

  \item{e}{data that needs to compute consistent long-run variance, normally, regression residuals}
}


\value{
Return the consistent estimate of long-run variance, Bartlett kernel proposed by Kurozumi (2002).
}

\references{
Kurozumi, E. (2002) Testing for stationarity with a break. \emph{Journal of Econometrics},108(1), 105-127.
}

\examples{
data(macro)
y=macro[,"INF"]
e=y-mean(y)
Kurozumi_Bartlett(e)
}


\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{kernel}
