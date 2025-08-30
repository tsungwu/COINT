\name{pp}

\alias{pp}

\encoding{latin1}

\title{Phillips and Perron Unit Root Test}
\description{
  Implement the Phillips-Perron (1988) PP unit root test, including both Za (Z-alpha) and Zt (Z-tau) statistics.
  This wrapper allows inputting additional deterministic part, for example season dummies, but the asymptotic critical values are not available.
}
\usage{

pp(y,type=c("none","const","trend"),d=NULL,lags=c("short","long","nill"),use=c("nw","ba"))

}
\arguments{
  \item{y}{Vector to be tested for a unit root.}
  \item{type}{The deterministic parts in the test regression.}
  \item{d}{Additional deterministic parts in addition to "type"" in the test regression.}
  \item{lags}{Lags used for correction of error term. See section "details" below.}
  \item{use}{User specified lags for correction of error term. See section "details" below.The default is the lag determined by Newey-West  bandwidth ("nw")  with Bartlett  kernel ("ba").}

}

\details{
  \code{lags="short"} sets the number of lags to
  \eqn{\sqrt[4]{4 \times (n/100)}}, whereas
  \code{lags="long"} sets the number of lags to
  \eqn{\sqrt[4]{12 \times (n/100)}}. If lags="nil" is choosen, no error correction is made.\cr

  Furthermore, "lags" and "use" are mutually exclusive: As long as "use" is not set to be NULL, its argument will be chosen first.
  One can specify a different number of maximum lags by setting "use" accordingly. Users can input number of your souce.
  This version suports two bandwidth functions: "nw" for Newey-West and "and" for Andrews. The kernel functions are supported: "ba"=Bartlett, "pa"=Parzen, "qs"=Qudratic Spectral
}

\value{

  \item{Zt}{The Z-tau test statistic.}
  \item{cvZt}{Critical values of Zt.}
  \item{Za}{The Z-alpha test statistic.}
  \item{cvZa}{Critical values of Za.}
  \item{lag}{Number of lags used for kernel function.}
  \item{resid}{Regression residuals.}


}

\references{
Phillips, P.C.B. and Perron, P. (1988) Testing for a unit root in  time series regression. \emph{Biometrika},  75(2), 335-346.\cr
MacKinnon, J.G. (1991) Critical Values for Cointegration Tests- Long-Run Economic Relationships, eds. R.F. Engle and C.W.J. Granger, London, Oxford, 267-276.

}

\examples{
data(macro)
y=macro[,"INF"]
pp(y,
   type=c("none","const","trend")[3],
   lags = c("short", "long", "nil")[2],
   use=c("nw","ba")) # If argument "use" is NOT NULL, argument lags will be ignored.

pp(y,lags = c("short", "long", "nil")[2],
   type=c("none","const","trend")[3],
   use=NULL)

pp(y,lags = c("short", "long", "nil")[2],
   type=c("none","const","trend")[3],
   use=18)
}


\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{unit root test}

\note{
This code modifies function \code{ur.pp} of package \code{urca}, which does not have relevant critical values for "Z-alpha" test statistic.
}
