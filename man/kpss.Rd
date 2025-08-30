\name{kpss}

\alias{kpss}

\encoding{latin1}

\title{KPSS Unit Root Test for the null of stationarity}
\description{
  Implement the KPSS unit root test for the null of I(0) stationarity.
  The test type as deterministic component is specified as x, see example below.
}
\usage{
kpss(y, x, lags = c("short", "long", "nil"), use=c("nw","ba"))
}
\arguments{
  \item{y}{Vector to be tested for a unit root.}
  \item{x}{data matrix for deterministic component. For example a vector of one for intercept, or trend. The default is "x=NULL", which is the same of a vector of one}
  \item{lags}{Lags used for correction of error term.}
  \item{use}{User specified lags for correction of error term. See section \code{details} below.The default is the lag determined by Newey-West  bandwidth ("nw") with Bartlettkernel ("ba").}
}

\details{
  \code{lags="short"} sets the number of lags to \eqn{\sqrt[4]{4 \times (n/100)}}, whereas lags="long" sets the number of lags to
  \eqn{\sqrt[4]{12 \times (n/100)}}. If lags="nil" is choosen, no error correction is made.\cr

  Furthermore, "lags" and "use" are mutually exclusive: As long as "use" is not NULL, its argument will be chosen first.
  One can specify a different number of maximum lags by setting "use" accordingly. Users can input number of your souce.
  This version suports two bandwidth functions: "nw" for Newey-West and "and" for Andrews. The kernel functions are supported: "ba"=Bartlett, "pa"=Parzen, "qs"=Qudratic Spectral
}

\value{

  \item{teststat}{The KPSS test statistic.}
  \item{cval}{Critical values.}
  \item{lag}{Number of lags used for kernel function.}
  \item{resid}{Regression residuals.}
}


\references{
Kwiatkowski, D., Phillips, P.C.B., Schmidt, P. and Shin, Y. (1992) Testing the Null Hypothesis of Stationarity Against the Alternative of a Unit Root: How Sure Are We That Economic Time Series Have a Unit Root? \emph{Journal of Econometrics}, 54,159-178.\cr
Phillips, P.C.B. and Sainan Jin (2002) The KPSS test with seasonal dummies. \emph{Economics Letters}, 77, 239-243.\cr
}


\examples{

data(macro)
y=macro[,"INF"]
const=rep(1,nrow(y))
trend=seq(nrow(y))/nrow(y)
D=cbind(const,trend) #seasonal dummies can be specified here
KPSS=kpss(y,x=D,lags = c("short", "long", "nil")[2],
          use=c("nw","ba")) # If argument use isn't NULL, the argument "lags" will be ignored.
KPSS$teststat
KPSS$cval
KPSS$lag

kpss(y,x=D,lags = c("short", "long", "nil")[2],use=15)

kpss(y,x=D,
    lags = c("short", "long", "nil")[2],
    use=NULL) #if "use=NULL", argument "lags"  will be chosen as input.



}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{unit root test}
