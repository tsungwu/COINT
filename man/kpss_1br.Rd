\name{kpss_1br}

\alias{kpss_1br}

\encoding{latin1}

\title{KPSS Unit Root Test with One Structural Break}

\description{
  Implement the Kurozumi (2002) sequential kpss test with one structural break.
}

\usage{
kpss_1br(y,
        lags = c("short", "long", "nil"),
        model=c("intercept","both"),
        use=c("nw","ba"),
        trim=0.1)
}

\arguments{

  \item{y}{Vector to be tested for a unit root.}
  \item{lags}{Lags used for correction of error term.}
  \item{model}{
Modelling where the unknown structural change occurs.\cr
intercept= structural change occurs at the intercept.\cr
both= structural changes occur at both the intercept and the trend.\cr
}
  \item{use}{User specified lags for correction of error term. The default is the lag determined by Newey-West bandwidth "nw" with Bartlett "ba" kernel. Users can input your own number.
This version suports two bandwidth functions: "nw" for Newey-West, "and" for Andrews. Three kernel functions are supported by both bandwidth functions: "ba"=Bartlett, "pa"=Parzen, "qs"=Qudratic Spectral}

  \item{trim}{
The trimming percentage. Default is 10%.
}

}


\details{
  \code{lags="short"} sets the number of lags to
  \eqn{\sqrt[4]{4 \times (n/100)}}, whereas
  \code{lags="long"} sets the number of lags to
  \eqn{\sqrt[4]{12 \times (n/100)}}. If \code{lags="nil"} is choosen, then no error correction is made. Furthermore, \code{lags} and \code{use} are mutually exclusive. As long as \code{use} is not set to be NULL, its argumenta will be chosen fisrt.
  One can specify a different number of maximum lags by setting \code{use} accordingly.
}


\value{

  \item{teststat}{The kpss test statistic with one structural break.}
  \item{cval}{Critical values.}
  \item{bpoint}{The breaking point that corresponds to the teststat.}
  \item{tstats}{The sequential KPSS test statistic.}
}

\references{
Kwiatkowski, D., Phillips, P.C.B., Schmidt, P. and Shin, Y. (1992) Testing the Null Hypothesis of Stationarity Against the Alternative of a Unit Root: How Sure Are We That Economic Time Series Have a Unit Root? \emph{Journal of Econometrics}, 54, 159-178.\cr
Kurozumi, E. (2002) Testing for stationarity with a break. \emph{Journal of Econometrics},108(1), 105-127.\cr
Phillips, P.C.B. and Sainan Jin (2002) The KPSS test with seasonal dummies. \emph{Economics Letters}, 77, 239-243.\cr

}


\examples{

data(macro)
y=macro[,"INF"]
KPSS1=kpss_1br(y,model=c("intercept","both")[2],use=c("nw","ba"))
KPSS1$teststat
KPSS1$cval
y[KPSS1$bpoint,]
#Plot
plot.ts(KPSS1$tstats,ylim=range(c(KPSS1$tstats,KPSS1$cval)));grid()
abline(h=KPSS1$cval[1],col="red")
abline(h=KPSS1$cval[2],col="blue")
abline(h=KPSS1$cval[3],col="green")
abline(v=KPSS1$bpoint,col="red",lty=2)

}


\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{unit root test}
