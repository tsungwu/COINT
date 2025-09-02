\name{kpss_2br}

\alias{kpss_2br}

\encoding{latin1}

\title{KPSS Unit Root Test with Two Structural Breaks}

\description{
  Implement the kpss unit root test with two unknown structural breaks. Carrion-i-Silvestre and Sanso (2007) extends Kurozumi (2002) to two breaks, and create critical values.
}

\usage{
kpss_2br(y, lags = c("short", "long", "nil"), model=1, use=c("nw","ba"),trace=TRUE)
}
\arguments{
  \item{y}{Vector to be tested for a unit root.}
  \item{lags}{Lags used for correction of error term.}
\item{model}{Modelling where the unknown structural change occurs.\cr
	1 = AAn (Level shift without trend);\cr
  2 = AA (Level shift with trend);\cr
  3 = BB (Trend shift);\cr
  4 = CC (Level and trend shift)\cr
  Default is 1.
}

  \item{use}{User specified lags for correction of error term. The default is the lag determined by Newey-West bandwidth "nw" with Bartlett "ba" kernel. Users can input your own number.
This version suports two bandwidth functions: "nw" for Newey-West, "and" for Andrews. Three kernel functions are supported by both bandwidth functions: "ba"=Bartlett, "pa"=Parzen, "qs"=Qudratic Spectral}

\item{trace}{
Logical. If TURE (default), then screen displays the sequential progress.
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

\item{teststat}{The test statistic.}
\item{cval}{The critical values that are tabulated in Carrion-i-Silvestre and Sanso, A. (2007).}
\item{bpoint1}{The first breaking point that corresponds to the teststat.}
\item{bpoint2}{The second breaking point that corresponds to the teststat.}


\item{timeElapse}{Time elapsed for sequential search.}

}
\references{
Carrion-i-Silvestre, J. L. and Sanso, A. (2006) A guide to the computation of stationarity tests. \emph{Empirical Economics}, 31(2), 433-448.\cr
Carrion-i-Silvestre, J. L. and Sanso, A. (2007) The KPSS test with two structural breaks,\emph{Spanish Economic Review}, 9(2), 105-127.\cr
Kwiatkowski, D., Phillips, P.C.B., Schmidt, P. and Shin, Y. (1992) Testing the Null Hypothesis of Stationarity Against the Alternative of a Unit Root: How Sure Are We That Economic Time Series Have a Unit Root? \emph{Journal of Econometrics}, 54, 159-178.\cr
Kurozumi, E. (2002) Testing for stationarity with a break. \emph{Journal of Econometrics},108(1), 105-127.\cr
}


\examples{

data(macro)

\donttest{
y=macro[1:200,"INF"]

KPSS2=kpss_2br(y,model=1,use=c("nw","ba"))
KPSS2$teststat
KPSS2$cval
y[KPSS2$bpoint,]
}


}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{unit root test}
