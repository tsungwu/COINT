\name{ZA_1br}

\alias{ZA_1br}

\encoding{latin1}

\title{
Zivot-Andrews unit root test with unknown one structural break.
}

\description{
This function implements Zivot-Andrews sequential ADF unit root test with unknown one structural break.
Handling two outlier models: "Innovational outlier" and "Additive outlier".
}

\usage{
ZA_1br(y,
      model=c("intercept", "trend", "both"),
      outlier=1,
      pmax=8,
      ic=c("AIC","BIC"),
      fixed=FALSE,
      trim=0.1,
      eq=1,
      season=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{

Univariate time series data, a preferable format is \code{ts} when season=TRUE.
}


\item{model}{
Modelling where the unknown structural change occurs.\cr
intercept= structural change occurs at the intercept.\cr
trend= structural change occurs at the trend.\cr
both= structural changes occur at both the intercept and the trend.\cr
}

\item{outlier}{
The outlier model. \cr
1=Innovational outlier model.\cr
2=Additive outlier model.\cr
}

\item{pmax}{
The maximal lags that are either included in the test regression or lag selection searches its optimal lag via "ic".
}

\item{ic}{
Information criteria, "AIC" or "BIC". The default is "AIC".
}

\item{fixed}{
Logical. If TURE, pmax is the fixed inputed lags.
If FALSE, pmax is the maximal lags where lag selection searches its optimal lag.
}

\item{trim}{
The trimming percentage. Default is 10%.
}
\item{eq}{
The type of dependent variable in ADF equation.\cr
1= original y.
2= differenced y.
}

\item{season}{
Logical. If TURE,then seasonal dummies will be included in the test regression.
}
}



\value{


\item{teststat}{The Zivot-Andrews test statistic, which is the \code{min(tstats)}.}

\item{cval}{The critical values that are tabulated in Zivot and Andrews(1992)}

\item{p}{The number of lags that are included in the test regression.}

\item{bpoint}{The breaking point that corresponds to the teststat.}

\item{tstats}{The sequential ADF test statistic.}

\item{testreg}{The \code{lm()} output that corresponds to teststat.}

\item{timeElapse}{Time elapsed for sequential search.}


}


\references{

Zivot,E. and Andrews, W.K. (1992) Further Evidence on the Great Crasch, the Oil-Price Shock, and the Unit-Root Hypothesis.
\emph{Journal of Business & Economic Statistics},10(3), 251-270.
}



\examples{
data(macro) #US inflation rate, 1958M1-2025M7
y=macro[1:200,"INF"]
za1=ZA_1br(y,
           ic=c("AIC","BIC")[2],
           outlier=1,
           pmax=8,
           fixed=TRUE,
           model=c("intercept","trend","both")[1],
           trim=0.01,
           eq=1,
           season=TRUE)
za1$timeElapse[3]
za1$teststat
za1$cval
y[za1$bpoint,]
za1$tstats
za1$p

#plotting
plot.ts(za1$tstats,ylim=range(c(za1$tstats,za1$cval)))
abline(h=za1$cval[1],col="red")
abline(h=za1$cval[2],col="blue")
abline(h=za1$cval[3],col="green")
abline(v=za1$bpoint,col="red",lty=2)


}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{unit root test}

\note{
This code modifies function \code{ur.za} of package \code{urca}. We add "season", "eq", "outlier",and "trim".
Specifically, "outlier" is crucial, "season" is left to advanced research.
}
