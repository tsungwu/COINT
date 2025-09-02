\name{ZA_2br}

\alias{ZA_2br}

\encoding{latin1}

\title{Zivot-Andrews unit root test with unknown one structural break.}

\description{
This function implements Zivot-Andrews sequential ADF unit root test with one unknown  structural break.
Handling two outlier models: "Innovational outlier" and "Additive outlier".
}

\usage{
ZA_2br(y,
      model=c("intercept", "both"),
      pmax=8,
      ic=c("AIC","BIC"),
      fixed=TRUE,
      trim=0.1,
      eq=1,
      trace=TRUE,
      season=FALSE)
}

\arguments{

  \item{y}{Univariate time series data, a preferable format is \code{ts} when "season=TRUE".
  }


  \item{model}{
Modelling where the unknown structural change occurs.\cr
intercept= structural change occurs at the intercept.\cr
trend= structural change occurs at the trend.\cr
both= structural changes occur at both the intercept and the trend.\cr
}

  \item{pmax}{
The maximal lags that are either included in the test regression or lag selection searches its optimal lag via "ic".
}

  \item{ic}{
Information criteria, "AIC" or "BIC". The default is "AIC".}

  \item{fixed}{
Logical. If TURE, pmax is the fixed inputed lags, and the default is TRUE.\cr
If FALSE, pmax is the maximal lags where lag selection searches its optimal lag.
}

  \item{trim}{
The trimming percentage. Default is 10%.
}

  \item{eq}{
The type of dependent variable in ADF equation.\cr
1= original y.\cr
2= differenced y.
}

\item{trace}{
Logical. If TURE, then screen displays the sequential progress.}

\item{season}{
Logical. If TURE,then seasonal dummies will be included in the test regression, and y must be in \code{ts} format.
}

}


\details{

This code entends Zivot-Andrews (1992) sequential procedure to two unknown structural changes. Critical values are from Narayan and Popp (2010).

}


\value{

\item{teststat}{The ADF test statistic in the presence of two structural breaks.}


\item{cval}{The critical values that are tabulated in Narayan and Popp (2010).}

\item{p}{
The number of lags that are included in the test regression.}


\item{bpoint1}{
The first breaking point that corresponds to the teststat.}

\item{bpoint2}{
The second breaking point that corresponds to the teststat.}


\item{timeElapse}{Time elapsed for sequential search.}

}



\references{
Narayan, P. K. and Popp, S. (2010) A new unit root test with two structural breaks in level and slope at unknown time.
\emph{Journal of Applied Statistics},37, 1425-1438. \cr

Zivot,E. and Andrews, W.K. (1992),Further Evidence on the Great Crasch, the Oil-Price Shock, and the Unit-Root Hypothesis.
\emph{Journal of Business & Economic Statistics},10(3), 251-270.
}


\examples{
data(macro) # US macro data, 1967M1-2025M7
# It takes time
\donttest{
y=macro[1:200,"INF"]
za2=ZA_2br(y,
           ic=c("AIC","BIC")[2],
           pmax=8,
           fixed=TRUE,
           model=c("intercept","trend","both")[1],
           trim=0.1,
           eq=1,
           season=TRUE)
za2$timeElapse[3]/60
za2$teststat
za2$cval
y[za2$bpoint1,] #The first dated strictural change
y[za2$bpoint2,] #The second dated strictural change
}
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{unit root test}

\note{
This code modifies function \code{ur.za} of package \code{urca}. We add "season", "eq", "outlier",and "trim".
Specifically, "outlier" is crucial, "season" is left to advanced research.
}

