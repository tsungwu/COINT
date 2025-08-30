\name{GHansen}

\alias{GHansen}

\encoding{latin1}

\title{Gregory-Hansen Test for Cointegration in Models with Regime Shifts}

\description{Conduct the cointegration analysis with regime shifts, proposed by Gregory-Hansen (1996A).
}

\usage{GHansen(y, x, model,trim=0.1, use=c("nw","ba"))}


\arguments{
  \item{y}{The data of dependent variable in a regression.}

  \item{x}{The data of independent variables in a regression.}

  \item{model}{Models for the structural breaks:\cr
  1 = C   (level shift);\cr
  2 = C/T (level shift with trend);\cr
  3 = C/S (regime shift);\cr
  4 = regime and trend shift.}

  \item{trim}{The trimming percentage. Default is 10 percent.}
  \item{use}{Conditions for \code{pp} test for the nonstationarity of regression residuals.}
}

\details{This function calculates three residual-based test for cointegration with regime shifts: ADF, and Za, Zt of \code{pp}.
Argument \code{use} is detailed by the example of \code{pp} documentation.
}

\value{
  \item{result}{Comprehensive results of three tests.}
  \item{teststat}{Time series of three sequential tests.}
  \item{test.reg.adf}{Final regression output for ADF test.}
  \item{test.reg.za}{Final regression output for Za test.}
  \item{test.reg.zt}{Final regression output for Zt test.}
}

\references{
Gregory, A.W. and Hansen, B. E. (1996A) Residual-based tests for cointegration in models with regime shifts.\emph{Journal of Econometrics}, 70, 99-126.\cr
Gregory, A.W. and Hansen, B. E.(1996B). Tests for Cointegration in Models with Regime and Trend Shifts. \emph{Oxford Bulletin Economics and Statistics}, 58(3), 555-560.
}


\examples{

data(macro)
y=macro[1:300,1]
x=macro[1:300,-1]
output=GHansen(y,x,model=1, use=c("nw","ba"))
output$result
summary(output$test.reg.adf)
head(output$teststat)

#Plotting

test.name=rownames(output$result)[1]
stat=output$teststat[,test.name]
CV=output$result[test.name,1:3]
bpoint=output$result[test.name,5]
main=paste(paste(unlist(strsplit(test.name,"_")),collapse = " "),"test")

plot(stat,main=main,ylab="",xlab="",ylim=range(c(max(stat)+3,min(stat)-1,CV)));grid()
  abline(h=CV[1],col="red")
  abline(h=CV[2],col="blue")
  abline(h=CV[3],col="seagreen")
  abline(v=as.POSIXct(time(y)[bpoint]),col="orange",lty=2)

#  legend(x=as.POSIXct("2010-01-01"), y=max(stat)+3, legend=c("1% cv" , "5% cv", "10% cv"),
#         col=c("red", "blue", "seagreen"),xjust=1, yjust=1, lty=1,
#         horiz=TRUE, cex=0.66, bty="n")



#plot(y,main=colnames(y),ylab="",xlab="");grid()
#abline(v=time(y[output$bpoint,]),col="orange",lty=2)


}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\keyword{cointegration}

