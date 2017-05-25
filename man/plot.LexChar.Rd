\name{plot.LexChar}
\encoding{latin1}
\alias{plot.LexChar}

\title{Plot LexChar objects}

\description{
Draws the characteristic and anti-characteristic words of documents from a LexChar object.
}

\usage{
\method{plot}{LexChar}(x, char.negat=TRUE, col.char.posit="blue", col.char.negat="red",
col.lines="black", theme=theme_bw(), text.size=12, numr=1, numc=2, top=NULL, 
max.posit=15, max.negat=15, ...)
}
\arguments{
  \item{x}{object of LexChar class}
  \item{char.negat}{if TRUE, the anti-characteristic words are plotted (by default TRUE)}
  \item{col.char.posit}{color for the characteristic words (by default "blue")}
  \item{col.char.negat}{color for the anti-characteristic words (by default "red")}
  \item{col.lines}{color for the lines of barplot (by default "black")}
  \item{theme}{used to modify the theme settings by ggplot2 package (by default theme_bw())}
  \item{text.size}{size of the font (by default 12)}
  \item{numr}{number of rows in each multiple graph (by default 1 row)}
  \item{numc}{number of columns in each multiple graph (by default 2 columns)}
  \item{top}{title of the graph (by default NULL)}
  \item{max.posit}{maximum number of characteristic words (by default 15)}
  \item{max.negat}{maximum number of anti-characteristic words (by default 15)}
  \item{\dots}{further arguments passed to or from other methods...}
}

\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica Bécue-Bertaut, Josep-Anton Sánchez-Espigares}

\seealso{ \code{\link{LexChar}}, \code{\link{print.LexChar}}}


\examples{
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), var.agg="Gen_Edu", Fmin=10, Dmin=10,
        remov.number=TRUE, stop.word.tm=TRUE)
LD<-LexChar(res.TD,maxCharDoc = 0)
plot(LD)
}

\keyword{plot}
