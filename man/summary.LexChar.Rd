\name{summary.LexChar}
\encoding{UTF-8}
\alias{summary.LexChar}

\title{Summary LexChar object}

\description{
Summarizes LexChar objects
}

\usage{
\method{summary}{LexChar}(object, CharWord=TRUE, stats=TRUE, CharDoc=TRUE, Vocab=TRUE,
    file = NULL, ...)

}

\arguments{
  \item{object}{object of TextData class}
  \item{CharWord}{if TRUE characteristic words of all the documents are shown (by default TRUE)}
  \item{stats}{if TRUE association statistics of lexical table are shown (by default TRUE)}
\item{CharDoc}{if TRUE characteristic source-documents of all the aggregate-documents are shown (by default TRUE)}
\item{Vocab}{if TRUE characteristic quantitative and qualitative variables of the words. CharWord and stats are provide}

  \item{file}{a connection, or a character string naming the file to print to in csv format. If NULL (the default), the results are not printed in a file}
  \item{\dots}{further arguments passed to or from other methods,...}
}

\details{
Vocab$quali$CharWord provides the qualitative variables and their categories.
Vocab$quali$stats provides association statistics for vocabulary and qualitative variables.
Vocab$quanti$CharWord provides characteristic quantitative variables for each word. This summary.LexChart function provides the characteristic words for each quantitative variable.
Vocab$quali$stats provides statistics for vocabulary and quantitative variables.
}


\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica Bécue-Bertaut, Josep-Antón Sánchez-Espigares}

\seealso{ \code{\link{LexChar}}, \code{\link{print.LexChar}}, \code{\link{plot.LexChar}}}




\examples{
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), var.agg="Gen_Edu", Fmin=10, Dmin=10, 
        remov.number=TRUE, stop.word.tm=TRUE)
res.LexChar <-LexChar(res.TD)
summary(res.LexChar)

}

\keyword{summary}

