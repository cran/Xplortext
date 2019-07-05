\name{summary.LexCA}
\encoding{UTF-8}
\alias{summary.LexCA}

\title{Summary LexCA object}

\description{
Summarizes LexCA objects
}

\usage{
\method{summary}{LexCA}(object, ncp=5, nb.dec = 3, ndoc=10, nword=10, nseg=10, 
 nsup=10, metaDocs=FALSE, metaWords=FALSE, file = NULL, ...)

}

\arguments{
  \item{object}{object of LexCA class}
  \item{ncp}{number of dimensions to be printed (by default 5)}
  \item{nb.dec}{number of decimal digits to be printed (by default 3)}
   \item{ndoc}{number of documents whose coordinates are listed (by default 10).
   Use ndoc="ALL" to have the results for all the documents. Use 
   ndoc=0 or ndoc=NULL if the results for documents are not wanted.}
  \item{nword}{number of words whose coordinates are listed (by default 10). 
   Use nword="ALL" to have the results for all the words. Use
   nword=0 or nword=NULL if the results for words are not wanted} 
  \item{nseg}{number of repeated segments whose coordinates are listed (by default 10). 
   Use nseg="ALL" to have the results for all the segments. Use 
   nseg=0 or nseg=NULL if the results for segments are not wanted} 
 \item{nsup}{number of supplementary elements whose coordinates are listed (by default 10). 
   Use nsup="ALL" to have the results for all the elements. Use 
   nsup=0 or nsup=NULL if the results for the supplementary elements are not wanted}
  \item{metaDocs}{axis by axis, the highest contributive documents are listed, separately for negative-part and positive-part documents; these documents have been identified in LexCA, taking into account lmd value (by default FALSE)}
  \item{metaWords}{axis by axis, the highest contributive words are listed, separately for negative-part and positive-part words; these words have been identified in LexCA, taking into account lmw value (by default FALSE)}
\item{file}{a connection, or a character string naming the file to print to (csv format). If NULL (the default), the results are not printed in a file}
 \item{\dots}{further arguments passed from other methods}
}

\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica Bécue-Bertaut, Josep-Antón Sánchez-Espigares}
\seealso{ \code{\link{LexCA}}, \code{\link{print.LexCA}},  \code{\link{plot.LexCA}}}

\examples{
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), Fmin=10, Dmin=10, stop.word.tm=TRUE)
res.LexCA<-LexCA(res.TD, lmd=1, lmw=1)
summary(res.LexCA)
}

\keyword{summary}

