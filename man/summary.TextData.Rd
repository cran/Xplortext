\name{summary.TextData}
\encoding{UTF-8}
\alias{summary.TextData}

\title{Summary of TextData objects}

\description{
Summarizes TextData objects.
}

\usage{
\method{summary}{TextData}(object, ndoc=10, nword=50, nseg=50, ordFreq = TRUE, file = NULL, sep=";", 
   ...) 
}
\arguments{
  \item{object}{object of TextData class}
  \item{ndoc}{statistical report on the first ndoc documents (by default 10). 
   Use ndoc="ALL" to have the results for all the documents. Use
   ndoc=0 or ndoc=NULL if the results on the documents are not wanted } 
  \item{nword}{index of the nword first words (by default 50). 
   Use nword="ALL" to have the complete index. Use
   nword=0 or nword=NULL if the results on the words are not wanted } 
  \item{nseg}{index of the nfirst nseg repeated segments (by default 50). 
   Use nseg="ALL" to have the complete list of segments. Use 
   nseg=0 or nseg=NULL if the results on the segments are not wanted } 
  \item{ordFreq}{if ordFreq=TRUE, glossaries of words and repeated segments, are listed in frequency order; 
if ordFreq=FALSE, glossaries are listed  in alphabetic order (by default TRUE)}
  \item{file}{a connection, or a character string naming the file to print to in csv format. If NULL (the default), the results are not printed in a file}
  \item{sep}{character string to insert between the objects to print (if the argument file is not NULL) (by default ";")}
  \item{\dots}{further arguments passed to or from other methods,...}
}

\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica Bécue-Bertaut, Josep-Antón Sánchez-Espigares}

\seealso{ \code{\link{TextData}}, \code{\link{print.TextData}},  \code{\link{plot.TextData}}}


\examples{
# Non aggregate analysis
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), remov.number=TRUE, Fmin=10, Dmin=10,  
 stop.word.tm=TRUE, context.quali=c("Gender","Age_Group","Education"), context.quanti=c("Age"))
summary(res.TD)

# Aggregate analysis and repeated segments
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), var.agg="Gen_Age", remov.number=TRUE, 
 Fmin=10, Dmin=10, stop.word.tm=TRUE, context.quali=c("Gender","Age_Group","Education"), 
 context.quanti=c("Age"), segment=TRUE)
summary(res.TD)
}

\keyword{summary}
