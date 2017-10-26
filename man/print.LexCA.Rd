\name{print.LexCA}
\encoding{latin1}
\alias{print.LexCA}

\title{Print LexCA objects}

\description{
Prints the Textual Correspondence Analysis (CA) results from a LexCA object

}

\usage{
\method{print}{LexCA}(x, file = NULL, sep=";", ...) 
}

\arguments{
  \item{x}{object of LexCA class}
  \item{file}{a connection, or a character string giving the name od the file to print to (in csv format). If NULL (the default), the results are not printed in a file}
  \item{sep}{character to insert between the objects to print (if the argument file is non-NULL) (by default ";")}
  \item{\dots}{further arguments passed from other methods}
}

\author{Ram�n Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica B�cue-Bertaut, Josep-Ant�n S�nchez-Espigares}

\seealso{ \code{\link{LexCA}}, \code{\link{plot.LexCA}},  \code{\link{summary.LexCA}}}

\examples{
data(open.question)
res.TD<-TextData(open.question,var.text=c(9,10), var.agg="Age_Group", Fmin=10, Dmin=10,
        remov.number=TRUE, stop.word.tm=TRUE)
res.LexCA<-LexCA(res.TD,lmd=0,lmw=1)
print(res.LexCA)
}

\keyword{print}