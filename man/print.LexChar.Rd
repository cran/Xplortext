\name{print.LexChar}
\encoding{UTF-8}
\alias{print.LexChar}

\title{Print LexChar objects}

\description{
Prints characteristic words and documents from LexChar objects
}

\usage{
\method{print}{LexChar}(x, file = NULL, sep=";", ...) 
}

\arguments{
  \item{x}{object of LexChar class}
  \item{file}{a connection, or a character string giving the name of the file to print to (in csv format). If NULL (the default), the results are not printed in a file}
  \item{sep}{character to insert between the objects to print (if the argument file is non-NULL) (by default ";")}
  \item{\dots}{further arguments passed to or from other methods}
}

\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Mónica Bécue-Bertaut, Josep-Antón Sánchez-Espigares}

\seealso{ \code{\link{LexChar}}, \code{\link{plot.LexChar}}}


\examples{
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), var.agg="Gen_Edu", Fmin=10, Dmin=10,
        stop.word.tm=TRUE)
LD<-LexChar(res.TD, maxCharDoc = 0)
print(LD)
}

\keyword{print}
