\name{print.TextData}
\encoding{latin1}
\alias{print.TextData}

\title{Print TextData objects}

\description{
Print statistical results for documents, words and segments from TextData objects, in alphabetical and frequency order.
}

\usage{
\method{print}{TextData}(x, file = NULL, sep=";", ...) 
}

\arguments{
  \item{x}{object of TextData class}
  \item{file}{connection, or character string giving the name of the file to print to (in csv format). If NULL (by default value), 
the results are not printed in a file}
  \item{sep}{character inserted between the objects to print (if file argument is non-NULL) (by default ";")}
  \item{\dots}{further arguments passed to or from other methods}
}

\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica Bécue-Bertaut, Josep-Antón Sánchez-Espigares}

\seealso{ \code{\link{TextData}}, \code{\link{plot.TextData}},  \code{\link{summary.TextData}}}


\examples{
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), remov.number=TRUE, Fmin=10, Dmin=10,  
stop.word.tm=TRUE, context.quali=c("Gender","Age_Group","Education"),
   context.quanti=c("Age"))
print(res.TD)
}

\keyword{print}
