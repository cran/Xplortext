\name{print.LexCA}
\encoding{UTF-8}
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
  \item{file}{a connection, or a character string giving the name of the file to print to (in csv format). If NULL (the default), the results are not printed in a file}
  \item{sep}{character to insert between the objects to print (if the argument file is non-NULL) (by default ";")}
  \item{\dots}{further arguments passed to or from other methods}
}


\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Mónica Bécue-Bertaut, Josep-Antón Sánchez-Espigares}

\seealso{ \code{\link{LexCA}}, \code{\link{plot.LexCA}}, \code{\link{summary.LexCA}}, \code{\link{TextData}}}


\examples{
data(open.question)
res.TD<-TextData(open.question,var.text=c(9,10), var.agg="Age_Group", Fmin=10, Dmin=10,
        remov.number=TRUE, stop.word.tm=TRUE)
res.LexCA<-LexCA(res.TD,lmd=0,lmw=1)
print(res.LexCA)
}

\keyword{print}

