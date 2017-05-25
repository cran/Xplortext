\name{LexChar}
\encoding{latin1}
\alias{LexChar}


\title{Characteristic words and documents (LexChar)}

\description{Characteristic words of documents from TextData objects.}

\usage{LexChar(object, proba=0.05, maxDocs=20, maxCharDoc=10, maxPrnDoc=100)}

\arguments{
 \item{object}{TextData object}
 \item{proba}{threshold on the p-value used when selecting the characteristic words (by default 0.05)}
 \item{maxDocs}{maximum number of documents in the working lexical table (by default 20). See details}
 \item{maxCharDoc}{maximum number of characteristic source-documents to extract (by default 10). See details}
 \item{maxPrnDoc}{maximum length to be printed for a characteristic document (by default 100 characters)}
}

\value{
Returns a list including:
 
\item{CharWord}{characteristic words of all the documents}
\item{CharDoc}{characteristic source-documents of all the aggregate-documents}
}

\details{
The lexical table provided by TextData can consider either source-documents or aggregate-documents, in accordance with the value of argument "var.agg" in TextData. 
Extracting the characteristic words for a too high number of documents is of no interest and time-consuming. So that, this function can be applied only when the number of documents in the lexical table is under or equal to maxDocs (by default 20). 
In the case of aggregate documents, extracting the characteristic source-documents is possible but of interest only if the source-documents are not too long. 
In any case, only the first maxPrnDoc characters of each characteristic document are printed (by default 100).
}

\references{
 
Lebart, L., Salem, A., & Berry, L. (1998). Exploring textual data. (D. Kluwer, Ed.).

}

\author{Monica Bécue-Bertaut, Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Josep-Antón
Sánchez-Espigares\cr}

\seealso{ \code{\link{TextData}}, \code{\link{print.LexChar}}, \code{\link{plot.LexChar}}}

\examples{
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), var.agg="Gen_Edu", Fmin=10, Dmin=10, 
        remov.number=TRUE, stop.word.tm=TRUE)
LexChar(res.TD)
}

\keyword{multivariate}

