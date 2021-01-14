\name{LexChar}
\encoding{UTF-8}
\alias{LexChar}

\title{Characteristic words and documents (LexChar)}

\description{Characteristic words of documents from TextData or DocumentTermMatrix objects.}

\usage{LexChar(object, proba=0.05, maxDocs=20, maxCharDoc=10, maxPrnDoc=100, 
              marg.doc="before", correct=TRUE,  context.sup="ALL", nbsample=500,
              seed=12345)}

\arguments{
 \item{object}{TextData, DocumentTermMatrix, dataframe or matrix object}
 \item{proba}{threshold on the p-value used when selecting the characteristic words (by default 0.05)}
 \item{maxDocs}{maximum number of documents in the working lexical table (by default 20). See details}
 \item{maxCharDoc}{maximum number of characteristic source-documents to extract (by default 10). See details}
 \item{maxPrnDoc}{maximum length to be printed for a characteristic document (by default 100 characters)}
 \item{marg.doc}{if before, frequencies before TextData selection are used as document weighting; if after, frequencies after TextData selection (by default before); if object is not TextData only after selection is used}
 \item{correct}{if TRUE, pvalue correction test is applied for WordsQuanti analysis (by default TRUE)}
 \item{context.sup}{column index(es) or name(s) of the contextual quantitative variables among those selected in TextData function (by default "ALL"); if object is not TextData context.sup must be a data frame with the same rows (documents or categories in the case of aggregated documents) as object}
 \item{nbsample}{number of samples drawn to evaluate the pvalues in quantitative contextual variables}
 \item{seed}{Seed to obtain the same results using permutation tests (by default 12345)}
}

\value{
Returns a list including:
 
\item{CharWord}{characteristic words of all the documents}
\item{stats}{association statistics of the lexical table}
\item{CharDoc}{characteristic source-documents of all the aggregate-documents}
\item{Vocab}{characteristic quantitative and qualitative variables of the words. CharWord and stats are provided}
}

\details{
The lexical table provided by TextData can consider either source-documents or aggregate-documents, in accordance with the value of argument "var.agg" in TextData.

If the object is a DocumentTermMatrix, dataframe or matrix object it is possible to build aggregate-documents, in accordance with the value of argument "context.sup" in TextData.

Extracting the characteristic words (CharWord) for a too high number of documents is of no interest and time-consuming. So that, this function can be applied only when the number of documents in the lexical table is under or equal to maxDocs (by default 20). 

In the case of aggregate documents, extracting the characteristic source-documents is possible but of interest only if the source-documents are not too long. 
In any case, only the first maxPrnDoc characters of each characteristic document are printed (by default 100).

Vocab$quali$CharWord provides the qualitative variables and their categories.
Vocab$quali$stats provides association statistics for vocabulary and qualitative variables.
Vocab$quanti$CharWord provides characteristic quantitative variables for each word. 
Vocab$quali$stats provides statistics for vocabulary and quantitative variables.
}

\references{
 
Lebart, L., Salem, A., & Berry, L. (1998). Exploring textual data. (D. Kluwer, Ed.). \doi{10.1007/978-94-017-1525-6}.
}

\author{Monica Bécue-Bertaut, Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Josep-Antón
Sánchez-Espigares, Belchin Kostov}

\seealso{ \code{\link{TextData}}, \code{\link{print.LexChar}}, \code{\link{plot.LexChar}}, \code{\link{summary.LexChar}}}

\examples{
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), var.agg="Gen_Edu", Fmin=10, Dmin=10, 
        remov.number=TRUE, stop.word.tm=TRUE)
# res.LexChar <-LexChar(res.TD)

}

\keyword{multivariate}

