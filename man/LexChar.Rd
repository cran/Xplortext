\name{LexChar}
\encoding{UTF-8}
\alias{LexChar}

\title{Characteristic words and documents (LexChar)}

\description{Measure of the association between vocabulary or words and quantitative or qualitative contextual variables.}

\usage{LexChar(object, proba=0.05, maxCharDoc=10, maxPrnDoc=100, 
              marg.doc="before",  context=NULL, correct=TRUE, nbsample=500,
              seed=12345,...)}

\arguments{
 \item{object}{TextData, DocumentTermMatrix, dataframe or matrix object}
 \item{proba}{threshold on the p-value used when selecting the characteristic words (by default 0.05)}
 \item{maxCharDoc}{maximum number of characteristic source-documents to extract (by default 10). See details}
 \item{maxPrnDoc}{maximum length to be printed for a characteristic document (by default 100 characters)}
 \item{marg.doc}{if after/before, frequencies after/before TextData selection are used as document weighting (by default "before"); 
  if before.RW all words under threshold in TextData function are included as a new word named RemovedWords}
 \item{context}{name of quantitative or qualitative variables}
 \item{correct}{if TRUE, pvalue correction test is applied for quantitative contextual variables (by default TRUE)}
 \item{nbsample}{number of samples drawn to evaluate the pvalues in quantitative contextual variables}
 \item{seed}{Seed to obtain the same results using permutation tests (by default 12345)}
 \item{\dots}{further arguments passed to or from other methods}
}

\value{
Returns a list including:
 
\item{CharWord}{characteristic words of all the documents}
\item{stats}{association statistics of the lexical table}
\item{CharDoc}{characteristic source-documents of all the aggregate-documents including qualitative contextual variables}
\item{Vocab}{characteristic quantitative and qualitative variables of the words. CharWord and stats are provided.}
}

\details{
The lexical table provided by TextData can consider either source-documents or aggregate-documents, in accordance with the value of argument "var.agg" in TextData.
Context cualitative variables allow to aggregate documents by combining the categories of the qualitative variables and the aggregation variable if any.

Extracting the characteristic words (CharWord) for a too high number of documents is of no interest and time-consuming. 

In any case, only the first maxPrnDoc characters of each characteristic document are printed (by default 100).

In the case of the association between words and qualitative variables, the usual characteristic words are provided.

quali$CharWord provides the qualitative variables (including the aggregation variable) and their categories.
quali$stats provides association statistics for vocabulary and qualitative variables (including the aggregation variable).
quali$CharDoc provides characteristic source-documents for the categories.
quanti$CharWord provides characteristic quantitative variables for each word. If there are aggregation variable and/or qualitative contextual variable, from aggregated lexical table.
quanti$stats provides statistics for vocabulary and quantitative variables. If there are aggregation variable and/or qualitative contextual variable, from aggregated lexical table.

If the lexical table (object) is not a TextData object, context argument can be columns of the same dataframe.
The aggregate lexical table is constructed from the combinations of the categories of the qualitative variables (including the aggregation variable).
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
 res.LexChar <-LexChar(res.TD)
 summary(res.LexChar)
}

\keyword{multivariate}

