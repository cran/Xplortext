\name{LexCA}
\encoding{latin1}
\alias{LexCA}

\title{Correspondence Analysis of a Lexical Table from a TextData object (LexCA)}

\description{
Performs Correspondence Analysis on the working lexical table contained in TextData object. Supplementary documents, words, segments, contextual quantitative and qualitative variables can be considered if previously selected in TextData function.
}

\usage{LexCA(object, ncp=5, context.sup="ALL", doc.sup=NULL, word.sup=NULL, 
  segment=FALSE, graph=TRUE, axes=c(1, 2), lmd=3, lmw=3)}

\arguments{
 \item{object}{object of TextData class}
 \item{ncp}{number of dimensions kept in the results (by default 5)}
 \item{context.sup}{column index(es) or name(s) of the contextual qualitative or quantitative variables among those selected in TextData function (by default "ALL")}
 \item{doc.sup}{vector indicating the index(es) or name(s) of the supplementary documents (rows) (by default NULL)}
 \item{word.sup}{vector indicating the index(es) or name(s) of the supplementary words (columns) (by default NULL)}
 \item{segment}{if TRUE, the repeated segments identified by TextData function will be considered as supplementary columns (by default FALSE)} 
 \item{graph}{if TRUE, basic graphs are displayed; use plot.LexCA to obtain more graphs (by default TRUE)}
 \item{axes}{length-2 vector indicating the axes to plot (by default axes=c(1,2))}
 \item{lmd}{only the documents whose contribution is over lmd times the average-document-contribution are plotted (by default lmd=3)}
 \item{lmw}{only the words whose contribution is over lmw times the average-word-contribution are plotted (by default lmw=3)}
}


\value{
Returns a list including:
\item{eig}{matrix with the eigenvalues, the percentages of inertia and the cumulative percentages of inertia}
\item{row}{list of matrices with all the results for the documents
(coordinates, square cosines, contributions, inertia)}
\item{col}{list of matrices with all the results for the words (coordinates, square cosines, contributions, inertia)}
\item{row.sup}{if row.sup is non-NULL, list of matrices with all the results for the supplementary documents (coordinates, square cosines)}
\item{col.sup}{if col.sup is non-NULL, list of matrices with all the results for the supplementary words (coordinates, square cosines)}
\item{quanti.sup}{if quanti.sup is non-NULL, list of matrices containing the results for the supplementary quantitative variables (coordinates, square cosines)}
\item{quali.sup}{if quali.sup is non-NULL, list of matrices with all the results for the supplementary categorical variables; see section details}
\item{meta}{list of the documents/words whose contribution is over lmd/lmw times the average document/word contribution}
\item{VCr}{Cramer's V coefficient}
\item{Inertia}{total inertia}
\item{info}{information about the corpus}
\item{segment}{if segment is TRUE, list of matrices with the results for the repeated segments (coordinates, square cosines)}
\item{var.agg}{name of the aggregation variable in the case of an aggregate correspondence analysis}
\item{call}{a list with some statistics}
}




\details{
In the case of a direct CA, DocTerm is a non-aggregate table and:
 \enumerate{
 \item the contextual quantitative variables are considered as supplementary quantitative columns in CA. 
 \item the categories of the contextual qualitative variables are considered as supplementary columns in CA.}
 
 In the case of an aggregate CA, DocTerm is an aggregate table and:
 \enumerate{
 \item the contextual quantitative variables are considered as supplementary quantitative columns in CA; the value of an active aggregate-document for a variable is the mean of the values corresponding to
 the source-documents belonging to this aggregate-document.
 \item the categories of the contextual qualitative variables are threatened as supplementary rows in CA; these rows contain the frequency with which each the set of documents belonging to this category has used the different words. }
}

\references{
Benzecri, J, P. (1981). Pratique de l'analyse des donnees. Linguistique & lexicologie (Vol.3). (P. Dunod., Ed).

Husson F., Le S., Pages J. (2011). Exploratory Multivariate Analysis by Example Using R. Chapman & Hall/CRC.

Lebart, L., Salem, A., & Berry, L. (1998). Exploring textual data. (D. Kluwer, Ed.).

Murtagh F. (2005). Correspondence Analysis and Data Coding with R and Java. Chapman & Hall/CRC.
}
\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica Bécue-Bertaut, Josep-Anton 
  Sánchez-Espigares\cr}

\seealso{ \code{\link{TextData}}, \code{\link{print.LexCA}},  \code{\link{plot.LexCA}}, \code{\link{summary.LexCA}}, \code{\link{ellipseLexCA}}}

\examples{
data(open.question)
\dontrun{
### non-aggregate CA
res.TD<-TextData(open.question, var.text=c(9,10), Fmin=10, Dmin=10,
        remov.number=TRUE, stop.word.tm=TRUE)
res.LexCA<-LexCA(res.TD, lmd=0, lmw=1)
}

### aggregate CA
res.TD<-TextData(open.question, var.text=c(9,10), var.agg="Age_Group", Fmin=10, Dmin=10,
        remov.number=TRUE, stop.word.tm=TRUE)
res.LexCA<-LexCA(res.TD, lmd=0, lmw=1)
}

\keyword{multivariate}




