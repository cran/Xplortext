\name{LabelTree}
\encoding{UTF-8} 
\alias{LabelTree}

\title{Hierarchical words (LabelTree)}

\description{Extracts the hierarchical characteristic words associated to the nodes of a hierarchical tree; the characteristic words of each node are extracted, then each word is associated to the node that it best characterizes.}

\usage{LabelTree(object, proba=0.05)}

\arguments{
 \item{object}{object of LexHCca or LexCHCca class}
 \item{proba}{threshold on the p-value when the characteristic words are computed (by default 0.05)}
}

\value{
Returns a list including:
 
\item{hierWord}{list of the characteristic words associated to the nodes of a hierarchical tree; 
only the non-empty nodes are included}

}

\references{
Bécue-Bertaut, M., Kostov, B., Morin, A., & Naro, G. (2014). Rhetorical Strategy in Forensic Speeches: 
Multidimensional Statistics-Based Methodology. Journal of Classification,31,85-106. \doi{10.1007/s00357-014-9148-9}.

Lebart, L., Salem, A., & Berry, L. (1998). Exploring textual data. (D. Kluwer, Ed.). \doi{10.1007/978-94-017-1525-6}.
}

\author{Monica Bécue-Bertaut, Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Josep-Anton 
Sánchez-Espigares, Belchin Kostov}

\seealso{ \code{\link{LexCA}}, \code{\link{LexCHCca}}}


\examples{
data(open.question)
res.TD<-TextData(open.question,var.text=c(9,10), var.agg="Age_Group", Fmin=10, Dmin=10,
        stop.word.tm=TRUE)
res.LexCA<-LexCA(res.TD, graph=FALSE)
res.LexCHCca<-LexCHCca(res.LexCA, nb.clust=4, min=3)
res.LabelTree<-LabelTree(res.LexCHCca)
}

\keyword{multivariate}

