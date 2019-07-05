\name{LexCHCca}
\encoding{UTF-8} 
\alias{LexCHCca}

\title{Chronogically Constrained Agglomerative Hierarchical Clustering on Correspondence Analysis Components (LexCHCca)}
\description{Chronogically constrained agglomerative hierarchical clustering on a corpus of documents.}

\usage{
LexCHCca (object, nb.clust=0, min=3, max=NULL, nb.par=5, graph=TRUE, proba=0.05) 
}

\arguments{
  \item{object}{object of LexCA class}
  \item{nb.clust}{number of clusters (see details). If 0, the tree is cut at the level the user clicks on. If -1, the tree is automatically 
cut at the suggested level. If a (positive) integer, the tree is cut with nb.clust clusters (by default 0)}
  \item{min}{minimum number of clusters (by default 3)}
  \item{max}{maximum number of clusters (by default NULL and then max is computed as the minimum between 10 and the number of documents divided by 2)}
  \item{nb.par}{number of edited paragons (para) and specific documents labels (dist) (by default 5)}
  \item{graph}{if TRUE, graphs are displayed (by default TRUE)}
  \item{proba}{threshold on the p-value used in selecting the characteristic words of the clusters and in selecting the axes when describing the clusters by the axes (by default 0.05)}
}

\value{
Returns a list including:
  \item{data.clust}{the original active lexical table with a supplementary column called clust containing the partition}
  \item{desc.word}{description of the clusters by their characteristic words}
  \item{desc.axes}{description of the clusters by the characteristic axes}
  \item{call}{list or parameters and internal objects}
  \item{desc.doc}{labels of the paragon (para) and specific documents (dist) of each cluster}
  \item{dendro}{list with the succession of nodes that are found when reading the tree downward}
Returns the graphs with the tree and the correspondence analysis map where the documents are colored according to the cluster they belong to (2D).
}

\details{
LexCHCca starts from the documents coordinates on textual correspondence analysis axes. 
The hierarchical tree is built taking into account that only chronological contiguous nodes can be grouped. 
The documents have to be ranked in the lexical table in the chronological order. 
Euclidean metric and complete linkage method are used.

The number of clusters is determined either a priori or from the constrained hierarchical tree structure. 
If nb.clust=0, a level for cutting the tree is automatically suggested. 
This is computed in the following way, reading the tree downward. 
At a given step, the tree could be cut into Q clusters (Q varying between min and max). 
The distance between the two nodes that are no longer grouped together using complete linkage method when passing from Q-1 to Q clusters and the distance between the two nodes that are no longer grouped together when passing from Q to Q+1 are computed. The suggested level corresponds to the maximum value of the ratio between the former and the latter of these values. These distances correspond to the criterion value when building the tree bottom up.
The user can choose to cut the tree at this level or at another one.

The results include a thorough description of the clusters. Graphs are provided.

The tree is plotted jointly with a barchart of the successive values of the aggregation criterion.
}


\references{
Bécue-Bertaut, M., Kostov, B., Morin, A., & Naro, G. (2014). Rhetorical Strategy in Forensic Speeches: Multidimensional Statistics-Based Methodology. Journal of Classification,31, 85-106. \doi{10.1007/s00357-014-9148-9}.

Lebart L. (1978). Programme d'agrégation avec contraintes. Les Cahiers de l'Analyse des Données, 3, pp. 275--288.

Legendre, P. & Legendre, L. (1998), Numerical Ecology (2nd ed.), Amsterdam: Elsevier Science.

Murtagh F. (1985). Multidimensional Clustering Algorithms. Vienna: Physica-Verlag, COMPSTAT Lectures.

}

\author{Monica Bécue-Bertaut, Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Josep-Antón Sánchez-Espigares, 
Belchin Kostov\cr}

\seealso{\code{\link{plot.LexCHCca}}, \code{\link{LabelTree}}, \code{\link{LexCA}}}

\examples{
data(open.question)
res.TD<-TextData(open.question,var.text=c(9,10), var.agg="Age_Group", Fmin=10, Dmin=10, 
        stop.word.tm=TRUE)
res.LexCA<-LexCA(res.TD, graph=FALSE)
res.ccah<-LexCHCca(res.LexCA, nb.clust=4, min=3)
}
\keyword{multivariate}

