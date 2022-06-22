\name{LexCHCca}
\encoding{UTF-8} 
\alias{LexCHCca}

\title{Chronological Constrained Hierarchical Clustering on Correspondence Analysis Components (LexCHCca)}
\description{Chronological constrained agglomerative hierarchical clustering on a corpus of documents}

\usage{
LexCHCca (object, nb.clust=0, min=2, max=NULL, nb.par=5, 
 graph=TRUE, proba=0.05, cut.test=FALSE, alpha.test =0.05, description=FALSE,
 nb.desc=5, size.desc=80)
}


\arguments{
  \item{object}{object of LexCA class}
  \item{nb.clust}{number of clusters only if no test (cut.test=FALSE). If 0 (or "click"), the tree is cut at the level the user clicks on. If -1 (or "auto"), the tree is automatically cut at the suggested level. If a (positive) integer, the tree is cut with nb.clust clusters (by default 0)}
  \item{min}{minimum number of clusters. Available only if cut.test=FALSE. (by default 3)}
  \item{max}{maximum number of clusters. Available only if cut.test=FALSE. (by default NULL; then max is computed as the minimum between 10 and the number of documents divided by 2)}
  \item{nb.par}{number of edited paragons (para) and specific documents labels (dist) (by default 5)}
  \item{graph}{if TRUE, graphs are displayed (by default TRUE)}
  \item{proba}{threshold on the p-value used to describe the clusters (by default 0.05)}
  \item{cut.test}{if FALSE (by default), Legendre test is not performed when joining two nodes. This test is used to determine whether two clusters should be joined or not; see details}
  \item{alpha.test}{threshold on the p-value used in selecting aggregation clusters for Legendre test (by default 0.05)}
  \item{description}{if TRUE, description of the clusters by the characteristic words/documents, paragon (para), specific documents (dist) and contextual variables if these latter have been selected in the previous LexCA function (by default FALSE)}
  \item{nb.desc}{number of paragons (para) and specific documents (dist) that are edited when describing the clusters (by default 5)}
  \item{size.desc}{maximum of characters when editing the paragons (para) and specific documents (dist) to describe the clusters (by default 80)}
}

\value{
Returns a list including:
 \item{data.clust}{the active lexical table used in LexCA plus a new column called Clust_ containing the partition}
 \item{coord.clust}{coordinates table issued from CA plus a new column called weigths and another column called Clust_, corresponds to the partition}
 \item{centers}{coordinates of the gravity centers of the clusters}
 \item{description}{
 $des.word for description of the clusters of documents by their characteristic words, the paragons (des.doc$para) and 
  specific documents (des.doc$dist) of each cluster; see details}
 \item{call}{list of internal objects. \code{call$t} giving the results for the hierarchical tree}
 \item{dendro}{hclust object. This allows for using the dendrogram in other packages}
\item{phases}{details of the tracking of the agglomerative hierarchical process. In particular, the cut points (joining documents not allowed) can be identified}
\item{sum.squares}{sum of squares decomposition for documents and clusters}
}

\details{
LexCHCca starts from the document coordinates issued from a textual correspondence analysis. 
The hierarchical tree is built in such a way that only chronological contiguous nodes can be joined. 
The documents have to be ranked in their chronological order in the source-base (data frame format) before to apply the function (TextData format).

Legendre test allows to determine whether the fusion between two nodes based on their contiguity lead to a heterogenous new node (no homogeneity-between-clusters).
If Legendre test is applied (cut.test=TRUE), the number of clusters is the number obtained by the test and nb.clust has not effects.

If no Legendre test is applied (cut.test= FALSE), the number of clusters is determined either a priori or from the constrained hierarchical tree structure.  

The object $para contains the distance between each document and the centroid of its class.
 
The object $dist contains the distance between each document and the centroid of the farthest cluster.

The results of the description of the clusters and graphs are provided.

}


\references{
Bécue-Bertaut, M., Kostov, B., Morin, A., & Naro, G. (2014). Rhetorical Strategy in Forensic Speeches: Multidimensional Statistics-Based Methodology. Journal of Classification,31, 85-106. \doi{10.1007/s00357-014-9148-9}.

Husson F., Lê S., Pagès J. (2017). Exploratory Multivariate Analysis by Example Using R. Chapman & Hall/CRC. \doi{10.1201/b21874}.

Lebart L. (1978). Programme d'agrégation avec contraintes. Les Cahiers de l'Analyse des Données, 3, pp. 275--288.

Legendre, P. & Legendre, L. (1998), Numerical Ecology (2nd ed.), Amsterdam: Elsevier Science.

Murtagh F. (1985). Multidimensional Clustering Algorithms. Vienna: Physica-Verlag, COMPSTAT Lectures.

}

\author{Monica Bécue-Bertaut, Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Josep-Antón Sánchez-Espigares, 
Belchin Kostov\cr}

\seealso{\code{\link{plot.LexCHCca}}, \code{\link{LexCA}}}

\examples{
data(open.question)
res.TD<-TextData(open.question,var.text=c(9,10), var.agg="Age_Group", Fmin=10, Dmin=10, 
        stop.word.tm=TRUE)
res.LexCA<-LexCA(res.TD, graph=FALSE)
res.ccah<-LexCHCca(res.LexCA, nb.clust=4, min=3)
}
\keyword{multivariate}

