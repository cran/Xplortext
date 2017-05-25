\name{LexHCca}
\encoding{latin1}
\alias{LexHCca}

\title{Hierarchical Clustering of Documents on Textual Correspondence Analysis Coordinates (LexHCca)}

\description{Agglomerative hierarchical clustering on a corpus of documents.}

\usage{LexHCca(object, nb.clust=0, consol=TRUE, iter.max=10, min=3, max=NULL, 
 order=TRUE, nb.par=5, edit.par=FALSE, graph=TRUE, proba=0.05,...)}
	
\arguments{
\item{object}{object of LexCA class}
  \item{nb.clust}{number of clusters (see details). If 0, the tree is cut at the level the user clicks on. If -1, the tree is automatically 
cut at the suggested level. If a (positive) integer, the tree is cut with nb.clust clusters (by default 0)}
\item{consol}{if TRUE, k-means consolidation is performed (by default TRUE)}
\item{iter.max}{maximum number of iterations in the consolidation step (by default 10)}
\item{min}{minimum number of clusters (by default 3)}
\item{max}{maximum number of clusters (by default NULL and then max is computed as the minimum between 10 and the number of documents divided by 2)}
\item{order}{if TRUE, the clusters are numbered depending on the coordinate of their centroid on the first axis (by default TRUE)}
\item{nb.par}{number of edited paragons (para) and specific documents (dist) (by default 5)}
\item{edit.par}{if TRUE, the literal text of the parangon and specific documents are listed in the results (by default FALSE)}
\item{graph}{if TRUE, graphs are displayed (by default TRUE)}
\item{proba}{threshold on the p-value used in selecting words, documents, axes and contextual variables when describing the clusters (by default 0.05)}
\item{...}{other arguments from other methods}
}

\value{
 Returns a list including:
 \item{data.clust}{the original active lexical table used in LexCA plus a new column called clust containing the partition}
 \item{desc.wordvar}{description of the clusters by their characteristic words and, if contextual variables were considered in LexCA, 
       description of the partition/clusters by these variables}
 \item{desc.axes}{description of the clusters by the characteristic axes}
 \item{call}{list of internal objects. \code{call$t} giving the results for the hierarchical tree;\cr See the first reference for more details}
 \item{desc.doc}{labels of the paragon (para) and specific documents (dist) of each cluster}
 \item{clust.count}{count of documents belonging to each cluster}
 \item{clust.content}{list of the document labels according to the cluster they belong to}
 \item{docspara}{if edit.par=TRUE, description of the clusters by the literal text of the nb.par "para" documents}
 \item{docsdist}{if edit.par=TRUE, description of the clusters by the literal text of the nb.par "dist" documents}
Returns the hierarchical tree with a barplot of the successive inertia gains, the 
CA map of the documents enriched by the tree (3D), the CA map with the document labels colored according to their cluster (2D).
}


\details{
LexHCca starts from the documents coordinates on textual correspondence analysis axes. 
Euclidean metric and Ward method are used. 

The number of clusters is determined either a priori or from the hierarchical tree structure. 
If nb.clust=0, a level for cutting the tree is automatically suggested. 
This is computed in the following way, reading the tree downward. 
At a given step, the tree could be cut into Q clusters (Q varying between min and max). 
The between-inertia gain when passing from Q-1 to Q clusters and the between-inertia gain when passing from Q to Q+1 clusters are computed. 
The suggested level corresponds to the maximum value of the ratio between the former and the latter of these inertia-gains.
Note that the between-inertia gain when passing from Q to Q+1 clusters is equal to the value of the Ward criterion when passing from Q+1 to Q clusters when building the tree bottom up.
In this latter case, a level where to cut the tree is suggested. The user can choose to cut the tree at this level or at another one.
 
The results include a thorough description of the clusters, taking into account contextual variables. Graphs are provided.

}

\references{ 
Husson F., Le S., Pages J. (2011). Exploratory Multivariate Analysis by Example Using R. Chapman & Hall/CRC.

Lebart, L., Salem, A., & Berry, L. (1998). Exploring textual data. (D. Kluwer, Ed.).
}

\author{Monica Bécue-Bertaut, Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Josep-Anton Sánchez-Espigares}

\seealso{ \code{\link{LexCA}}}

\examples{
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), Fmin=10, Dmin=10, stop.word.tm=TRUE,
        context.quali=c("Gender","Age_Group","Education"), context.quanti=c("Age"))
res.LexCA<-LexCA(res.TD, graph=FALSE, ncp=8)
res.hcca<-LexHCca(res.LexCA, graph=TRUE, nb.clust=5, order=TRUE)
}


\keyword{multivariate}


