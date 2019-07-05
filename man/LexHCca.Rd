\name{LexHCca}
\encoding{UTF-8}
\alias{LexHCca}

\title{Hierarchical Clustering of Documents on Textual Correspondence Analysis Coordinates (LexHCca)}

\description{Agglomerative hierarchical clustering on a corpus of documents.}

\usage{LexHCca(object, nb.clust=0, consol=FALSE, iter.max=10, min=3, max=NULL, 
   kk=Inf, order=TRUE, graph=TRUE, proba=0.05,cluster.CA="rows",description=TRUE,
   nb.par=0,size.par=80,marg.doc=FALSE,seed=12345,...)}


\arguments{
\item{object}{object of LexCA class}
\item{nb.clust}{number of clusters (see details). If 0, the tree is cut at the level the user clicks on. If -1, the tree is automatically 
cut at the suggested level. If a (positive) integer, the tree is cut with nb.clust clusters (by default 0)}
\item{consol}{if TRUE, consolidation is performed after hierachical clustering (by default FALSE)}
\item{iter.max}{maximum number of iterations in the consolidation step (by default 10)}
\item{min}{minimum number of clusters (by default 3)}
\item{max}{maximum number of clusters (by default NULL and then max is computed as the minimum between 10 and the number of documents divided by 2)}
\item{kk}{An integer corresponding to the number of clusters used in a Kmeans preprocessing before the hierarchical clustering; the top of the hierarchical tree is then constructed from this partition. This is very useful if the number of individuals is high. Note that consolidation cannot be performed if kk is different from Inf and some graphics are not drawn. Inf is used by default and no preprocessing is done, all the graphical outputs are then given.}
\item{order}{if TRUE, the clusters are numbered depending on the coordinate of their centroid on the first axis (by default TRUE)}
\item{graph}{if TRUE, graphs are displayed (by default TRUE)}
\item{proba}{threshold on the p-value used in selecting words, documents, axes and contextual variables when describing the clusters (by default 0.05)}
\item{cluster.CA}{if 'rows' or 'docs' cluster is performed with documents; 'columns' or 'words' with words (by default 'rows')}
\item{description}{if TRUE, description of the clusters by their characteristic words/documents, by the characteristic axes and 
by contextual variables if considered in LexCA (by default TRUE)}
\item{nb.par}{number of edited paragons (para) and specific documents (dist) (by default 0)}
\item{size.par}{text size of edited paragons (para) and specific documents (dist) (by default 80)}
\item{marg.doc}{if FALSE, frequencies before TextData selection are the marginal frequencies for documents in description analysis, 
 TRUE if frequencies after TextData selection (by default FALSE)}
\item{seed}{Seed to obtain the same results using k-means (by default 12345)}
\item{...}{other arguments from other methods}
}

\value{
 Returns a list including:
 \item{data.clust}{the original active lexical table used in LexCA plus a new column called clust containing the partition}
 \item{centers}{coordinates of centers from LexCA results for each cluster} 
 \item{clust.count}{count of documents/words belonging to each cluster and some statistics}
 \item{clust.content}{list of the document/word labels according to the cluster they belong to}
 \item{ss}{total sum of squares}
 \item{call}{list of internal objects. \code{call$t} giving the results for the hierarchical tree;\cr See the first reference for more details}
 \item{desc.axes}{description of the clusters by the characteristic axes}
 \item{desc.wordvar}{if description=TRUE, description of the clusters by their characteristic words, supplementary words and, if contextual variables were considered in LexCA, 
       description of the partition/clusters by these variables}
 \item{desc.doc}{if description=TRUE, description of the clusters by their characteristic documents}
 \item{wordslabels}{labels of the paragon (para) and specific words (dist) of each cluster}
 \item{docslabels}{labels of the paragon (para) and specific documents (dist) of each cluster}
 \item{docspara}{if nb.par>0, description of the clusters by the nb.par "para" documents writing the first size.par characters of the
 literal text}
 \item{docsdist}{if nb.par>0, description of the clusters by the nb.par "dist" documents writing the first size.par characters of the
 literal text}

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
Husson F., Lê S., Pagès J. (2011). Exploratory Multivariate Analysis by Example Using R. Chapman & Hall/CRC. \doi{10.1201/b10345}.

Lebart, L., Salem, A., & Berry, L. (1998). Exploring textual data. (D. Kluwer, Ed.). \doi{10.1007/978-94-017-1525-6}.
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


