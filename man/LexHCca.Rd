\name{LexHCca}
\encoding{UTF-8}
\alias{LexHCca}

\title{Hierarchical Clustering on Textual Correspondence Analysis Coordinates (LexHCca)}

\description{Agglomerative hierarchical clustering of documents or words issued from correspondence analysis coordinates}


\usage{LexHCca(x, cluster.CA="docs", nb.clust="click", min=2, max=NULL, kk=Inf, 
   consol=FALSE, iter.max=500, graph=TRUE, description=TRUE, 
   proba=0.05, nb.desc=5, size.desc=80, seed=12345,...)}
   

\arguments{
\item{x}{object of LexCA class}
\item{cluster.CA}{if "rows" or "docs" cluster analysis is performed on documents; if "columns" or "words", cluster analysis is performed on words (by default "docs")}
\item{nb.clust}{number of clusters. If 0 (or "click"), the tree is cut at the level the user clicks on.
If -1 (or "auto"), the tree is automatically cut at the suggested level. 
If a (positive) integer, the tree is cut with nb.clust clusters (by default "click")}
\item{min}{minimum number of clusters (by default 2)}
\item{max}{maximum number of clusters (by default NULL, then max is computed as the minimum between 10 and the number of documents divided by 2)}
\item{kk}{in case the user wants to perform a Kmeans clustering previously to the hierarchical clustering (preprocessing step), kk is an integer corresponding to the number of clusters of this previous partition. Further, the hierarchical tree is constructed starting from the nodes of this partition as terminal elements. This is very useful when the number of elements to be classified is very large. By default, the value is Inf and no Kmeans preprocessing is performed}
\item{consol}{if TRUE, a Kmeans consolidation step is performed after the hierarchical clustering (consolidation cannot be performed if kk is used and equals a number) (by default FALSE)}
\item{iter.max}{maximum number of iterations in the consolidation step (by default 500)}
\item{graph}{if TRUE, graphs are displayed (by default TRUE)}
\item{description}{if TRUE, description of the clusters of documents or words by the axes, the characteristic words in the case of clustering documents or the characteristic documents in the case of clustering words. The documents or words considered as paragon (para) or specific (dist) are identified. In the case of clustering documents, contextual variables also characterize the clusters. These variables have to be selected in LexCA (by default TRUE)}
\item{proba}{threshold on the p-value used in selecting the elements characterizing significantly the clusters (by default 0.05)}
\item{nb.desc}{Maximum of characters when editing the paragons (para) and specific documents (dist) to describe the clusters (by default 80))}
\item{size.desc}{text size of edited paragons (para) and specific documents (dist) when describing the clusters of documents (by default 80)}
\item{seed}{Seed to obtain the same results in successive Kmeans (by default 12345)}
\item{...}{other arguments from other methods}
}

\value{
 Returns a list including:
 \item{data.clust}{the active lexical table used in LexCA plus a new column called Clust_ containing the partition}
 \item{coord.clust}{coordinates table issued from CA plus a new column called Clust_ containing the partition}
 \item{centers}{coordinates of the gravity centers of the clusters} 
  \item{clust.count}{counts of documents/words belonging to each cluster and contribution of the clusters to the variability decomposition}
 \item{clust.content}{list of the document/word labels according to the cluster they belong to}
 \item{call}{list of internal objects. \code{call$t} giving the results for the hierarchical tree. See the second reference for more details}
 \item{description}{$desc.axes for description of the clusters by the characteristic axes ($axes) and eta-squared between axes and clusters ($quanti.var).

 $des.cluster.doc for description of the clusters by their characteristic words ($word), supplementary words ($wordsup) and, 
 if contextual variables were considered in LexCA, description of the partition/clusters by qualitative ($qualisup) and quantitative ($quantisup) variables, paragons ($para) and specific words ($dist) of each cluster.

  $des.word.doc description of the clusters of words by their characteristic documents ($docs), paragons ($para) and 
  specific documents ($dist) of each cluster. 

}


Returns the hierarchical tree with a barplot of the successive inertia gains, and the 
first CA map of the documents/words. The labels are colored according to the cluster.
}


\details{
LexHCca starts from the documents/words coordinates issued from correspondence analysis axes. 
Euclidean metric and Ward method are used. 

If the agglomerative clustering starts from many elements (documents or words), it is possible to previously perform a Kmeans partition with kk clusters to further build the tree from these (weighted) kk clusters. 
 
The object $para contains the distance between each document and the centroid of its class. 

The object $dist contains the distance between each document and the centroid of the farthest cluster.

The results include a thorough description of the clusters. Graphs are provided.

}



\references{ 

Bécue-Bertaut M. Textual Data Science with R. Chapman & Hall/CRC.  \doi{10.1201/9781315212661}.

Husson F., Lê S., Pagès J. (2017). Exploratory Multivariate Analysis by Example Using R. Chapman & Hall/CRC. \doi{10.1201/b21874}.

Lebart, L., Salem, A., & Berry, L. (1998). Exploring textual data. (D. Kluwer, Ed.). \doi{10.1007/978-94-017-1525-6}.
}

\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica Bécue-Bertaut, Josep-Anton Sánchez-Espigares}

\seealso{ \code{\link{LexCA}}, \code{\link{plot.LexHCca}}}

\examples{
data(open.question)	
res.TD<-TextData(open.question, var.text=c(9,10), Fmin=10, Dmin=10, stop.word.tm=TRUE,	
        context.quali=c("Gender","Age_Group","Education"), context.quanti=c("Age"))	
res.LexCA<-LexCA(res.TD, graph=FALSE, ncp=8)	
res.hcca<-LexHCca(res.LexCA, graph=FALSE, nb.clust=5)	
}


\keyword{multivariate}


