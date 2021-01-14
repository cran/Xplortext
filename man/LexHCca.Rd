\name{LexHCca}
\encoding{UTF-8}
\alias{LexHCca}

\title{Hierarchical Clustering of Documents on Textual Correspondence Analysis Coordinates (LexHCca)}

\description{Agglomerative hierarchical clustering on documents or words}

\usage{LexHCca(x, cluster.CA="docs", marg.doc="before", nb.clust="click", min=2, 
   max=NULL, kk=Inf, consol=FALSE, iter.max=500, graph=TRUE, description=TRUE, 
   proba=0.05, nb.desc=5, size.desc=80, seed=12345,...)}
   

\arguments{
\item{x}{object of LexCA class}
\item{cluster.CA}{if "rows" or "docs" cluster analysis is performed on documents; "columns" or "words" and cluster is performed on words (by default "docs")}
\item{marg.doc}{This parameter concerns only the description of the clusters by the characteristic words in case of using TextData function to build the lexical table. To compute the corresponding p-value, document weightings correspond to their relative lengths either "before" or "after" word selection in TextData (by default before)}
\item{nb.clust}{number of clusters (see details). If 0 (or "click"), the tree is cut at the level the user clicks on.
If -1 (or "auto"), the tree is automatically cut at the suggested level. 
If a (positive) integer, the tree is cut with nb.clust clusters (by default "click")}
\item{min}{minimum number of clusters (by default 2)}
\item{max}{maximum number of clusters (by default NULL and then max is computed as the minimum between 10 and the number of documents divided by 2)}
\item{kk}{An integer corresponding to the number of clusters used in a Kmeans preprocessing before the hierarchical clustering; the top of the hierarchical tree is then constructed from this partition. This is very useful if the number of individuals is high. Inf is used by default and no Kmeans preprocessing is perdformed. It can be used with nb.clust values}
\item{consol}{if TRUE, consolidation is performed after hierachical clustering (by default FALSE)}
\item{iter.max}{maximum number of iterations in the consolidation step (by default 500)}
\item{graph}{if TRUE, graphs are displayed (by default TRUE)}
\item{description}{if TRUE, description of the clusters by characteristic words/documents,
characteristic axes, paragon (para), specific words/documents (dist) and contextual variables if considered in LexCA (by default TRUE)}
\item{proba}{threshold on the p-value used in selecting words, documents, axes and contextual variables when describing the clusters (by default 0.05)}
\item{nb.desc}{number of edited paragons (para) and specific documents/words (dist) when describing the clusters (by default 5)}
\item{size.desc}{text size of edited paragons (para) and specific documents (dist) when describing the clusters (by default 80)}
\item{seed}{Seed to obtain the same results in successive Kmeans (by default 12345)}
\item{...}{other arguments from other methods}
}

\value{
 Returns a list including:
 \item{data.clust}{the active lexical table used in LexCA plus a new column called Clust_ containing the partition}
 \item{coord.clust}{coordinates table issued from CA plus a new column called Clust_ containing the partition}
 \item{centers}{coordinates of centers} 
  \item{clust.count}{counts of documents/words belonging to each cluster and contribution of the clusters to the variability decomposition}
 \item{clust.content}{list of the document/word labels according to the cluster they belong to}
 \item{call}{list of internal objects. \code{call$t} giving the results for the hierarchical tree;\cr See the first reference for more details}
 \item{description}{$desc.axes for description of the clusters by the characteristic axes ($axes) and eta-squared between axes and clusters ($quanti.var).

 $des.cluster.doc for description of the clusters of documents by their characteristic words ($word), supplementary words ($wordsup) and, 
 if contextual variables were considered in LexCA, description of the partition/clusters by qualitative ($qualisup) and quantitative ($quantisup) variables, paragons ($para) and specific words ($dist) of each cluster.

  $des.word.doc description of the clusters of words by their characteristic documents ($docs), paragons ($para) and 
  specific documents ($dist) of each cluster.
}


Returns the hierarchical tree with a barplot of the successive inertia gains, and the 
first CA map of the documents/words with the labels colored according to their cluster.
}


\details{
LexHCca starts from the documents/words coordinates issued from correspondence analysis axes. 
Euclidean metric and Ward method are used. 

The number of clusters is determined either a priori or from the hierarchical tree structure. 
If nb.clust=0 (or 'auto'), a level for cutting the tree is automatically suggested. 
This is computed in the following way, reading the tree downward. 
At a given step, the tree could be cut into Q clusters (Q varying between min and max). 
The between-inertia gain when passing from Q-1 to Q clusters and the between-inertia gain when passing from Q to Q+1 clusters are computed. 
The suggested level corresponds to the maximum value of the ratio between the former and the latter of these inertia-gains.
Note that the between-inertia gain when passing from Q to Q+1 clusters is equal to the value of the Ward criterion when passing from Q+1 to Q clusters when building the tree bottom up.
In this latter case, a level where to cut the tree is suggested. The user can choose to cut the tree at this level or at another one.
 
The results include a thorough description of the clusters, taking into account contextual variables. Graphs are provided.

Frequencies before TextData selection are the marginal frequencies for documents in description analysis.

}



\references{ 

Bécue-Bertaut M. Textual Data Science with R. Chapman & Hall/CRC.  \doi{10.1201/9781315212661}.

Husson F., Lê S., Pagès J. (2011). Exploratory Multivariate Analysis by Example Using R. Chapman & Hall/CRC. \doi{10.1201/b10345}.

Lebart, L., Salem, A., & Berry, L. (1998). Exploring textual data. (D. Kluwer, Ed.). \doi{10.1007/978-94-017-1525-6}.
}

\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Monica Bécue-Bertaut, Josep-Anton Sánchez-Espigares}

\seealso{ \code{\link{LexCA}}, \code{\link{plot.LexHCca}}}

\examples{
data(open.question)
res.TD<-TextData(open.question, var.text=c(9,10), Fmin=10, Dmin=10, stop.word.tm=TRUE,
        context.quali=c("Gender","Age_Group","Education"), context.quanti=c("Age"))
res.LexCA<-LexCA(res.TD, graph=FALSE, ncp=8)
res.hcca<-LexHCca(res.LexCA, graph=TRUE, nb.clust=5)
}


\keyword{multivariate}


