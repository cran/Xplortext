\name{plot.LexCHCca}
\encoding{UTF-8}
\alias{plot.LexCHCca}
\title{Plots for Chronological Constrained Hierarchical Clustering from LexCHCca Objects}

\description{Plots graphs from LexCHCca results: tree, barplot of the aggregation criterion values and first CA map with the documents colored in accordance with the cluster.

}

\usage{
\method{plot}{LexCHCca}(x, axes=c(1, 2), type=c("tree","map","bar"), rect=TRUE, title=NULL, 
  ind.names=TRUE, new.plot=FALSE, max.plot=15, tree.barplot=TRUE,...)  
}
\arguments{
  \item{x}{object of LexCHCca class}
  \item{axes}{length-2 vector defining the axes of the CA map to plot (by default (1,2))}
  \item{type}{type of graph. "tree" plots the tree; "bar" plots the barplot of the successive 
    values of the aggregation criterion (downward reading of the tree); "map" plots the CA map where the individuals 
    are colored in accordances with the cluster of belonging (by default "tree")}
  \item{rect}{if TRUE, when choice="tree" rectangles are drawn around the clusters (by default TRUE)}
  \item{title}{title of the graph. If NULL, a title is automatically defined (by default NULL)}
  \item{ind.names}{if TRUE, the document labels are written on the CA map (by default TRUE)}
  \item{new.plot}{if TRUE, a new window is opened (by default FALSE)}
  \item{max.plot}{maximum of bars in the bar plot of the aggregation criterion (by default 15)}
  \item{tree.barplot}{if TRUE, the barplot of intra inertia losses is added on the tree graph (by default TRUE)}
  \item{\dots}{further arguments passed from other methods...}
}


\value{
Returns the chosen plot
}

\author{Mónica Bécue-Bertaut, Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Josep-Anton Sánchez-Espigares}
\seealso{\code{\link{LexCHCca}} }


\examples{
\dontrun{
data(open.question)
res.TD<-TextData(open.question,var.text=c(9,10), var.agg="Age_Group", Fmin=10, Dmin=10,
        stop.word.tm=TRUE)
res.LexCA<-LexCA(res.TD, graph=FALSE)
res.chcca<-LexCHCca(res.LexCA, nb.clust=4, min=3, graph=FALSE)
plot(res.chcca, type="tree")
plot(res.chcca, type="map")
plot(res.chcca, type="bar", max.plot=5)
}
}

\keyword{plot}
