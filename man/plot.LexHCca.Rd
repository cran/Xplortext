\name{plot.LexHCca}
\encoding{UTF-8}
\alias{plot.LexHCca}
\title{Plots for Hierarchical Clustering from LexHCca Objects}

\description{Plots graphs from LexHCca results: tree and CA maps with the documents or words colored in accordance with the cluster.}


\usage{
\method{plot}{LexHCca}(x, type="map", plot=c("points", "labels", "centers"), selClust="ALL",
     selInd="ALL",axes=c(1, 2), theme=theme_bw(), palette=NULL, title=NULL, 
     axis.title=NULL, axis.text=NULL, points=NULL, labels=NULL,centers=NULL, 
     traject=NULL, hull=NULL, xlim=NULL, ylim=NULL, hvline=TRUE,...) 
}

\arguments{
  \item{x}{object of LexHCca class}
  \item{type}{type of graph. "map" plots the CA map where the individuals are colored in accordance with the cluster of belonging (by default); "tree" plots the dendrogram if hierarchical method without consolidation is performed from LexHCca. See details}
  \item{plot}{elements to plot for map graph: points, labels, centers, hull or traject; by default "ALL" and points, labels and centers are plotted. Also combinations are allowed, i.e: plot=c(points,centers)}
  \item{selClust}{vector indexes with the numbers of the clusters to plot (by default "ALL")}
  \item{selInd}{vector with the active documents/words to plot (indexes, names or rules; see details; by default "ALL")}
  \item{axes}{length-2 vector indicating the axes of the CA map to plot; by default (1,2)}
  \item{theme}{used to modify the theme settings by ggplot2 package (by default theme_bw())}
  \item{palette}{color palette used to draw the clusters. As many numbers as clusters. See details}
  \item{title}{title of the graph. If NULL, a title is automatically defined (by default NULL).
  Other parameters can be chosem using a list: text, color, size, family, face, just; See details}
  \item{axis.title}{axis titles parameters can be used: text.x, text.y, color, size, family, face,
  just;  If text.x and text.y are NULL automatic texts are plotted (by default NULL). See details}
  \item{axis.text}{format of numbers can be chosen: color, size, family, face; See details}
  \item{hvline}{horizontal (intercept.y) and vertical line (intercept.x) added by default at (0,0) position in map. Parameters: intercept.y, intercept.x, linetype (by default "dashed"), color, size, alpha.t}
   \item{points}{format of points. Parameters: size (if size=0 the points are no plotted), shape (by default 21),
   fill (if a color, the same for all the points, if color is NULL palette colors used for the clusters are applied; if more than one color use palette argument; only for shapes from 21 to 25 to fill the point), stroke (controls the edge of the point (by default 0 no edge), border (color of the border, same specifications than fill), alpha.t (by default 1). See geom_point() in ggplot2 library. See details}
   \item{labels}{format of labels. Parameters: size (if size=0 the labels are not plotted; by default 4), family, face, hjust, vjust, color.text, alpha.t.text, numbers(if TRUE the label will be replaced by the number of the cluster to which it belongs, by default FALSE), rect (if TRUE a rectangle will be drawn around the label, by default FALSE), color.fill (color into the rectangle, by default FALSE is transparent), alpha.t.fill, force (to do repulsive textual annotations and make it easier to read), max.overlaps (maximum number of overlapped points, by default 10, can be Inf)
   set.seed (by default a new seed for each plot draws different positions, for the same seed i.e: set.seed=1234)}
 \item{traject}{draws trajectory arrows in accordance with the order of clusters or in the selInd order. Parameters: color (by default blue), linetype (by default 1 solid), space (by default 0 and no space is added from point to arrow, be careful with this value), size (width,by default 1), arrow.length (of the arrow, by default .3), arrow.type (by default "closed"), arrow.angle (by default 30), alpha.t. See geom_segment for details}
   \item{centers}{draws the barycenter of the clusters. Parameters: size (by default 5), family, face, color (of the border, only one), fill, alpha.t, labels (string vector with the names of the clusters)}
 \item{hull}{draws a hull containing all the elements of each cluster. Parameters: type (ellipse, by default, hull), alpha.t, color, linetype (by default "dotted"))}
 \item{xlim}{pair of values xlim=c(xmin,xmax). If a NA value, this limit is automatically calculated}
 \item{ylim}{pair of values ylim=c(ymin,ymax). If a NA value, this limit is automatically calculated}
 \item{...}{other arguments from other methods}
}

\details{

Parameter type="tree" shows the dendrogram
\preformatted{
 - if hierarchical cluster without consolidation is performed.
 - if hierarchical cluster with consolidation before the consolidation.
 - if kmeans the hierarchical tree with the output of kmeans.
}
You can make customer dendrograms by accessing the hclust format object located inside the object in
hclust format from object$call$t$tree

Selection of individuals (documents or words) to plot:

1. Using labels:
\preformatted{
selInd = c("doc1","doc5"): only the documents with labels doc1 and doc5 are plotted.
}

2. Using indexes:
\preformatted{
selInd = c(1:5): cases 1 to 5 are plotted.
   }


3. Using rules:
\preformatted{
 Rules are based on the coordinates (coord), the contribution (contrib or meta; 
 concerning only active elements) or the square cosine (cos2). }
 
 Somes examples hereafter:

\preformatted{
selInd="coord 10": only the 10 cases with the highest coordinates, as globally
   computed on the 2 axes, are plotted.
selInd="contrib 10": the cases with a contribution to the inertia, of any of 
   the 2 axes over 10 percent.
selInd="meta 3": the cases with a contribution over 3 times the average word/document 
   contribution on any of the two axes are plotted.
selInd="cos2 .85": the documents with a cos2 over 0.85, as summed on the 2 axes, 
   are plotted.
}

Parameters can be used in combination, e.g.: title=c("text"="CA", "color"="red").

See grDevices package (The R Graphics Devices and Support for Colours and Fonts).

palette, the color of the palette used to draw the points. 
   By default colors are chosen. If you want to define the colors for three clusters : palette=c("black","red","blue"); 
   or you can use: palette=palette(rainbow(30)); or in black and white for example: palette=palette(gray(seq(0,.9,len=25))).

Family Fonts (family). Also see the extrafont package for a much better support of fonts: library(extrafont); font_import(). By default "family"='serif'.

Face fonts (face). Can be 'plain', 'bold', 'italic', 'bold.italic', 'symbol'. By default 'plain'.

alpha.t is the level of transparency for some objects. 0 value means full transparency and 1 opacity. By default 1.

Values for horizontal justification hjust, vertical vjust and both hvjust can be (c,centered or 0.5 if centered; l,left or 0 if left; r, right or 1 if right)

\preformatted{
By default in:
* title: text="Clusters on the CA map"; color=black; size=18; familiy=serif; face=plain; 
      hjust=0.5.

* axis titles:  text.x=Dim x (\%), text.y=Dim y (\%), color=black, size=12, family=serif, 
      face=plain, just=centered.
  
* axis.text: color=black, size=8, family=serif, face=plain.

* hvline: intercept.x=0, intercept.y=0, linetype=dashed, color=gray, size=0.5, alpha.t=1.
  
* points: size=2, shape=21, border:automatic cluster color, fill:automatic cluster color, 
      stroke=0, border: automatic cluster color, alpha.t=1.

* labels: size=4, family=serif, face=plain, hjust=1, vjust=1, color.text=same of points, 
      alpha.t.text=1, numbers=FALSE, rect=FALSE, color.fill=transparent, alpha.t.fill=1,
      force=1, max.overlaps=10.

* traject: color=blue, linetype=solid, space=1, arrow.length=.3, arrow.type= closed, 
      arrow.angle=30, alpha.t=1. 

* centers: size=5, family=serif, face=italic, color, fill=automatic cluster color,
      alpha.t=1, labels=automatic strig vector with the names of the clusters.

* hull: type=ellipse, alpha.t=0.1, color=black, linetype=dotted }
}


\author{Ramón Alvarez-Esteban \email{ramon.alvarez@unileon.es}, Mónica Bécue-Bertaut, Josep-Anton Sánchez-Espigares}

\seealso{\code{\link{LexHCca}}}

\references{
The Xplortext web site provides several examples at <https://xplortext.unileon.es/?page_id=766>.
}

\examples{
data(open.question)
res.TD<-TextData(open.question,var.text=c(9,10), var.agg="Age_Group", Fmin=10, Dmin=10,
        stop.word.tm=TRUE)
res.LexCA<-LexCA(res.TD, graph=FALSE)
res.chcca<-LexHCca(res.LexCA, nb.clust=4, min=3, graph=FALSE)
plot(res.chcca, type="tree")
plot(res.chcca, type="map")
}

\keyword{plot}
