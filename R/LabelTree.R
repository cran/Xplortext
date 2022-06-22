#' @export

LabelTree<-function(object,proba=0.05){
  options(stringsAsFactors = FALSE)
  
    if (!inherits(object, "LexHCca") & !inherits(object, "LexCHCca"))
      stop("Object should be LexHCca or LexCHCca class")
mots.h=function(LTab){
	if (length(LTab)==0) stop("No hierarchical words are identified")
	table=NULL
	for (i in 1:length(LTab)){
		tab=cbind(LTab[[i]],rep(i,nrow(LTab[[i]])))
		tab=tab[tab[,6]>0,,drop=F]
		nuevas=!row.names(tab)%in%row.names(table)
		for (pal in row.names(tab)[!nuevas]){
			if (table[pal,6]<tab[pal,6]) table[pal,]=tab[pal,]
		                                     }
		table=rbind(table,tab[nuevas,,drop=F])            
	                        }
	
	nam=factor(table[,7],levels=1:length(LTab),labels=names(LTab))
      b=split(as.data.frame(table[,1:6]),nam,drop=F)
      lisa<-lapply(b, function(el) if(nrow(el)==0) NULL else el[order(el[,6],decreasing=TRUE),])    
      for(i in sort(which(unlist(lapply(lisa,is.null))),decreasing=T)) lisa[[i]]<-NULL
      return(lisa) 
}
       

## 


###############  <----------------- Comprobar la siguiente, no funciona

if (inherits(object, "LexCHCca")) dendro<-object$dendro$clust
else {
  mat.merge <- object$call$t$tree$merge
  dendro<- list()
  for(i in 1:nrow(mat.merge)) {
    if(mat.merge[i,1]<0) V1 <- -mat.merge[i,1] else V1 <- unlist(dendro[mat.merge[i,1]])
    if(mat.merge[i,2]<0) V2 <- -mat.merge[i,2] else V2 <-unlist(dendro[mat.merge[i,2]])
    a1 <- list(V1) ; a2 <- list(V2)
    dendro[[i]] <- c(a1,a2)
  }
}

DocTerm<-object$data.clust[,1:(dim(object$data.clust)[2]-1)]
nodt<-nrow(DocTerm)  
lista_tabmots<-list()
#lista_tabmots<-structure((vector(mode= "list", length=(nodt*2-2)),
#              names =rep(" ",(nodt*2-2) ))

# nodt number of documents or categories

inodA<-0
for (inode in 1:(nodt-2)) 
 {
  nodecours<-(nodt)-(inode+1)   
  don<-unlist(dendro[[nodecours]])
    label.node<-paste(rownames(DocTerm)[don],collapse=" ")
    by.quali<-rep(2,dim(DocTerm)[1])
    by.quali[don]<-1
    res<-descfreq(DocTerm,by.quali=by.quali,proba=proba)
    if (!is.null(res[[1]]))  {
    inodA<-inodA+1
    lista_tabmots[[inodA]]<- as.data.frame(res[[1]])
    names(lista_tabmots)[inodA]<-label.node 
                              } 
   }
res<-descfreq(DocTerm,proba=proba)


for (inode in 1:nodt) {                         

   label.node<-rownames(DocTerm)[inode]
   if( !is.null(res[[inode]]))  {
   inodA<-inodA+1  
   lista_tabmots[[inodA]]<- as.data.frame(res[[inode]])
   names(lista_tabmots)[inodA]<-label.node  
                                  }    
} 



#  only the non-doubled are kept
res<-NULL
   hierWord <- mots.h(lista_tabmots)
   class(hierWord) = c("LabelTree")
      return(hierWord)
}
