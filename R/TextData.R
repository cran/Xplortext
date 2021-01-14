#' @import FactoMineR slam stringr tm graphics gridExtra utils stringi
#' @rawNamespace import(stats, except = c(hclust))
#' @importFrom utils packageDescription
#' @export
TextData <- function (base, var.text=NULL, var.agg=NULL, context.quali=NULL, context.quanti= NULL, 
    selDoc="ALL", lower=TRUE, remov.number=TRUE, lminword=1, Fmin=Dmin, Dmin=1, Fmax=Inf,
    stop.word.tm=FALSE, idiom="en", stop.word.user=NULL, segment=FALSE, 
    sep.weak="(['?]|[[:punct:]]|[[:space:]]|[[:cntrl:]])+",
    sep.strong="\u005B()\u00BF?./:\u00A1!=+;{}-\u005D", 
    seg.nfreq=10, seg.nfreq2=10, seg.nfreq3=10,
    graph=FALSE)
{
  
  dfold <- deparse(substitute(base))

  
# library(SnowballC)   docs <- tm_map(docs, stemDocument)
## REvisar el print algunas palabras tienen sin cabecera
# filt = "(['?]|[[:punct:]]|[[:space:]]|[[:cntrl:]])+"
filt=sep.weak

if(!is.null(var.agg)) if(is.character(base[,var.agg])) base[,var.agg] <- as.factor(base[,var.agg] ) # version 1.3.1	

#---------------------------------------------------
plotTextData <- function()
{
# if(dev.interactive()) dev.new()
plot.TextData(y)
}

#---------------------------------------------------
# Count occurrences		
occurrFunc <- function(z,title, dOcc, bagg){		
DT <- data.frame(z$i,z$j,z$v)	
colnames(DT) <- c("i","j","v") 	
a1 <- data.frame(sapply(split(DT, DT$i, drop=FALSE), function(w) sum(w$v)))	
q <-4			
if(bagg==TRUE) if(!is.null(dOcc )){q <- q+1 }		
if(is.null(dOcc)){ dOcc <- cbind(data.frame(z$dimnames$Docs,0)); q<-0}		
dOcc[,(q+2)] <- 0	
dOcc[rownames(a1),(q+2)] <- a1[1]	
a2 <- data.frame(sapply(split(DT, DT$i, drop=FALSE), function(w) length(w$v))) 	
dOcc <- cbind(dOcc,0)		
dOcc[,(q+3)] <- 0		
dOcc[rownames(a1),(q+3)] <- a2[1]		
if(q==0) colnames(dOcc)[1] <- c("DocName")	
colnames(dOcc)[q+2]<- c( paste0("Occurrences.",title,sep="" ))		
colnames(dOcc)[q+3]<- c( paste0("DistinctWords.",title,sep=""))		
return(dOcc)}		

#---------------------------------------------------			
aggreg <- function(text.var, grouping.var){			
 G <- as.character(substitute(grouping.var))			
 grouping <- unlist(grouping.var)			
 y <- rle(as.character(as.vector(grouping)))			
 lens <- y$lengths			
 group <- y$values			
 x <- cumsum(lens)			
 st <- c(1, x[-length(x)] + 1)			
 end <- c(x)			
 L1 <- invisible(lapply(seq_along(st), function(i) {			
	pasteagg(text.var[st[i]:end[i]], sep = " . ")}))			
 names(L1) <- group			
 DF <- data.frame(x = names(L1), text.var = unlist(L1), row.names = NULL)			
 colnames(DF)[1] <- "Group"			
return(DF)}			

#---------------------------------------------------
pasteagg <- function(mc,sep = ".") 
{
mc <- lapply(mc, function(x) {
 gsub("^\\s+|\\s+$", "", x)})
mc <- do.call("cbind", mc)
m <- { apply(mc, 1, function(x) {
 if (any(is.na(x))) { NA }
   else {
    paste(x, collapse = sep)}})}
names(m) <- NULL
return(m)
}

#---------------------------------------------------		
# Function to recodify the position of the words		
recoderFunc <- function(z, from, to) return(to[match(z, from)])		

#---------------------------------------------------		
#Function to select words		
selectFunc <- function(z, selwords) {		
pos<- which(z$j%in%sel.words)		
z$j<- z$j[pos]; z$i<- z$i[pos]; z$v<- z$v[pos]		
z$j<-recoderFunc(z$j,selwords,1:length(selwords))		
z$j<- as.numeric(factor(z$j,labels=c(1:length(selwords))))		
z$dimnames$Terms<-z$dimnames$Terms[selwords]		
z$ncol<-length(z$dimnames$Terms)		
z$nrow<-length(z$dimnames$Docs)		
return(z)}		

#---------------------------------------------------					
infoNew <- function(){
mbase <- list(dfold, "name of input R data.frame")
menvir <- list(globalenv(),"name of environment")
mvartext <- list(var.text, "names of textual columns")
midiom <- list(idiom, "idiom of the corpus, (by default en)")
mlminword <- list(lminword, "minimum length of a word (by default 1)")					
mlower <- list(lower, "converting the corpus into lowercase (by default TRUE)")					
mremnum <- list(remov.number, "removing the numbers (by default TRUE)")            					
mFmin <- list(Fmin, "minimum frequency of a word (by default Dmin)")					
if(is.null(Fmax)) Fmax <- Inf					
mFmax <- list(Fmax, "maximum frequency of a word (by default Inf)")					
mDmin <-  list(Dmin, "minimum number of documents using a word (by default 1)")		
mndoc <- list(dtm$nrow, "number of non-empty non-aggregate documents") 	
mnlength <- list(sum(Nfreqword), "corpus size (total number of occurences)")	
mnWord <- list(dtm$ncol, "vocabulary size (total number of words)")
if(is.null(var.agg)) var.agg <- "" 	
mndocsagg <- list(var.agg, "name of the aggregation variable")
stop_word_tm <- stop.word.tm
if(stop.word.tm==TRUE) 			
 stop_word_tm <- stopwords(kind=idiom)
mSwtm <-  list(stop_word_tm, "stopword list provided by tm package")
mSwuser <-  list(stop.word.user, "stopword list provided by the user")
mbsegm <- list(segment, "searching repeated segments (by default FALSE)")
zinfo <- list(base=mbase,menvir=menvir, var.text=mvartext, idiom = midiom, lminword=mlminword, lower=mlower,remov.number=mremnum, Fmin=mFmin,
Fmax=mFmax,Dmin=mDmin,Ndoc=mndoc, LengthW= mnlength, Nword=mnWord,
name.var.agg=mndocsagg ,stop.word.tm=mSwtm, stop.word.user =mSwuser, segment.searched = mbsegm)


if(segment==TRUE) { 
 mNseg <- list( nbseg, "number of segments")
 mseg.nfreq <- list(nfreq, "minimum frequency of a more-than-3-words repeated segment (by default 10)")			
 mseg.nfreq2 <- list(nfreq2, "minimum frequency of a length-two repeated segment (by default 10)")					
 mseg.nfreq3 <- list(nfreq3, "minimum frequency of a length-three repeated segment (by default 10)")					
 segments <- list(Nseg=mNseg, seg.nfreq=mseg.nfreq, seg.nfreq2=mseg.nfreq2 ,seg.nfreq3=mseg.nfreq3) 
 zinfo<- c(zinfo, segments=list(segments))
 }
return(zinfo)
}

 nxlon<-20
 nfreq <- seg.nfreq; nfreq2<- seg.nfreq2; nfreq3<-seg.nfreq3 					
if(segment==FALSE) { nfreq <- ""; nfreq2<-""; nfreq3<-""}	

blongErr <- FALSE				
								
if(segment==TRUE) {					
# auxiliary functions
REPWEAK <-function(chaine,sep.weak) res<-stringr::str_replace_all(chaine,sep.weak, " ")					
REPSTRONG <-function(chaine,sep.strong) res<-stringr::str_replace_all(chaine,sep.strong, " zzwwxxyystr ")					
PROCHE <-function(ideb,ifin,ITEX,ITDR,ITRE,nfreq,nfreq2,nfreq3,long,nxlon,nbseg)					
 {   					
# the function proche detects the first sublist of adresses in ITDR corresponding a same successor					
# if this successor is not "end of answer" or "strong separator", we have located an admissible segment					
 list.segment<-list()				
 ad.segment<-vector()					
 te.segment<-NULL					
 rep.segment<-vector()					
 mfrec<-0					
 ipunt<-ideb-1					
 isucc<-ITEX[ITDR[ideb]+long-1] 					
 while(  (ITEX[ITDR[ipunt+1]+long-1]==isucc) &  (ipunt < ifin) )					
   {					
    if (!( (isucc=="zzwwxxyystr") | (isucc=="zzwwxxyyendrep")))					
      {					
       mfrec<-mfrec+1					
       ad.segment[mfrec]<-ITDR[ipunt+1]					
       rep.segment[mfrec]<-min(which(ITRE>ad.segment[mfrec]))					
       }					
   ipunt<-ipunt + 1					
    }									
 ifin<-ipunt					
 nfreq.threshold<-nfreq					
 if (long==1)  nfreq.threshold<-999999999999					
 if (long==2)  nfreq.threshold<-nfreq2					
 if (long==3)  nfreq.threshold<-nfreq3					
 ltrou<-( !( (isucc=="zzwwxxyystr") | (isucc=="zzwwxxyyendrep"))  & (mfrec >= nfreq)) 					
 ltrouseg<-( ltrou & (mfrec >= nfreq.threshold))  					
					
### here we have to control that it is not a "constrained segment", either at the left part or at the right part					
 if (ltrouseg) 					
  { 					
   contraintG<-TRUE  					
   contraintD<-TRUE								
   IGAUC<-ad.segment-1					
   special<-c("zzwwxxyystr", "zzwwxxyyendrep")					
   if (0 %in% IGAUC) contraintG<-FALSE 					
   if (max(special %in% ITEX[IGAUC])) contraintG<-FALSE    					
   if (contraintG) {					
     SUCC<-as.factor(ITEX[IGAUC])					
     if (nlevels(SUCC) > 1) contraintG <- FALSE 					
                    }					
     IDROI<-ad.segment+long					
     if (max(special %in% ITEX[IDROI])) contraintD<-FALSE    					
     if (contraintD) {					
         SUCC<-as.factor(ITEX[IDROI])					
         if (nlevels(SUCC) > 1) contraintD <- FALSE 					
                       }					
     if (contraintD | contraintG) ltrouseg<-FALSE					
     if (contraintG) ltrou<-FALSE      					
  }					
###  final contrained segments removing					
 if (ltrouseg) 					
   { 					
    for (i in 1:long){					
        te.segment<-paste(te.segment,ITEX[ITDR[ideb]+(i-1)],sep=" ")				
               te.segment<-stringr::str_trim(te.segment)}					
   lo.segment<-long				
   fr.segment<-mfrec
   nr.segment<-nbseg+1					
   list.segment<-list(te.segment,fr.segment,ad.segment,rep.segment,lo.segment,nr.segment)					
       names(list.segment)<-c("text","frequency","adresses","documents","length","nr.seg")					
          }       					
       return(list(ifin=ifin,ltrou=ltrou,ltrouseg=ltrouseg,list.segment=list.segment))					
   }                                                                                                                                                                                                                                                                                                                                                                                                                                                					
					
ORD.EXT<-function(ICRIT,ADR,long1)					
 {					
# Ordering adresses from successors in text					
 ICRIT_ord<-order(ICRIT)					
 ADR<-ADR[c(ICRIT_ord)]					
 return(list(ADR=ADR))					
 }					
}					
##### Final internal functions
# dfold <- deparse(substitute(base))  # Movido al principio 15/06/2020






if(!is.null(var.agg)){
  if(is.numeric(base[,var.agg])) base[,var.agg] <- as.factor(base[,var.agg])
# Añadida la siguiente el 10/01/2019 para eliminar niveles de factores no utilizados en variable de agrupación
  base[,var.agg] <- factor(base[,var.agg])
  var.agg.seg <- data.frame(base[,var.agg,drop=FALSE])
   }
  remov.docs <- rownames(base)
  
#--------- Selecting docs by rownumber or rowname -------------------					
if(selDoc!="ALL") {					
 if (!is.character(selDoc)) 					
  selDoc <- rownames(base)[selDoc]					
  selDoc <- which(rownames(base) %in% selDoc)					
  base <- base[selDoc,]}					
					
#--------- Save corpus var.text  -------------------					
 if(!is.character(var.text)) var.text <- colnames(base)[var.text]
 var.text <- var.text[which(var.text %in% colnames(base))]
 if(min(var.text)<1) stop("Error in var.text")					
    if(length(var.text) == 0) stop("You must define var.text")					
       corpus <- base[, var.text[1]]					
 if(length(var.text) > 1) {					
   for (i in 2:length(var.text)){					
      corpus <- paste(corpus, base[, var.text[i]], sep = ".")}}					
   corpus <- data.frame(corpus, stringsAsFactors = FALSE)					
 rownames(corpus) <- rownames(base)					


#--------- Save context.quanti  -------------------
data.context.quanti <- NULL
if(!is.null(context.quanti)){ 
 if(length(context.quanti)==1) {
   if(context.quanti=="ALL") 
    context.quanti <- names(which(sapply(base,is.numeric)))  }
 if (!is.character(context.quanti)) 
   context.quanti <- colnames(base[context.quanti])
 context.quanti <- context.quanti[which(context.quanti %in% colnames(base))]
 tq <- names(which(sapply(base,is.numeric)))
 pos <- which(context.quanti %in% tq)
 data.context.quanti <- data.frame(base[,context.quanti[pos]])
 colnames(data.context.quanti) <- context.quanti[pos]
 qf <- context.quanti[-pos]
   if(!is.null(qf))
    if(length(qf)>0)
      { # There are quantitative variables as factors
    for (i in 1:length(qf)) {
    levels(base[qf[i]]) <- levels(droplevels(base[qf[i]]))  
    rdo <- suppressWarnings(as.numeric(levels(base[,qf[i]])))
    valrdo <- rdo[!is.na(rdo)]
    if(length(valrdo)>0) {
      rdo <- suppressWarnings(as.numeric(as.character(base[,context.quanti[i]])))
        if(is.null(data.context.quanti)) {	
          data.context.quanti <- data.frame(rdo) 	
          colnames(data.context.quanti) <- context.quanti[i]	
        } else { 	
        data.context.quanti <- cbind(data.context.quanti, rdo) 	
        colnames(data.context.quanti)[ncol(data.context.quanti)] <- context.quanti[i]	
}}}} 
if(!is.null(data.context.quanti)) rownames(data.context.quanti) <- rownames(base)
 }

var.check <- NULL

#--------- Save context.quali  -------------------
tmpquali <- NULL
if(!is.null(context.quali)) {
 if(length(context.quali)==1) {
   if(context.quali=="ALL") 
    context.quali <- names(which(sapply(base,is.factor)))
}
 if (!is.character(context.quali)) 	
  context.quali <- colnames(base[context.quali])	
  context.quali <- context.quali[which(context.quali %in% colnames(base))]
   tq <- names(which(sapply(base,is.numeric)))
   pos <- which(context.quali %in% tq)
   if(length(pos)>0)  
    base[,context.quali[pos]] <- as.factor(as.character(base[,context.quali[pos]]))

# Remove selected text variables, var.text from context.quali
  pos <- which(context.quali %in% var.text)
  if(length(pos)>0) context.quali <- context.quali[-pos]
}

# Remove var.agg context.quali
if(!is.null(var.agg)){
  if(is.numeric(base[,var.agg])) base[,var.agg] <- as.factor(base[,var.agg])
if(!is.null(context.quali)) {  		
# Remove var.agg from qualitative context		
  pos <- which(context.quali %in% var.agg)		
  if(length(pos)>0) context.quali <- context.quali[-pos]		
  var.check <- c(context.quali, var.agg)			
}}		

if(!is.null(context.quali)) var.check <- context.quali
if(!is.null(var.agg)) var.check <- c(context.quali, var.agg)
if(length(var.agg)>1) stop("Only one variable for aggregation")

#--------- Rename NA levels in factor var.context and var.agg variables To Missing
nvcheck <- length(context.quali)
for(i in 1:nvcheck) {
 labi <- levels(base[,context.quali[i]])  # Levels of factors
# If any value is NA but there is not a level NA
 if(any(is.na(base[,context.quali[i]]))){
  levels(base[,context.quali[i]]) <- c(labi,"Missing")
  pos <- which(is.na(base[,context.quali[i]]))
  base[pos,context.quali[i]] <- "Missing"}
# If any value is '' but there is not a level NA
 pos <- which(labi %in% '')
 if(!is.null(pos)>0) levels(base[,context.quali[i]])[pos] <- "Missing"
# If any value is <NA> but there is not a level <NA>
 pos <- which(labi %in% '<NA>')
 if(length(pos)>0) levels(base[,context.quali[i]])[pos] <- "Missing"
}

#--------- Rename repeated levels in factor var.context and var.agg variables
if(nvcheck >1){
for(i in 1:(nvcheck -1)) {
 strnamei <- var.check[i]
 levi <- levels(base[,strnamei])
 for(j in (i+1):nvcheck) {
  strnamej <- var.check[j]
  levj <- levels(base[,var.check[j]])
  repetij <- levi[(which(levi %in% levj))]
  nrep <- length(repetij)
    if(nrep>0){
     missrep <- which("Missing" %in% repetij)
   if(missrep==1) nrep <- nrep-1
      if(nrep>0){
        levels(base[,strnamei]) <- paste0(strnamei,"_",levi)
        levels(base[,strnamej]) <- paste0(strnamej,"_",levj)
}}}}}


# Rename "Missing" to "Missing" & variable name
for(i in 1:(nvcheck)) {
  strnamei <- var.check[i]
  levi <- levels(base[,strnamei])
  miss <- which(levi %in% "Missing")
  levels(base[,strnamei])[miss] <- paste0("Missing","_",strnamei)
}

if(!is.null(var.agg)){
  dfvaragg <- data.frame(base[,var.agg,drop=FALSE])
}

#if(packageDescription("tm")$Version >"0.7-1") {
 colnames(corpus)[1] <- "text"
 corpus$doc_id <- rownames(corpus)
#}

#--------- Read texts from tm -------------------
# dtmCorpus <- Corpus(DataframeSource(corpus), readerControl = list(language = idiom))
dtmCorpus <- VCorpus(DataframeSource(corpus), readerControl = list(language = idiom)) 
dtmCorpus <- tm_map(dtmCorpus, content_transformer(function(x) gsub(filt, " ", x)))
dtmCorpus <- tm_map(dtmCorpus, stripWhitespace)
dtm <- DocumentTermMatrix(dtmCorpus, control = list(tolower = lower, wordLengths = c(lminword, Inf)))
rownames(dtm) <- rownames(base)

if(!is.null(var.agg)) SourceTerm <- dtm

# ---------------- If aggregation ---------
if(!is.null(var.agg)){
# To build a data frame with 3 columns (rows, columns and frequency)
  # This is a compressed table
 agg <- data.frame(base[dtm$i,var.agg], dtm$j,dtm$v,drop=FALSE)
 agg <-aggregate(agg[,3], by=list(agg[,1], agg[,2]),FUN=sum, na.rm=TRUE)
 agg <- agg[order(agg[,1],agg[,2]),]
 dtmagg <- dtm

 agg[,1] <- droplevels(agg[,1])
 dtmagg$nrow <- length(levels(agg[,1]))
 dtmagg$i <- as.numeric(agg[,1])
 dtmagg$j <- agg[,2]
 dtmagg$v <- agg[,3]
 dtmagg$dimnames$Docs <- levels(agg[,1])
 detOccAgg <- occurrFunc(dtmagg, "before",NULL, TRUE) 
 }
#--------- ------------------			
Nfreqword <-tapply(dtm$v,dtm$j,sum)			
Ndocword  <-tapply(dtm$v>0,dtm$j,sum)			
Table <- cbind(Nfreqword,Ndocword)			
rownames(Table) <- dtm$dimnames$Terms			
colnames(Table) <- c("Frequency", "N.Documents")			
TFreq <- Table[order(Nfreqword, Ndocword, decreasing = TRUE), ]					
detOcc <- occurrFunc(dtm, "before",NULL, FALSE) 			
ndocIni <- nrow(detOcc)			
ndocIniEmpty <- ndocIni - length(unique(dtm$i))			
rownamesdocs.no.empty <- rownames(base)[unique(dtm$i)]			
rownamesdocs.empty <- rownames(base)[-unique(dtm$i)]			

N <- ndocIni
detOcc$PctLength.before <- 100*detOcc[,2]/sum(detOcc[,2])
detOcc$MeanLength100.before <- round(N*100*detOcc[,2]/sum(detOcc[,2]),2)
detOcc$PctLength.before <- round(detOcc$PctLength.before,2)				
wordsafter <- dtm$ncol			

Docs.before <- detOcc[,c("DocName", "Occurrences.before")]

if(!is.null(var.agg)){			
  numberdocs <- table(base[,var.agg])			
  posic <- which(names(numberdocs) %in% dtmagg$dimnames$Docs) 			
  numberdocs <- numberdocs[posic]			
  detOccAgg$NumberDocs <- numberdocs			
  detOccAgg$PctLength.before <- 100*detOccAgg[,2]/sum(detOccAgg[,2])			
  detOccAgg$MeanLength100.before <- round((detOccAgg[,2]*100/numberdocs)/(sum(detOccAgg[,2])/sum(numberdocs)),2)			
  detOccAgg$PctLength.before <- round(detOccAgg$PctLength.before,2)			
}			


if(!is.null(var.agg)) {	
if(!is.null(corpus$doc_id)) corpus$doc_id <- NULL
 corpusSeg <- corpus
 corpus <- cbind(corpus,base[,var.agg]) 		
 corpus <- corpus[order(corpus[,2]), ]		
 names(corpus)[1] <- "my_text"		
 names(corpus)[2] <- "my_agg"	
 corpus <- aggreg(text.var=corpus$my_text, grouping.var = corpus$my_agg)		
 rownames(corpus) <- corpus[,1]		
 names(corpus)[2] <- "corpus"		
 corpus[,1] <- NULL 
}		

if(segment==TRUE) {						
 maj.in.min = lower # y$info[3,2]						
 sep.weak ="([\u0027\u02BC\u0060]|[,;'?\n\u202F\u2009\u0028]|[[:punct:]]|[[:space:]]|[[:cntrl:]])+"						
 if (nfreq2<nfreq) nfreq2<-nfreq		
 if (nfreq3<nfreq) nfreq3<-nfreq		
 if (nfreq3>nfreq2)nfreq3<-nfreq2	
 text1<-apply(as.matrix(apply(as.matrix(corpus),1,FUN=REPSTRONG,sep.strong)),1,FUN=REPWEAK,sep.weak)						
 text3<-apply(matrix(text1),1,(stringr::str_c),"zzwwxxyyendrep",sep=" ") 						
 nrep<-NROW(text3)						
 listrep<-strsplit(as.character(text3),split=" ") 						
 ITEX <- unlist(listrep)						
 # ITEX is a vector of occurrences of words 						
 if (maj.in.min == TRUE)  ITEX <- tolower(ITEX)						
 if (remov.number == TRUE) ITEX<- removeNumbers(ITEX)						
 # ITEX is a vector of occurrences of words, but with fictitious "empty" words because of the multiple spaces						
 # these fictitious words have to be eliminated						
	sel <- which(ITEX=="") 					
	if (length(sel)!=0){	 				
    	ITEX <- ITEX[-sel]                					
	}					
 # The text is in the form of a vector of occurrences  						
 ITEX.f<-as.factor(ITEX)						
 FREQ.mots<-table(ITEX.f)						
 FREQ.cum<-cumsum(FREQ.mots)						
 Vplus<-dim(FREQ.mots)										
 # To conserve the addresses when ordering ITEX						
 ITDR<-order(ITEX)
 # adress of the answers (=adress of the first word corresponding to the answer) in ITEX						
 ITRE<-vector()						
 ITRE<-which(ITEX=="zzwwxxyyendrep")											
#######  the data structures are built						
 Nplus<-length(ITEX)						
 ITDR<-seq(1,Nplus,1)						
 lpil<-vector()						
 list.tot.segment<-list()						
 # global initialisations						
 ideb<-1						
 ifin<-Nplus						
 long<-0						
 nbseg<-0						
# for all the distinct words, we have to detect the segments beginning with this word						
 ltrou<-((ifin-ideb+1) >= nfreq)  											
while(ltrou)      						
 {                						
  while (ltrou & (long<=nxlon))       #exploration of the possible segments issued from word_in_course	
     {						
  if(long>nxlon) blongErr <- TRUE	
      long1<-long						
      long<-long+1					
      lpil[long]<-ifin						
      res.ORD.EXT<-ORD.EXT(ITEX[ITDR[ideb:ifin]+long1],ITDR[ideb:ifin],long1)						
      ITDR[ideb:ifin]<-res.ORD.EXT$ADR						
      ltrou<-FALSE						
      res.proch<-PROCHE(ideb,ifin,ITEX,ITDR,ITRE,nfreq,nfreq2,nfreq3,long,nxlon,nbseg)						
      ifin<-res.proch$ifin						
      ltrou<-res.proch$ltrou						
      ltrouseg<-res.proch$ltrouseg						
      if (ltrouseg)						
       {						
        nbseg<-res.proch$list.segment[[6]]						
        list.segment<-res.proch$list.segment						
        list.tot.segment[[nbseg]]<-list.segment 						
        }         						
      }				
       ltrou<-FALSE						
       while (!ltrou & (long>=1) & (ifin<Nplus))						
        {						
         ideb=ifin+1						
         ifin=lpil[long]						
         while (  ( (ideb+nfreq)>ifin) & (long > 1) )						
          {						
              ideb=ifin+1						
              long<-long-1						
              ifin=lpil[long]						
           }						
						
        if (long>=1)						
            {						
              ltrou<-FALSE						
              res.proch<-PROCHE(ideb,ifin,ITEX,ITDR,ITRE,nfreq,nfreq2,nfreq3,long,nxlon,nbseg)						
              ifin<-res.proch$ifin						
              ltrou<-res.proch$ltrou						
              ltrouseg<-res.proch$ltrouseg						
              if (ltrouseg)						
              {						
                 nbseg<-res.proch$list.segment[[6]]						
                 list.segment<-res.proch$list.segment						
                 list.tot.segment[[nbseg]]<-list.segment						
               }                						
             }     						
           }						
        }   						
 # all the segments have been detected and the doc_segments (tab.seg) will be created						
						
 tab.seg<-matrix(0,nrow=nrep,ncol=nbseg)						
 rownames(tab.seg)<-rownames(dtm$DocTerm)						
 if (nbseg==0) print ("\nno segments fullfil the conditions\n")						
						
if (nbseg>0)						
     {						
      for (iseg in 1:nbseg)						
        {						
         list.segment<-list.tot.segment[[iseg]]						
         mfreq<-list.segment[[2]]						
         long.seg<-list.segment[[5]]						
         nseg<-list.segment[[6]]						
         for (i in 1:mfreq) 						
           {						
            rep<-list.segment[[4]][i]						
            tab.seg[rep,nseg]<-tab.seg[rep,nseg]+1						
            }						
         }						
      row.names(tab.seg)<-row.names(dtm$DocTerm$dimnames$Docs)						
      nom.col<-vector()						
      for (iseg in 1:nbseg) nom.col[iseg]<-(list.tot.segment[[iseg]]$text)						
      colnames(tab.seg)<-nom.col						
     }						
						
impri.segment<-data.frame(ncol=3)						
# Segment list in alphabetical ordre						
for (iseg in 1:nbseg) 						
    {						
      impri.segment[iseg,1]<-list.tot.segment[[iseg]]$text						
      impri.segment[iseg,2]<-list.tot.segment[[iseg]]$frequency						
      impri.segment[iseg,3]<-list.tot.segment[[iseg]]$length						
     }						
colnames(impri.segment)<-c("segment","frequency","long")						
segOrderFreq<-with(impri.segment,impri.segment[order(frequency,long,decreasing=TRUE),])						
segOrderlist<-impri.segment						
Index.segments<-list(segOrderFreq=segOrderFreq, segOrderlist=segOrderlist)						
namesSeg<-colnames(tab.seg)						
numSeg<-rep(1:ncol(tab.seg),1)						
colnames(tab.seg) = paste(numSeg, namesSeg, sep=":")						
rownames(tab.seg)<-rownames(dtm$DocTerm)						
}  # Final segments						




#--------- Remove the numbers  ------------------						
# To Detect if the colname is a number and remov.number=TRUE we must remove the column						
if(remov.number == TRUE) {						
 sel.words <- dtm$dimnames$Terms[suppressWarnings(is.na(as.numeric(dtm$dimnames$Terms)))]						
 sel.words <- which(dtm$dimnames$Terms%in%sel.words)						
 if(length(sel.words)>0){ 						
  dtm <- selectFunc(dtm,sel.words)						
  Nfreqword <- Nfreqword[sel.words]}}						
						
#--------- Removing words with low length lminword ------------------						
if (lminword > 1) {
#  sel.words <- which(nchar(dtm$dimnames$Terms) > (lminword-1)) 										
 sel.words <- which(stringi::stri_length(dtm$dimnames$Terms) > (lminword-1)) 
 if(length(sel.words)>0){						
  dtm <- selectFunc(dtm,sel.words)						
  Nfreqword <- Nfreqword[sel.words] }}						
					
#--------- Removing words with low frequency "Fmin" times ------------------						
if (Fmin > 1) {						
 sel.words <- which(Nfreqword > (Fmin-1)) 						
 if(length(sel.words)>0){						
   dtm <- selectFunc(dtm,sel.words)						
   Nfreqword <- Nfreqword[sel.words] }}						
						
#--------- Selecting words appearing with a minimum frequency of "Fmin" times 						
#--------- in a minimum of "Dmin" documents											
Ndocword <-tapply(dtm$v>0,dtm$j,sum)						
if(Fmin>1 | Dmin>1) {						
 sel.words <- which(Nfreqword >= Fmin & Ndocword >= Dmin)						
 if(length(sel.words)>0) dtm <- selectFunc(dtm,sel.words)						
}						
  Nfreqword<-tapply(dtm$v,dtm$j,sum)											
						
#--------- Removing words appearing a maximum frecuency of Fmax times						
if(!is.null(Fmax)) {						
 sel.words <- which(Nfreqword < (Fmax+1))						
 if(length(sel.words)>0) {						
 dtm <- selectFunc(dtm,sel.words)						
 Nfreqword <- Nfreqword[sel.words]}}
						
#--------- Removing stopwords tm (defined in tm package)only if not previously removed						
if(stop.word.tm==TRUE){						
 stop.word <- stopwords(idiom)						
 sel.words <- which(!dtm$dimnames$Terms%in%stop.word)						
 if(length(sel.words)>0) {						
 dtm <- selectFunc(dtm,sel.words)						
 Nfreqword <- Nfreqword[sel.words]}}						
						
#--------- Removing user stopwords						
if(!is.null(stop.word.user)) {						
 if(is.data.frame(stop.word.user)) stop.word.user <- t(stop.word.user)						
  stop.word.user <- stop.word.user[order(stop.word.user)] 						
  sel.words <- which(!dtm$dimnames$Terms%in%stop.word.user)						
  if(length(sel.words)>0){						
  dtm <- selectFunc(dtm,sel.words)						
  Nfreqword <- Nfreqword[sel.words]}}										

docsbefore <- dtm$dimnames$Docs
docsafter <- docsbefore[unique(dtm$i)]
dseldoc <- which(docsbefore %in% docsafter) 
remov.docs <- docsbefore[-dseldoc]


# ---------------- If aggregation ---------				
if(!is.null(var.agg)){				
 agg <- data.frame(base[dtm$i,var.agg], dtm$j,dtm$v)				
 agg <- aggregate(agg[,3], by=list(agg[,1], agg[,2]),FUN=sum, na.rm=TRUE)				
 agg <- agg[order(agg[,1],agg[,2]),]				
 dtmagg <- dtm				

 
# Changed 10/Jan/2019. When aggregated documents are wide, showed
#  agg[,1] <- droplevels(agg[,1])				
 dtmagg$nrow <- length(levels(agg[,1]))				
 dtmagg$i <- as.numeric(agg[,1])				
 dtmagg$j <- agg[,2]				
 dtmagg$v <- agg[,3]				
 dtmagg$dimnames$Docs <- levels(agg[,1])				
 detOccAgg <- occurrFunc(dtmagg, "after",detOccAgg, TRUE) 				
 detOcc <- detOccAgg			

 # Changed 10/Jan/2019. When aggregated documents are wide, showed
 if(length(levels(droplevels(agg[,1])))!= dtmagg$nrow){
   agg[,1] <- droplevels(agg[,1])				
   dtmagg$nrow <- length(levels(agg[,1]))				
   dtmagg$i <- as.numeric(agg[,1])				
   dtmagg$j <- agg[,2]				
   dtmagg$v <- agg[,3]				
   dtmagg$dimnames$Docs <- levels(agg[,1])				
   }
 
 
  } else {  detOcc <- occurrFunc(dtm, "after",detOcc, FALSE)}				

#--------- If there is aggregation with supplementary variables			
if (!is.null(var.agg)) {  			
 qualivar <- NULL; qualitable <- NULL; qualincat <-NULL			
 T2 <- as.matrix(dtm)			
 if(!is.null(context.quali))			
 for (i in 1:length(context.quali)) {			
   dis.X <- tab.disjonctif(base[,context.quali[i]])			
   T1 <- t(T2)%*% dis.X			
   sumcateg <- which(colSums(T1)==0)			
  if (length(sumcateg)>0) T1 <- T1[,-sumcateg]			
   acpos <- ncol(T1)			
   qualitable <- cbind(qualitable, T1)			
   qualivar <- rbind(qualivar,context.quali[i])			
   qualincat <- rbind(qualincat,acpos)			
}			
# Numerical variables		
  nnum <- ncol(data.context.quanti)		
  quantivar <- NULL; qcolname <- NULL		
 if(!is.null(nnum)) {		
 for (i in 1:nnum) {		
  qcolname <- c(qcolname,colnames(data.context.quanti[i]))	
     if(any(is.na(data.context.quanti[,i])))
       warning("\n", colnames(data.context.quanti[i]), 
         " variable has missing values.\n They will be replaced by the mean of the category\n
         \n Consider to use missMDA R package") 
  qcateg <- aggregate(data.context.quanti[,i], by=list(base[,var.agg]), FUN=mean, na.rm=TRUE)		
  acpos <- which(qcateg[,1]%in% dtmagg$dimnames$Docs)		
  qcateg <- qcateg[acpos,]		
  acpos <- qcateg[,1]		
  qcateg <- data.frame(qcateg[,-1])		
  rownames(qcateg) <- acpos ; colnames(qcateg) <- colnames(data.context.quanti[i])		
  quantivar <- cbind(quantivar,qcateg[,1])		
  rownames(quantivar) <- acpos 		
  colnames(quantivar) <- qcolname		
}}		
 dtm <- dtmagg 		
} # Final aggregation		

# -------------  Computing final frequencies and tables				
 Nfreqword<-tapply(dtm$v,dtm$j,sum)				
 Ndocword<-tapply(dtm$v>0,dtm$j,sum)				
 Table <- cbind(Nfreqword,Ndocword)				
 rownames(Table) <- dtm$dimnames$Terms				
 colnames(Table) <- c("Frequency", "N.Documents")				
 TFreq <- Table[order(Nfreqword, Ndocword, decreasing = TRUE), ]					


 if(is.null(var.agg)) {	
  dsel <- unique(dtm$i)
  dtm$nrow <- length(dsel)	
  rownamesdocs.no.empty <- rownames(base)[dsel]
  dsel1 <- which(dtm$dimnames$Docs %in% rownamesdocs.no.empty) 
  dtm$dimnames$Docs <- dtm$dimnames$Docs[dsel1]			
  dtm$i<- as.numeric(factor(dtm$i,labels=c(1:length(dtm$dimnames$Docs))))
if(segment==TRUE) {
  tcoln <-   colnames(tab.seg)
  tab.seg <- data.frame(tab.seg[dsel,,drop=FALSE])
  rownames(tab.seg) <- dtm$dimnames$Docs
  colnames(tab.seg) <- tcoln
    }
} 

 
 
 			
if(!is.null(var.agg)) { 
# Remove words in SourceTerm supressed in DocTerm
  wordsafteragg <- dtm$dimnames$Terms
  sel.words <- which(SourceTerm$dimnames$Terms%in%wordsafteragg)	
  if(length(sel.words)>0) SourceTerm <- selectFunc(SourceTerm,sel.words)	
  SourceTerm$j <- as.integer(SourceTerm$j)
  SourceTerm$i <- as.integer(SourceTerm$i)
  sumwdocs <- slam::row_sums(SourceTerm)
  pos.noempty <- which(sumwdocs>0)
  SourceTerm$nrow <- length(pos.noempty)
  SourceTerm$dimnames$Docs <- SourceTerm$dimnames$Docs[pos.noempty]	
  dsel <- which(SourceTerm$dimnames$Docs %in% rownamesdocs.no.empty) 		
  SourceTerm$dimnames$Docs <- SourceTerm$dimnames$Docs[dsel]		
  SourceTerm$i<- as.numeric(factor(SourceTerm$i,labels=c(1:length(SourceTerm$dimnames$Docs))))		
  dfvaragg <- dfvaragg[rownamesdocs.no.empty,,drop=FALSE]
  SourceTerm.freq <- Docs.before[rownames(SourceTerm),]
}

 
if(!is.null(var.agg)) { 	
if(!is.null(qualincat)){	
 qualincat <- data.frame(qualincat, row.names=NULL)	
 qualivar <- data.frame(qualivar)
 rownames(qualincat) <- rownames(qualivar)		
 qualivar <- cbind(qualivar, qualincat)
 coltmp <- colnames(qualitable) 
 qualitable <- t(qualitable)
  rownames(qualitable) <- coltmp	
} else {qualitable <- NULL; qualivar <- NULL}
 quali <- list(qualitable=qualitable, qualivar=qualivar)	
 context <- list(quali=quali ,quanti=quantivar)		
} else {	
  context <- list(quali=data.frame(base[,context.quali,drop=FALSE]), quanti=data.context.quanti)		
}		

 
 
 

#--------- Compute results for the total of docs  ------------------				
seqDoc <- c(N, sum(detOcc[,2]), wordsafter, round(sum(detOcc[,2])/N,2))				
rwnDoc <- c("Documents","Occurrences","Words","Mean-length")				
detOcc$PctLength.after <- 100*detOcc[,"Occurrences.after"]/sum(detOcc[,"Occurrences.after"])				

if(is.null(var.agg)){				
  detOcc$MeanLength100.after <- round(N*100*detOcc[,"Occurrences.after"]/sum(detOcc[,"Occurrences.after"]),2)				
} else {				
 detOcc$MeanLength100.after <- round((100*detOcc[,"Occurrences.after"]/detOcc[,"NumberDocs"])				
      /(sum(detOcc[,"Occurrences.after"])/sum(detOcc[,"NumberDocs"])),2)				
}				
 detOcc$PctLength.after <- round(detOcc$PctLength.after,2)				
				
 
 
 
 
 
# ------------------------   Print summary for Tfreqdoc	
seqDocAf <- c(dtm$nrow, sum(detOcc[,"Occurrences.after"]), dtm$ncol, round(sum(detOcc[,"Occurrences.after"])
   /dtm$nrow,2))				
ndocAft <- length(unique(dtm$i))
ndocFEmpty <- ndocIni - ndocAft		
if(ndocFEmpty>0) {		
rwnDoc <-c(rwnDoc, "NonEmpty.Docs", "NonEmpty.Mean-length")		
seqDoc <- c(seqDoc, (ndocIni-ndocIniEmpty),round(sum(detOcc[,"Occurrences.before"])/(ndocIni-ndocIniEmpty),2))		
seqDocAf <- c(seqDocAf, ndocAft,round(sum(detOcc[,"Occurrences.after"])/(ndocAft),2))		
}
		
seqDoc <- c(seqDoc, seqDocAf) 		
 mTfreqdoc <- matrix(seqDoc, ncol=2, byrow=FALSE)		
 rownames(mTfreqdoc) <- rwnDoc		
 colnames(mTfreqdoc) <- c("Before", "After")	
 info <- infoNew()		
 

 attr(dtm, "language") <- info$idiom[[1]] # Version 1.3.1
y <- list(summGen=mTfreqdoc,summDoc=detOcc, indexW = TFreq, DocTerm =dtm) 

if(segment==TRUE) {
  y$indexS <- Index.segments
  y$DocSeg <- slam::as.simple_triplet_matrix(tab.seg)
}

y$info <- info

if(!is.null(var.agg)) { 

  
  
 y$SourceTerm <- SourceTerm 
 y$SourceTerm.freq <- SourceTerm.freq
if(!is.null(context.quali)) y$SourceTerm.qual <- base[,context.quali]
 

 y$var.agg <- var.agg.seg[rownames(y$SourceTerm),,drop=FALSE]
 if(!is.null(context$quali)) y$context$quali <- context$quali
 if(!is.null(context$quanti)) y$context$quanti <- context$quanti
}
 y$remov.docs <- remov.docs

if(is.null(var.agg)) { 
 if(!is.null(context$quanti))  {
   if(length(remov.docs)==0) {y$context$quanti <- context$quanti}
   else {
  pos <- which(rownames(context$quanti) %in% remov.docs)
  y$context$quanti <- context$quanti[-pos,,drop=FALSE] 
        } 
 } # Final !is.null(context$quanti)

 if(!is.null(context$quali))  {
   if(length(remov.docs)==0) {y$context$quali <- context$quali}
   else {
  pos <- which(rownames(context$quali) %in% remov.docs)
  y$context$quali <- context$quali[-pos,,drop=FALSE] 
   } 
   i <- sapply(y$context$quali, is.character)
   if(length(i)>0) y$context$quali[i] <- lapply(y$context$quali[i], as.factor) # version 1.3.1
   
 } # Final !is.null(context$quali)  
} # Final if(is.null(var.agg))


 df <- y$summDoc[,2,drop=FALSE]
 rownames(df) <- y$summDoc[,1]
 y$rowINIT <- df
 
 
  class(y) <- c("TextData", "list")
if(blongErr==TRUE) warning("Only repeated segments < 20 words have been computed")	

if(graph==TRUE) plotTextData() 

return(y)
}
