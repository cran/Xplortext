#' @import ggplot2 
#' @importFrom grDevices dev.new dev.interactive
#' @importFrom utils head write.table
#' @importFrom plotly plotly_build
#' @export
plot.TextData <- function (x, ndoc=25, nword=25, nseg=25, sel=NULL, ordFreq=TRUE,
  stop.word.tm=FALSE, idiom= "en", stop.word.user=NULL, theme=theme_bw(), title=NULL, 
  xtitle=NULL, col.fill="grey", col.lines="black", text.size=12, freq=NULL, vline=NULL, 
  interact=FALSE, round.dec=4,...)
{

  if(!is.null(vline))
  if(vline=="median")
  {
    warning("median is deprecated, it is changed by the mean")
    vline <- "mean"
  }

  options(stringsAsFactors = FALSE)
if (!inherits(x, "TextData")) 
 stop("Object x should be TextData class")
  pdoc <- pword <- pseg <- NULL 
 if(!is.null(vline)) if(vline==FALSE) vline <-NULL
  
idiom <- x$info$idiom[[1]]
theme$text$size <- text.size
words<-iFreq<-docnames<-rsegment<-frequency<-DocName<-y.var<-NULL

if(!is.null(freq)) if(freq=="YES"| freq==TRUE) freq <- 5


###############################################################################################################
## sel. will have the selected words/docs/segments, sel.tot. all the possible elements
#  1) bsel.seg <- bsel.word <- bsel.doc <- FALSE
bsel.word <- bsel.doc <- bsel.seg <- FALSE

n.Init.word <-length(x$DocTerm$dimnames$Terms)
n.Init.doc <- length(x$DocTerm$dimnames$Docs)


# 2) if("indexS" %in% names(x)) bSegment <- TRUE else bSegment <- FALSE
# 3) if(bSegment) n.Init.seg <- length(x$indexS$segOrderlist$segment)

if("indexS" %in% names(x)) bSegment <- TRUE else bSegment <- FALSE
if(bSegment) n.Init.seg <- length(x$indexS$segOrderlist$segment)


########################################################################################

if(is.null(sel)) {
  if(interact==TRUE) stop("seg=NULL is not allowed with interact=TRUE. \n Select word, doc or seg in sel argumment")
  bsel.seg <- bsel.word <- bsel.doc <- TRUE
  
  if(bSegment) s.init.seg <- s.sel.seg <- x$indexS$segOrderlist$segment
  s.init.word <- s.sel.word <- x$DocTerm$dimnames$Terms
  s.init.doc <- s.sel.doc <- x$DocTerm$dimnames$Docs

} else { # No null sel
  
  if(length(sel)==1){
    
    if(bSegment) if(sel=="seg") {ndoc<-0; nword<-0;  s.sel.seg <- s.init.seg <- x$indexS$segOrderlist$segment
    n.Init.seg <- length(s.init.seg); bseg <- TRUE
    }

    if(sel=="word") {ndoc<-0; nseg<-0;  s.sel.word <- s.init.word <- x$DocTerm$dimnames$Terms
    n.Init.word <- length(s.init.word); bword <- TRUE
    }
    
    if(sel=="doc") {nword<-0; nseg<-0;  s.sel.doc <- s.init.doc <- x$DocTerm$dimnames$Docs
    n.Init.word <- length(s.init.doc); bdoc <- TRUE
    } 
    
    
  } else {   # no length(sel)==1
    
    if(bSegment) if(sel$type=="seg") { 
      ndoc <- nword <- 0
      s.init.seg <- x$indexS$segOrderlist$segment
      
      if(is.numeric(sel$select)) s.sel.seg <- s.init.seg[sel$select]
      else {
        if(length(sel$select)==1) {
          s.sel.seg<- s.init.seg[grepl(sel$select, s.init.seg)]
        } else
          s.sel.seg <- s.init.seg[s.init.seg %in% sel$select] 
      }
      
    }
    
    if(sel$type=="word"){ 
      ndoc <- nseg <- 0
      s.init.word <- x$DocTerm$dimnames$Terms

      if(is.numeric(sel$select)) s.sel.word <- s.init.word[sel$select]
      else  s.sel.word <- s.init.word[s.init.word %in% sel$select] }
    
    if(sel$type=="doc"){ 
      nword <- nseg <- 0
      s.init.doc <- x$DocTerm$dimnames$Docs
      if(is.numeric(sel$select)) s.sel.doc <- s.init.doc[sel$select]
      else  s.sel.doc <- s.init.doc[s.init.doc %in% sel$select] }
  } # End   length(sel)
} # End sel
#################################



########################  Segments   ####################################
if(bSegment) if(nseg>0){
  df <- data.frame(x$indexS$segOrderlist)  # All the segments from TextData
  colnames(df)[1] <- "rsegment"
  #  nseg.TD <- nrow(df)
  
  df.rs <- df
  df.rs$rank <- rank(-df.rs$frequency, ties.method = "average")
  df.rs <- df.rs[df.rs$rsegment %in%  s.sel.seg,]
  df.rs <- df.rs[with(df.rs, order(-frequency,rsegment)), ]
  
  bsel <- ifelse(nrow(df.rs) == n.Init.seg, FALSE, TRUE)  # For titles
  nseg <- min(nseg,nrow(df.rs))
  # Order by freq
  

  if(nseg<nrow(df.rs))  df.rs <- df.rs[c(1:nseg),]
  
  if(is.null(title)){
    if(bsel) ttitle.s <- paste0( nseg," selected segments" )
    else ttitle.s <- paste0( nseg," most frequent segments" )
  } else  ttitle.s <- title
  
  ifelse(is.null(xtitle), txtitle.s <- "Repeated segment frequency" ,txtitle.s <- xtitle)
  
  #  df <- merge(df, df.rs, by="rsegment")
  #  df <- df[, c(1,2,6)]
  colnames(df.rs)[2] <- "iFreq"
  
  df.rs <- df.rs[order(-df.rs$iFreq,df.rs$rsegment), ]
  

  if(ordFreq) level_order <- rev(df.rs$rsegment)
  else level_order <-  df[rev(order(df$rsegment)),]$rsegment
  
  pseg <- ggplot(df.rs)+ geom_bar(aes(x=factor(rsegment, levels=level_order), y=iFreq),stat = "identity", color = col.lines, fill = col.fill) + 
    coord_flip() +
    ylab(txtitle.s) + xlab("") + ggtitle(ttitle.s) + theme + theme(plot.title = element_text(hjust = 0.5))
  
  if(!is.null(freq))
    pseg <- pseg + geom_text(aes(x=df.rs$rsegment, y=iFreq+ freq, label = iFreq))
  
  if(!is.null(vline)) if(is.numeric(vline)) 
    pseg<- pseg +  geom_hline(yintercept=vline, linetype="dashed", color = "red")
  

  if(interact==TRUE) {
    my.text <-  paste("Rank", pseg$data$rank, "", pseg$data$rsegment, "<br>", "Freq: ", pseg$data$iFreq, "<br>")
    g <- plotly::plotly_build(pseg)
    g$x$data[[1]]$text <- my.text
    pseg <- g
  }
}


##################################################  End Segments ####################################



##############################  stop words tm and users  #######################################

if(ndoc>0 | nword > 0) {
  
  # It is not a segment (it is word or doc)
  df <- data.frame(x$indexW[,1], rownames(x$indexW))
  colnames(df) <- c("iFreq", "words") 
  # nword.TD <- nrow(df)
  
  
  # Remove stopwords
  if(stop.word.tm!=FALSE) {
    stop.word <- tm::stopwords(idiom)
    pos <- which(!df$words%in%stop.word)
    df <-  df[pos,]
  }
  
  # Remove user words
  if(!is.null(stop.word.user)) {
    pos <- which(!df$words %in% stop.word.user)
    df <-  df[pos,]}
  
  final.words <- rownames(df)
  #  nword <- min(nword,nrow(df))
  
}
################################  End stop words tm and users ######################################


if(nword>0){
  ##################################
  # Selection of words by number and name
  df.rw <- df
  df.rw$rank <- rank(-df.rw$iFreq, ties.method = "average")
  df.rw <- df.rw[df.rw$words %in%  s.sel.word,]
  df.rw <- df.rw[with(df.rw, order(-iFreq,words)), ]

  
  # Mean words before
  vlineword.Init.w <- x$summGen[2,2] / x$summGen[3,2]
  
  b.sel.word <- ifelse(nrow(df.rw) == n.Init.word, FALSE, TRUE)  # For titles
  
  nword <- min(nword,nrow(df.rw))
  # Order by freq
  
  if(nword<nrow(df.rw))  df.rw <- df.rw[c(1:nword),]
  

  if(is.null(title)){
    if(b.sel.word) ttitle.w <- paste0(nword," selected words" )
    else ttitle.w <- paste0(nword," most frequent words" )
  } else  ttitle.w <- title
  
  
  ifelse(is.null(xtitle), txtitle.w <- "Word frequency" ,txtitle.w <- xtitle)

  
  
  # Mean  for indexW
  if(!is.null(vline)) {
    if(is.numeric(vline)) {vlineword.Final.w <-  vlineword.Init.w <- vline}
    if(vline=="YES"|vline==TRUE|vline=="mean" ) {
      vlineword.Final.w <- mean(df$iFreq)
    } }
  
  df.rw <- df.rw[with(df.rw, order(-iFreq,words)), ]
  
  
  if(ordFreq) level_order.w <- rev(df.rw$words)
  else level_order.w <-  df.rw[rev(order(df.rw$words)),]$words
  

  pword <- ggplot(df.rw)+ geom_col(aes(x=factor(words, levels=level_order.w), y=iFreq),color = col.lines, fill = col.fill) + 
     coord_flip() +
    ylab(txtitle.w) + xlab("") + ggtitle(ttitle.w) + theme + theme(plot.title = element_text(hjust = 0.5))

  if(!is.null(freq))
    pword <- pword + geom_text(aes(x=df.rw$words, y=iFreq+ freq, label = iFreq))
  

  if(!is.null(vline)) {
    pword <- pword+ geom_hline(yintercept=vlineword.Final.w, linetype="dashed", color = "blue")
    if(vlineword.Final.w != vlineword.Init.w)
     pword <- pword+ geom_hline(yintercept=vlineword.Init.w, linetype="dashed", color = "red")

  }

  


    if(interact==TRUE) {
   my.text <-  paste("Rank", pword$data$rank, "Word: ", pword$data$words, "<br>", "Freq: ", pword$data$iFreq, "<br>",
                      "%Before " , round(100*pword$data$iFreq/ x$summGen[2,1],round.dec), "<br>",
                      "%After " , round(100*pword$data$iFreq/ x$summGen[2,2],round.dec))
    

    suppressWarnings(g <- plotly::plotly_build(pword))
    g$x$data[[1]]$text <- my.text
    pword <- suppressWarnings(g) 
  }
} # End Words



####################################################   Docs   #################################################

if(ndoc>0){
  # If aggegate documents by qualitative variable mean length is shown. In summary and print both mean and total are shown
  var.agg <- ifelse(is.null(x$var.agg), FALSE, TRUE)
  
  mean.before.d <- x$summGen[2,2]/x$summGen[1,1]

  if(var.agg==FALSE) {
    df.doc <- as.matrix(x$DocTerm)[,df$words]
    df.doc <- data.frame("DocName"=rownames(df.doc), "iFreq"=rowSums(df.doc))
    mean.after.d <- sum(df.doc$iFreq)/ x$summGen[1,1]
   
    a <- suppressWarnings(sum(as.numeric(df.doc$DocName)))  # If a is NA, non numeric rownames
    if(!is.na(a)) df.doc$DocName <- as.numeric(df.doc$DocName)
    df.doc <- df.doc[with(df.doc, order(-iFreq, DocName)), ]
    colnames(df.doc) <- c("DocName", "Number.Words")
    # Ranks
    df.rd <- df.doc
    df.rd$rank <- rank(-df.rd$Number.Words, ties.method = "average")
    df.rd <- df.rd[df.rd$DocName %in% s.sel.doc,]
    mean.after.d <- sum(df.rd$Number.Words) / nrow(df.rd)
  }
  

  
   if(var.agg) {
    df.before <- x$summDoc[, c("DocName", "Occurrences.after", "NumberDocs")]
    df.after <- as.matrix(x$DocTerm[,final.words])
    df.doc <- data.frame("DocName"=rownames(df.after), "iFreq"=rowSums(df.after))
    df.doc <- merge(df.doc, df.before, by=c("DocName"))
    df.doc$MeanWords <- df.doc$iFreq / df.doc$NumberDocs
    df.doc <- df.doc[, c("DocName",  "iFreq", "NumberDocs", "MeanWords")]
    colnames(df.doc) <- c("DocName",  "Number.Words", "Number.Docs", "Mean.Words")
    df.rd <- df.doc
    
    df.rd<- df.rd[order(-df.rd$Mean.Words, df.rd$DocName),]
    
    df.rd$rank <- rank(-df.rd$Mean.Words, ties.method = "average")
    df.rd <- df.rd[df.rd$DocName %in% s.sel.doc,]
    df.rd$Mean.Words <-  as.numeric(df.rd$Mean.Words)
    df.rd$Number.Docs <- as.integer(df.rd$Number.Docs)
    mean.after.d <- sum(df.rd$Number.Words) / sum(df.rd$Number.Docs)
  }

  ndoc <- min(ndoc,nrow(df.rd))

  
  if(is.null(title))
    if(var.agg){
      ttitle.d<- ifelse(nrow(df.rd)==ndoc, paste0("Mean length of each category"), 
                        paste0("Mean length of ", ndoc, " longest categories"))
    } else {  # no varagg
      ttitle.d<- ifelse(nrow(df.rd)==ndoc, paste0("Mean length of each document"), 
                        paste0("Mean length of ", ndoc, " longest documents"))
    } else ttitle.d <- title
  
  if(is.null(xtitle)) {
    txtitle.d <- ifelse(var.agg, "Mean words of the category", "Number of words") 
  } else txtitle.d <- xtitle
  
  
  df.rd <- df.rd[1:ndoc,]
  
  df.rd$DocName <- as.character(df.rd$DocName)
  
  if(ordFreq==FALSE) {
    if(is.numeric(df.rd$DocName))  level_order <- as.character(sort(df.rd$DocName, decreasing=TRUE))
    else  level_order.d <- sort(df.rd$DocName, decreasing=TRUE)
  } else level_order.d <- rev(df.rd$DocName)
  

  if(var.agg) {
    df.rd$y.var <- df.rd$Mean.Words
    
  } else {
    df.rd$y.var <- df.rd$Number.Words
    df.rd$DocName <- factor(df.rd$DocName, levels = df.rd$DocName[order(df.rd$DocName)])
  }
 
  
   pdoc <- ggplot(df.rd)+ geom_bar(aes(x=factor(DocName, levels=level_order.d), y=y.var),stat = "identity",
                                  color = col.lines, fill = col.fill) +   coord_flip() +
    ggtitle(ttitle.d) + theme + theme(plot.title = element_text(hjust = 0.5)) + ylab(txtitle.d) + xlab("") 


  if(!is.null(freq))
    pdoc <- pdoc + geom_text(aes(x=df.rd$DocName, y=df.rd$y.var + freq, label = round(df.rd$y.var,round.dec)))  
  

  
  if(!is.null(vline)) {
    if(is.numeric(vline)) mean.before.d <- mean.after.d  <- vline
    pdoc <- pdoc + geom_hline(yintercept=mean.before.d, linetype="dashed", color = "red") +
     geom_hline(yintercept=mean.after.d, inetype="dashed", color = "blue")
  } 

  
  if(interact==TRUE) {
    if(var.agg) 
      my.text <-  paste("Cat: ",  pdoc$data$DocName," Number.Docs: ", pdoc$data$Number.Docs,  
                        "<br>", "Rank", pdoc$data$rank,  " Mean.Words: ", round(pdoc$data$Mean.Words,round.dec)  )
    else my.text <-  paste("Doc: ",  pdoc$data$DocName," Words: ", pdoc$data$y.var,  
                           "<br>", "Rank", pdoc$data$rank)
  #  g <- plotly::plotly_build(pdoc)

    g <- plotly::ggplotly(pdoc)
    g$x$data[[1]]$text <- my.text
    pdoc <- g

  }
} # End ndoc


if(!is.null(pseg)) {
  print(pseg)
  if(interact==TRUE) return(pseg) } 

if(!is.null(pword)) {
  suppressWarnings(print(pword))
  if(interact==TRUE) suppressWarnings(return(pword))
  }

if(!is.null(pdoc)) {
  print(pdoc)
  if(interact==TRUE) return(pdoc)}
}
