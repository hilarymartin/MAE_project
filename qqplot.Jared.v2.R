cols <- paste0(c("#E41A1C","#377EB8","#4DAF4A","#984EA3"),"CC")
argv <- commandArgs(TRUE)

thr <- 5e5 #recombinations within this distance of one another are filtered
male.lengths <- c(1.9512,1.8955,1.6071,1.4654,1.512,1.3762,1.2835,1.0794,1.1725,1.3389,1.0936,1.3554,1.0131,0.9462,1.0257,1.081,1.0856,0.9862,0.9264,0.7472,0.4731,0.4896)
female.lengths <- c(3.4541,3.2541,2.7564,2.5906,2.6019,2.4159,2.3033,2.0994,1.982,2.1813,1.9553,2.0664,1.5588,1.4236,1.5496,1.4962,1.6153,1.4257,1.2682,1.2197,0.764,0.8276)
malesize <- sum(male.lengths)
femalesize <- sum(female.lengths)

qqr <- function(k1,rate,col1,makeplot=TRUE,...) {
  k1 <- sort(k1)
  n <- length(k1)
  e <- qpois(1:n/(n+1),rate)
  l <- c(min(k1),min(max(k1),rate+6*sqrt(rate)))
  etmp <- e[k1>=l[2]]
  ret <- cbind(jitter(etmp),l[2])

  if(makeplot) {
      plot(e,k1,type='n',xlab="",ylab="",ylim=l,...)
      abline(0,1,col=1,lty=2);grid()
      points(jitter(e[k1<l[2]]),jitter((k1[k1<l[2]])),col=col1,pch=16)
      points(ret,col=col1,pch=17)
  }

  return(  cbind(e,k1))
}

tabulateRecombination <- function(r) {
    rec <- data.frame(r[!duplicated(r[,c('CHILD','PARENT')]),1:2],stringsAsFactors=FALSE)
    rec$duo <- paste(rec$CHILD,rec$PARENT,sep="-")
    rec$sex <- as.factor(c("Male","Female")[as.integer(fam[match(rec$PARENT,fam$V2),5])] )
    rec$informative <- fam[match(rec$PARENT,fam$V2),]$grandparent | fam[match(rec$PARENT,fam$V2),]$nkid>2
    rec$nkid <-  fam[match(rec$PARENT,fam$V2),]$nkid
    rec$nrec <- table(r$duo)[rec$duo]
    rec <- subset(rec,informative)
    print(tapply(rec$nrec,rec$sex,mean))
    rec
}


fam <- read.table(paste0(argv[1],".fam"),colClasses="character")
fam$nkid <- sapply(fam$V2,function(x) sum(x==fam$V3)+sum(x==fam$V4))
fam$grandparent <- fam$V3 %in% fam$V2 | fam$V4 %in% fam$V2

rall <- list()
for(i in 1:22) {
    rall[[i]] <- read.table(paste0(argv[1],"-chr",i,"-recombinations.txt"),header=TRUE,as.is=TRUE,colClasses=c("character","character","integer","integer","numeric"))
#    rall[[i]] <- read.table(paste0(argv[1],"-chr",i,"-recombinations.txt"),header=FALSE,as.is=TRUE,colClasses=c("character","character","integer","integer","numeric"),col.names=c("CHILD","PARENT","START","END","PROB_RECOMBINATION"))
    rall[[i]]$chr <- i
}

r <- subset(do.call("rbind",rall),PROB_RECOMBINATION>0.5)
r$duo <- paste0(r$CHILD,"-",r$PARENT)
r$sex <- as.factor(c("Male","Female")[as.integer(fam[match(r$PARENT,fam$V2),5])] )
r$informative <- fam[match(r$PARENT,fam$V2),]$grandparent | fam[match(r$PARENT,fam$V2),]$nkid>2
r$diff <- c(NA,tail(r$END,-1) - head(r$START,-1))
r$diff[c(FALSE,tail(r$chr,-1) != head(r$chr,-1)) | c(FALSE,tail(r$CHILD,-1) != head(r$CHILD,-1)) | c(FALSE,tail(r$PARENT,-1) != head(r$PARENT,-1))] <- NA
r.flt <- r[-c(which(r$diff<thr),which(r$diff<thr)-1),]

dup.pos <- r[which(duplicated(r[,c("START","END")])),c("START","END")]
r.flt2 <- r[-which(r$START%in%dup.pos[,1] & r$END%in%dup.pos[,2] ),]

#print(head(r))
stopifnot(nrow(r)>0)
###takes the proper expectation - this is the right thing to do *if* our probs are well calibrated.
## r <- subset(do.call("rbind",rall))
## r$pair <- paste0(r$CHILD,"-",r$PARENT)
#x <- tapply(r$PROB_RECOMBINATION,r$pair,sum)

## rec <- data.frame(r[!duplicated(r[,c('CHILD','PARENT')]),1:2],stringsAsFactors=FALSE)
## rec$duo <- paste(rec$CHILD,rec$PARENT,sep="-")
## rec$sex <- as.factor(c("Male","Female")[as.integer(fam[match(rec$PARENT,fam$V2),5])] )
## rec$informative <- fam[match(rec$PARENT,fam$V2),]$grandparent | fam[match(rec$PARENT,fam$V2),]$nkid>2
## rec$nkid <-  fam[match(rec$PARENT,fam$V2),]$nkid
## rec$nrec <- x[rec$duo]
## rec <- subset(rec,informative)

rec <- tabulateRecombination(r)
rec.flt <- tabulateRecombination(r.flt)
rec.flt2 <- tabulateRecombination(r.flt2)

rec$nrec_filtered <- rec.flt$nrec[match(rec$duo,rec.flt$duo)]
rec$nrec_filtered2 <- rec.flt2$nrec[match(rec$duo,rec.flt2$duo)]


#print(tapply(rec$nrec,rec$sex,mean))
stopifnot(nrow(rec)>0)


print("QQ-PLOT")
pdf(paste0(argv[1],"-qq.pdf"),width=12,height=6)

clook=c("paternal"=cols[3],"maternal"=cols[4])
par(mar=c(2,2,1,1),mfrow=c(1,2),cex=1.6,oma=c(2,2,0,0))
tmp <- qqr(subset(rec,sex=="Male")$nrec,malesize,cols[1])
points(qqr(subset(rec.flt,sex=="Male")$nrec,malesize,makeplot=0),pch=16,col=cols[2])
points(qqr(subset(rec.flt2,sex=="Male")$nrec,malesize,makeplot=0),pch=16,col=cols[3])
legend("topleft",c("Unfiltered",paste0("<",thr/1000,"kb filtered"),"Duplicate positions removed"),col=cols[1:3],pch=16,bg="white")
tmp <- qqr(subset(rec,sex=="Female")$nrec,femalesize,cols[1])
points(qqr(subset(rec.flt,sex=="Female")$nrec,femalesize,makeplot=0),pch=16,col=cols[2])
points(qqr(subset(rec.flt2,sex=="Female")$nrec,femalesize,makeplot=0),pch=16,col=cols[3])

mtext("Expected recombinations genome-wide",1,outer=T,cex=1.6,padj=1.7)
mtext("Observed recombinations genome-wide",2,outer=T,cex=1.6,padj=-1)

dev.off()

#print(head(rec))
print(paste0(argv[1],".rec"))
rec <- rec[order(rec$sex,rec$nrec),]
write.table(x=rec,file=paste0(argv[1],".rec"),quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")


r <- subset(r,informative)
nmale   <- sum(r[!duplicated(r[,1:2]),"sex"]=="Male")
nfemale <- sum(r[!duplicated(r[,1:2]),"sex"]=="Female")

pdf(paste0(argv[1],"-perchromosome.pdf"),width=12,height=6)

par(mar=c(2,2,1,1),mfrow=c(1,2),cex=1.6,oma=c(2,2,0,0))
xhat <- seq(0,10,.01)
y <- table(subset(r,sex=="Male")$chr)/nmale
l <- c(0,max(male.lengths,y))
plot(male.lengths,y,ylab="Mean per chromosome",xlab="Expected per chromosome",type='n',xlim=l,ylim=l)
mtext(paste("Paternal (n =",nmale,")"),3,cex=1.6)
y1 <- xhat - 1.96*sqrt(xhat)/sqrt(nmale)
y2 <- xhat + 1.96*sqrt(xhat)/sqrt(nmale)
polygon(c(xhat,rev(xhat)),c(y1,rev(y2)),col='light grey',border=NA)
abline(0,1,col=1,lty=2)
text(male.lengths,y,labels=paste(1:22),col=cols[1])

y <- table(subset(r,sex=="Female")$chr)/nfemale
l <- c(0,max(female.lengths,y))
plot(female.lengths,y,ylab="Mean per chromosome",xlab="Expected per chromosome",type='n',xlim=l,ylim=l)
mtext(paste("Maternal (n =",nfemale,")"),3,cex=1.6)
y1 <- xhat - 1.96*sqrt(xhat)/sqrt(nfemale)
y2 <- xhat + 1.96*sqrt(xhat)/sqrt(nfemale)
polygon(c(xhat,rev(xhat)),c(y1,rev(y2)),col='light grey',border=NA)
abline(0,1,col=1,lty=2)
text(female.lengths,y,labels=paste(1:22),col=cols[1])

mtext("Expected recombinations per chromosome",1,outer=T,cex=1.6,padj=1.7)
mtext("Observed recombinations per chromosome",2,outer=T,cex=1.6,padj=-1)

#legend("topleft",c("SHAPEIT2","Merlin"),col=cols[1:2],pch=16,bg='white')
dev.off()

save.image(paste0(argv[1],".RData"))

