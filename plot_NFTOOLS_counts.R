setwd("/gpfs1/well/donnelly/hilary/maternal_age_and_recombination")

mydir="/home/hilary/maternal_age_recombination/NFTOOLS_results"

QTR.data=read.delim("/well/donnelly/hilary/TwinRecombination/Bayesian_modeling/input_data/all_Brisbane_crossovers.k5.chips_separated.excluding_triplets_and_aneuploidies.txt",header=T,stringsAsFactors=F)
rownames(QTR.data)=QTR.data$child
QTR610K.data=QTR.data[QTR.data$chip =="610K",]
QTR370K.data=QTR.data[QTR.data$chip =="370K",]

#NTR
NTR.data=read.delim("/well/donnelly/hilary/TwinRecombination/Bayesian_modeling/input_data/crossover_counts_HMM_and_NFTOOLS.families_with_min_3_unique_kids.all.txt",header=T,stringsAsFactors=F)
NTR.data$chip="NTR"
NTR.data=NTR.data[,c("GID", "FATH","MOTH","chip", "NFTOOLS_mat" ,"NFTOOLS_pat","mat.age.birth","pat.age.birth")]
colnames(NTR.data)=colnames(QTR.data)
rownames(NTR.data)=NTR.data$child
#FC
fc.data=read.delim("/well/donnelly/hilary/TwinRecombination/Bayesian_modeling/input_data/FC_NFTOOLS_counts.k5.txt",header=T,stringsAsFactors=F)
fc.data$chip="FC"
fc.data$child=rownames(fc.data)
fc.data=fc.data[,c("child","father","mother","chip","mat.count","pat.count","mat.age","pat.age")]
colnames(fc.data)=colnames(QTR.data)

#fc.duoHMM=read.delim("FC/duoHMM_results/FC.generr_removed.total_recombination_counts.maternal.p_0.5.txt",header=T,stringsAsFactors=F)
#table(fc.duoHMM[fc.duoHMM$child %in% fc.data$child,"total_same_mother"])


rownames(fc.data)=fc.data$child
families=as.integer(unlist(lapply(strsplit(rownames(fc.data),"_"),function(x){return(x[2])})))
rownames(fc.data)=paste(as.character(families),"-",rownames(fc.data),sep="")
fc.data$child=paste(as.character(families),"-",fc.data$child,sep="")
fc.data$mother=paste(as.character(families),"-",fc.data$mother,sep="")
fc.data$father=paste(as.character(families),"-",fc.data$father,sep="")

#Hutterites
hutt.pat.data=read.delim("/well/donnelly/hilary/TwinRecombination/Bayesian_modeling/input_data/Hutterites_patageeffectTable_282.txt",header=T,stringsAsFactors=F)[,1:4]
hutt.mat.data=read.delim("/well/donnelly/hilary/TwinRecombination/Bayesian_modeling/input_data/Hutterites_matageeffectTable_282.txt",header=T,stringsAsFactors=F)[,1:4]
colnames(hutt.mat.data)=c("family","child","mat.count","mat.age")
colnames(hutt.pat.data)=c("family","pat.age","pat.count","child")
rownames(hutt.mat.data)=hutt.mat.data$child
rownames(hutt.pat.data)=hutt.pat.data$child

hutt.data=data.frame(hutt.mat.data$child,paste0(hutt.mat.data$child,"_mother"),paste0(hutt.mat.data$child,"_father"),rep("HUTT",nrow(hutt.mat.data)),hutt.mat.data$mat.count,hutt.pat.data$pat.count,hutt.mat.data$mat.age,hutt.pat.data$pat.age)
colnames(hutt.data)=colnames(fc.data)


library(scales)
library(ggplot2)

plot.NFTOOLS.data=function(NFTOOLS,mydir,prefix){
NFTOOLS.mat=data.frame(NFTOOLS[,c("child","mat.count","mat.age")])
NFTOOLS.pat=data.frame(NFTOOLS[,c("child","pat.count","pat.age")])

##QQ plots
random.mat=rpois(n=nrow(NFTOOLS.mat),lambda=42.81)
random.pat=rpois(n=nrow(NFTOOLS.pat),lambda=25.9)

pdf(paste0(mydir,"/","qqplots_for_NFTOOLS_k_5_crossovers.",prefix,".pdf"),height=5,width=10)
par(mfrow=c(1,2))
plot(sort(random.mat),sort(NFTOOLS.mat[,"mat.count"]),ylab="Observed crossovers",xlab="Expected crossovers",main=paste0("maternal crossovers ",prefix),xlim=range(c(random.mat,NFTOOLS.mat$mat.count))
,ylim=range(c(random.mat,NFTOOLS.mat$mat.count)),col=alpha("red",0.4),bg=alpha("red",0.4),pch=21)
abline(a=0,b=1,col="black")

plot(sort(random.pat),sort(NFTOOLS.pat[,"pat.count"]),ylab="Observed crossovers",xlab="Expected crossovers",main=paste0("paternal crossovers ",prefix),xlim=range(c(random.pat,NFTOOLS.pat$pat.count))
,ylim=range(c(random.pat,NFTOOLS.pat$pat.count)),col=alpha("blue",0.4),bg=alpha("blue",0.4),pch=21)
abline(a=0,b=1,col="black")
dev.off()

### Histogram
random.mat=rpois(n=10000,lambda=42.81)
random.pat=rpois(n=10000,lambda=25.9)
pdf(paste(mydir,"/","distribution_of_NFTOOLS_k_5_maternal_crossovers.",prefix,".pdf",sep=""),height=5,width=5)
hist(NFTOOLS.mat[,"mat.count"],main=paste("maternal crossovers ",prefix,sep=""),xlab="NFTOOLS count",prob=T,col="grey",ylim=c(0,0.06))
lines(density(random.mat),col="darkgreen")
abline(v=42.81,col="darkgreen",lwd=2)
abline(v=mean(NFTOOLS.mat$mat.count),col="black",lwd=2)
dev.off()

pdf(paste(mydir,"/","distribution_of_NFTOOLS_k_5_paternal_crossovers.",prefix,".pdf",sep=""),height=5,width=5)
hist(NFTOOLS.pat[,"pat.count"],main=paste("paternal crossovers ",prefix,sep=""),xlab="NFTOOLS count",prob=T,col="grey")
lines(density(random.pat))
abline(v=25.9,col="darkgreen",lwd=2)
abline(v=mean(NFTOOLS.pat$pat.count),col="black",lwd=2)
dev.off()

### # crossovers vs age
mat.regress = lm(as.numeric(NFTOOLS.mat$mat.count)~NFTOOLS.mat$mat.age)
mat.predict = predict(mat.regress, interval="confidence") 
pdf(paste(mydir,"/","NFTOOLS_k_5_maternal_crossovers_vs_age.",prefix,".pdf",sep=""),height=5,width=5)
plot(as.numeric(NFTOOLS.mat$mat.age),NFTOOLS.mat$mat.count,main=paste("maternal crossovers ",prefix,sep=""),xlab="maternal age",ylab="maternal crossovers",col=alpha("blue",0.4),bg=alpha("blue",0.4),pch=21)
abline(mat.regress,col="red")
#lines(x=NFTOOLS.mat$mat.age,y=mat.predict[,2],lty=2)
#lines(x=NFTOOLS.mat$mat.age,y=mat.predict[,3],lty=2)
legend("topleft",c(paste0("slope = ",round(mat.regress$coefficients[2],digits=2)),paste0("p = ",round(summary(mat.regress)$coefficients[2,4],digits=3))))

dev.off()

}

plot.NFTOOLS.data(QTR610K.data,mydir,"QTR 610K")
plot.NFTOOLS.data(QTR370K.data,mydir,"QTR 370K")
plot.NFTOOLS.data(NTR.data,mydir,"NTR")
plot.NFTOOLS.data(fc.data,mydir,"FC")
plot.NFTOOLS.data(hutt.data,mydir,"Hutterites")


all.datasets=list(QTR610K.data,QTR370K.data,NTR.data,fc.data,hutt.data)
names(all.datasets)=c("QTR 610K","QTR 370K","NTR","FC","Hutterites")

pdf(paste0(mydir,"/","qqplots_for_NFTOOLS_k_5_maternal_crossovers.pdf"),height=3,width=15)
par(mfrow=c(1,5))
for(i in 1:5){
NFTOOLS=all.datasets[[i]]
prefix=names(all.datasets)[i]
NFTOOLS.mat=data.frame(NFTOOLS[,c("child","mat.count","mat.age")])
random.mat=rpois(n=nrow(NFTOOLS.mat),lambda=42.81)
plot(sort(random.mat),sort(NFTOOLS.mat[,"mat.count"]),ylab="Observed crossovers",xlab="Expected crossovers",main=prefix,xlim=range(c(random.mat,NFTOOLS.mat$mat.count))
,ylim=range(c(random.mat,NFTOOLS.mat$mat.count)),col=alpha("red",0.4),bg=alpha("red",0.4),pch=21,cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
abline(a=0,b=1,col="black")
}
dev.off()

pdf(paste0(mydir,"/","qqplots_for_NFTOOLS_k_5_paternal_crossovers.pdf"),height=3,width=15)
par(mfrow=c(1,5))
for(i in 1:5){
NFTOOLS=all.datasets[[i]]
prefix=names(all.datasets)[i]
NFTOOLS.pat=data.frame(NFTOOLS[,c("child","pat.count","pat.age")])
random.pat=rpois(n=nrow(NFTOOLS.pat),lambda=25.9)
plot(sort(random.pat),sort(NFTOOLS.pat[,"pat.count"]),ylab="Observed crossovers",xlab="Expected crossovers",main=prefix,xlim=range(c(random.pat,NFTOOLS.pat$pat.count))
,ylim=range(c(random.pat,NFTOOLS.pat$pat.count)),col=alpha("blue",0.4),bg=alpha("blue",0.4),pch=21,cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
abline(a=0,b=1,col="black")
}
dev.off()

pdf(paste0(mydir,"/","histograms_of_NFTOOLS_k_5_maternal_crossovers.pdf"),height=3,width=15)
par(mfrow=c(1,5))
random.mat=rpois(n=100000,lambda=42.81)
for(i in 1:5){
NFTOOLS=all.datasets[[i]]
prefix=names(all.datasets)[i]
NFTOOLS.mat=data.frame(NFTOOLS[,c("child","mat.count","mat.age")])

hist(NFTOOLS.mat[,"mat.count"],main=prefix,xlab="NFTOOLS count",prob=T,col=alpha("red",0.4),ylim=c(0,0.06),cex.main=1.5,cex.axis=1.5,cex.lab=1.5,xlim=c(19,74))
lines(density(random.mat, n = 2*512),col="darkgreen")
#par(new=F)
#hist(random.mat,col=alpha("darkgreen",0.4),add=T,prob=T)
abline(v=42.81,col="darkgreen",lwd=2)
abline(v=mean(NFTOOLS.mat$mat.count),col="black",lwd=2)
}
dev.off()


pdf(paste0(mydir,"/","histograms_of_NFTOOLS_k_5_paternal_crossovers.pdf"),height=3,width=15)
random.pat=rpois(n=20000,lambda=25.9)
par(mfrow=c(1,5))
for(i in 1:5){
NFTOOLS=all.datasets[[i]]
prefix=names(all.datasets)[i]
NFTOOLS.pat=data.frame(NFTOOLS[,c("child","pat.count","pat.age")])
hist(NFTOOLS.pat[,"pat.count"],main=prefix,xlab="NFTOOLS count",prob=T,col=alpha("blue",0.4),cex.main=1.5,cex.axis=1.5,cex.lab=1.5,xlim=c(12,41),ylim=c(0,0.12))
lines(density(random.pat, n = 512))
abline(v=25.9,col="darkgreen",lwd=2)
abline(v=mean(NFTOOLS.pat$pat.count),col="black",lwd=2)

}
dev.off()


pdf(paste(mydir,"/","NFTOOLS_k_5_maternal_crossovers_vs_age.pdf",sep=""),height=3,width=15)
par(mfrow=c(1,5))
for(i in 1:5){
NFTOOLS=all.datasets[[i]]
prefix=names(all.datasets)[i]
NFTOOLS.mat=data.frame(NFTOOLS[,c("child","mat.count","mat.age")])
mat.regress = lm(as.numeric(NFTOOLS.mat$mat.count)~NFTOOLS.mat$mat.age)
mat.predict = predict(mat.regress, interval="confidence")
plot(as.numeric(NFTOOLS.mat$mat.age),NFTOOLS.mat$mat.count,main=prefix,xlab="maternal age",ylab="maternal crossovers",col=alpha("red",0.4),bg=alpha("red",0.4),pch=21,ylim=c(19,80))
abline(mat.regress,col="black")
legend("topleft",paste0("slope = ",round(mat.regress$coefficients[2],digits=2)),cex=1.2)
legend("topright",paste0("p = ",round(summary(mat.regress)$coefficients[2,4],digits=3)),cex=1.2)
}
dev.off()



pdf(paste(mydir,"/","NFTOOLS_k_5_paternal_crossovers_vs_age.pdf",sep=""),height=3,width=15)
par(mfrow=c(1,5))
for(i in 1:5){
NFTOOLS=all.datasets[[i]]
prefix=names(all.datasets)[i]
NFTOOLS.pat=data.frame(NFTOOLS[,c("child","pat.count","pat.age")])
pat.regress = lm(as.numeric(NFTOOLS.pat$pat.count)~NFTOOLS.pat$pat.age)
pat.predict = predict(mat.regress, interval="confidence")
plot(as.numeric(NFTOOLS.pat$pat.age),NFTOOLS.pat$pat.count,main=prefix,xlab="paternal age",ylab="paternal crossovers",col=alpha("blue",0.4),bg=alpha("blue",0.4),pch=21,ylim=c(12,50))
abline(pat.regress,col="black")
legend("topleft",paste0("slope = ",round(pat.regress$coefficients[2],digits=2)),cex=1.2)
legend("topright",paste0("p = ",round(summary(pat.regress)$coefficients[2,4],digits=3)),cex=1.2)
#if(i==1){
#NFTOOLS.pat=NFTOOLS.pat[NFTOOLS.pat$pat.age<50,]
#pat.regress = lm(as.numeric(NFTOOLS.pat$pat.count)~NFTOOLS.pat$pat.age)
#pat.predict = predict(mat.regress, interval="confidence")
#legend("topright",c("Without outliers:",paste0("slope = ",round(pat.regress$coefficients[2],digits=2)),paste0("p = ",round(summary(pat.regress)$coefficients[2,4],digits=3))))
#}

}
dev.off()

