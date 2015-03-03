options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

recomb.file=args[1]#recomb.file="FC/duoHMM_results/FC.1.more_stringent.generr_removed.recombinations.txt"

#recomb.file="QTR/duoHMM_results/QTR_minus_MZs.370K.1.post_MERLIN.more_stringent.generr_removed.recombinations.txt"
sample.file = args[2]#sample.file="FC/SHAPEIT_results/FC.1.more_stringent.generr_removed.SHAPEIT_phased.sample"
#sample.file="QTR/SHAPEIT_results/QTR_minus_MZs.370K.1.post_MERLIN.more_stringent.generr_removed.SHAPEIT_phased.sample"
out.prefix = args[3]#out.prefix="qqplots_duoHMM_v3/FC.more_stringent.informative_meioses.qqplots."
#out.prefix="QTR/duoHMM_results/QTR370.more_stringent.informative_meioses.qqplots."
cohort=args[4]#cohort="FC"

samples=read.delim(sample.file,header=F,sep="",stringsAsFactors=F,skip=2)
samples=samples[,-c(3,7)]
samples2=cbind(as.character(samples[,1]),as.character(samples[,2]),as.character(samples[,3]),as.character(samples[,4]),as.character(samples[,5]))
rownames(samples2) = samples2[,2]
samples2=as.data.frame(samples2,stringsAsFactors=F)
#remove parents not genotyped
samples2[! samples[,3] %in% samples[,2],3]="0"
samples2[! samples[,4] %in% samples[,2],4]="0"

samples2$total_fam_size=sapply(samples2[,1],function(x){sum(samples2[,1]==x)})
samples2$total_same_father=sapply(samples2[,2],function(x){sum(samples2[,3]==samples2[x,3])})
samples2$total_same_mother=sapply(samples2[,2],function(x){sum(samples2[,4]==samples2[x,4])})
samples2$total_same_parents=sapply(samples2[,2],function(x){sum(samples2[,4]==samples2[x,4] & samples2[,3]==samples2[x,3])})
samples2$total_same_father[samples2[,3]=="0"]=0
samples2$total_same_mother[samples2[,4]=="0"]=0
samples2$total_same_parents[samples2[,4]=="0" & samples2[,3]=="0"]=0

samples2$maternal_grandmother="0"
samples2$maternal_grandmother[samples2[,4] %in% samples2[,2]] = samples2[samples2[samples2[,4]%in% samples2[,2],4],4]
samples2$maternal_grandfather="0"
samples2$maternal_grandfather[samples2[,4] %in% samples2[,2]] = samples2[samples2[samples2[,4]%in% samples2[,2],4],3]

samples2$paternal_grandmother="0"
samples2$paternal_grandmother[samples2[,3] %in% samples2[,2]] = samples2[samples2[samples2[,3]%in% samples2[,2],3],4]
samples2$paternal_grandfather="0"
samples2$paternal_grandfather[samples2[,3] %in% samples2[,2]] = samples2[samples2[samples2[,3]%in% samples2[,2],3],3]

samples2$n_genotyped_maternal_grandparents =rowSums(samples2[,c("maternal_grandmother","maternal_grandfather")]!="0")
samples2$n_genotyped_paternal_grandparents =rowSums(samples2[,c("paternal_grandmother","paternal_grandfather")]!="0")
samples2$n_genotyped_grandparents =samples2$n_genotyped_maternal_grandparents+samples2$n_genotyped_paternal_grandparents
samples2$n_genotyped_children = sapply(samples2[,2],function(x){sum(samples[,3] == x | samples[,4] == x)})

halfsibs=samples2[(samples2$total_same_parents!=samples2$total_same_mother&samples2[,4]!="0")|(samples2$total_same_parents!=samples2$total_same_father & samples2[,3]!="0"),]
samples2$possible_halfsib=((samples2$total_same_parents!=samples2$total_same_mother&samples2[,4]!="0")|(samples2$total_same_parents!=samples2$total_same_father & samples2[,3]!="0"))

all.recomb=NULL
for(i in 1:22){
recomb.file2=gsub("1.",paste(i,".",sep=""),recomb.file,fixed=T)

recomb=read.delim(recomb.file2,header=T,sep="\t",stringsAsFactors=F)
recomb$chr=i
all.recomb=rbind(all.recomb,recomb)

}

rec=read.delim( gsub("recombinations","recombination_count.all_duos.p_0.5.annotated.double_xovers_within_X_SNPs_removed",gsub("1.","all.",recomb.file,fixed=T)),header=T)
rec=rec[!is.na(rec$informative),]

all.recomb$CHILD=gsub("-","_",all.recomb$CHILD)
all.recomb$PARENT=gsub("-","_",all.recomb$PARENT)

all.recomb$duo=paste0(all.recomb$CHILD,"-",all.recomb$PARENT)
all.recomb=all.recomb[all.recomb$duo %in% rec$duo,]
all.recomb=cbind(all.recomb,rec[match(all.recomb$duo,rec$duo),c("informative" , "informative.2gen" ,  "informative.2gen.2parents" ,"informative.2gen.1parent"   , "noninf.2kids.2gen.2parents","noninf.2kids.2gen.1parent"
   ,"noninf.1kid.2gen.2parents" ,"noninf.1kid.2gen.1parent"  , "informative.3gen"  ,"informative.3gen.2parents" ,"informative.3gen.1parent"  )])


qqr <- function(k1,rate,col1,makeplot=TRUE,title,limits,xlowerlimit,...) {
      k1 <- sort(k1)
        n <- length(k1)
        e <- qpois(1:n/(n+1),rate)
        l <- c(min(k1),min(max(k1),rate+6*sqrt(rate)))
        etmp <- e[k1>=l[2]]
        ret <- cbind(jitter(etmp),l[2])

        if(makeplot) {
            #       plot(e,k1,type='n',xlab="",ylab="",ylim=limits,xlim=c(xlowerlimit,max(limits)),main=title,...)
                         plot(e,k1,type='n',xlab="",ylab="",ylim=limits,xlim=limits,main=title,...)
                               abline(0,1,col=1,lty=2);grid()
                               points(jitter(e[k1<l[2]]),jitter((k1[k1<l[2]])),col=col1,pch=16)
                         #      points(ret,col=col1,pch=17)
                     }

        return(  cbind(e,k1))
  }

samples2$V2=gsub("-","_",samples2$V2)
samples2$V3=gsub("-","_",samples2$V3)
samples2$V4=gsub("-","_",samples2$V4)

maternal=samples2[samples2[,4] != "0",c(1,2,4,3,6:ncol(samples2))]
paternal=samples2[samples2[,3] != "0",c(1,2,3,4,6:ncol(samples2))]

rownames(maternal)=maternal[,2]
rownames(paternal)=paternal[,2]

colnames(maternal)[1:4]=c("family","child","mother","father")
colnames(paternal)[1:4]=c("family","child","father","mother")

maternal$duo=paste0(maternal$child,"-",maternal$mother)

paternal$duo=paste0(paternal$child,"-",paternal$father)

maternal=maternal[maternal$duo %in% all.recomb$duo,]
paternal=paternal[paternal$duo %in% all.recomb$duo,]


maternal=maternal[maternal$duo %in% rec$duo,]
paternal=paternal[paternal$duo %in% rec$duo,]

maternal=cbind(maternal,rec[match(maternal$duo,rec$duo),c("informative" , "informative.2gen" ,  "informative.2gen.2parents" ,"informative.2gen.1parent"   , "noninf.2kids.2gen.2parents","noninf.2kids.2gen.1parent"
   ,"noninf.1kid.2gen.2parents" ,"noninf.1kid.2gen.1parent"  , "informative.3gen"  ,"informative.3gen.2parents" ,"informative.3gen.1parent"  )])

paternal=cbind(paternal,rec[match(paternal$duo,rec$duo),c("informative" , "informative.2gen" ,  "informative.2gen.2parents" ,"informative.2gen.1parent"   , "noninf.2kids.2gen.2parents","noninf.2kids.2gen.1parent"
   ,"noninf.1kid.2gen.2parents" ,"noninf.1kid.2gen.1parent"  , "informative.3gen"  ,"informative.3gen.2parents" ,"informative.3gen.1parent"  )])


colours=c("black","blue","red","green","darkgreen","skyblue","red4","purple","orange","pink3")
male.lengths <- c(1.9512,1.8955,1.6071,1.4654,1.512,1.3762,1.2835,1.0794,1.1725,1.3389,1.0936,1.3554,1.0131,0.9462,1.0257,1.081,1.0856,0.9862,0.9264,0.7472,0.4731,0.4896)
female.lengths <- c(3.4541,3.2541,2.7564,2.5906,2.6019,2.4159,2.3033,2.0994,1.982,2.1813,1.9553,2.0664,1.5588,1.4236,1.5496,1.4962,1.6153,1.4257,1.2682,1.2197,0.764,0.8276)
malesize <- sum(male.lengths)
femalesize <- sum(female.lengths)


summarise.counts=function(counts){
    mysummary=c(length(counts),mean(counts),quantile(counts,c(0,0.25,0.5,0.75,1)),sd(counts))
    names(mysummary)=c("N","mean","min","25th_pc","median","75th_pc","max")
    return(mysummary)
}
summary.counts=NULL
#cutoffs= seq(from=0,to=0.9,by=0.1)
#cutoffs= seq(from=0,to=0.9,by=0.3)
cutoffs=c(0,0.3,0.5,0.7,0.9)
for(i in 1:length(cutoffs)){
    p=cutoffs[i]
    recomb2=all.recomb[all.recomb$PROB_RECOMBINATION > p,]
    maternal=data.frame(maternal,rep(0,nrow(maternal)),rep(0,nrow(maternal)),stringsAsFactors=F)
    colnames(maternal)[ncol(maternal)-1]=paste("E.nrec.p.",p,sep="")
    colnames(maternal)[ncol(maternal)]=paste("nrec.p.",p,sep="")
    maternal.counts=apply(maternal[,2:3],1,function(indivs){sum(recomb2[recomb2[,1]==indivs[1] & recomb2[,2]==indivs[2],"PROB_RECOMBINATION"])})
    maternal.counts2=apply(maternal[,2:3],1,function(indivs){nrow(recomb2[recomb2[,1]==indivs[1] & recomb2[,2]==indivs[2],])})
    maternal[names(maternal.counts),ncol(maternal)-1]=maternal.counts
    maternal[names(maternal.counts),ncol(maternal)]=maternal.counts2
    summary.counts=rbind(summary.counts,summarise.counts(maternal.counts[maternal$informative]),summarise.counts(maternal.counts2[maternal$informative]))
}

for(i in 1:length(cutoffs)){
    p=cutoffs[i]
    recomb2=all.recomb[all.recomb$PROB_RECOMBINATION > p,]
    paternal=data.frame(paternal,rep(0,nrow(paternal)),rep(0,nrow(paternal)),stringsAsFactors=F)
    colnames(paternal)[ncol(paternal)-1]=paste("E.nrec.p.",p,sep="")
    colnames(paternal)[ncol(paternal)]=paste("nrec.p.",p,sep="")
    paternal.counts=apply(paternal[,2:3],1,function(indivs){sum(recomb2[recomb2[,1]==indivs[1] & recomb2[,2]==indivs[2],"PROB_RECOMBINATION"])})
    paternal.counts2=apply(paternal[,2:3],1,function(indivs){nrow(recomb2[recomb2[,1]==indivs[1] & recomb2[,2]==indivs[2],])})
    paternal[names(paternal.counts),ncol(paternal)-1]=paternal.counts
    paternal[names(paternal.counts),ncol(paternal)]=paternal.counts2
    summary.counts=rbind(summary.counts,summarise.counts(paternal.counts[paternal$informative]),summarise.counts(paternal.counts2[paternal$informative]))
}

write.table(maternal,paste0(out.prefix,".maternal_recombination_counts_different_p_thresholds.txt"),quote=F,sep="\t")
write.table(paternal,paste0(out.prefix,".paternal_recombination_counts_different_p_thresholds.txt"),quote=F,sep="\t")

rownames(summary.counts) =paste0(c(rep("mat",2*length(cutoffs)),rep("pat",2*length(cutoffs))),".",rep(c("E.nrec.p.","nrec.p."),2*length(cutoffs)),c(rep(cutoffs,each=2),rep(cutoffs,each=2)))
summary.counts=data.frame(summary.counts)
summary.counts$cohort=cohort
if(FALSE){
write.table(summary.counts,paste0(out.prefix,"_summary.txt"),quote=F,sep="\t")
cutoffs=c(0,0.3,0.5,0.7,0.9)
pdf(paste0(out.prefix,"qqplots.summing_probabilities.pdf"),height=5,width=10)
par(mar=c(2,2,1,1),mfrow=c(1,2),cex=1,oma=c(2,2,1,0))
for(i in 1:length(cutoffs)){
    p=cutoffs[i]
    maternal.counts=maternal[maternal$informative,paste("E.nrec.p.",p,sep="")]
    if(p==0){
        tmp <- qqr(maternal.counts,femalesize,colours[i],title=paste0(cohort, " maternal"),limits=range(20,76),xlowerlimit=20)
    } else{
        points(jitter(qqr(maternal.counts,femalesize,makeplot=0)),pch=16,col=colours[i])
    }
}
legend("topleft",paste0("Expected number for p > ",cutoffs),col=colours,pch=16,  bg="white",cex=0.8)
mtext("Expected recombinations genome-wide",1,outer=T,cex=1.3,padj=1.2)
mtext("Observed recombinations genome-wide",2,outer=T,cex=1.3,padj=-0.5)
for(i in 1:length(cutoffs)){
    p=cutoffs[i]
    paternal.counts=paternal[paternal$informative,paste("E.nrec.p.",p,sep="")]
    if(p==0){
        tmp <- qqr(paternal.counts,malesize,colours[i],title=paste0(cohort," paternal"),limits=range(9,45),xlowerlimit=10)
    } else {
        points(jitter(qqr(paternal.counts,malesize,makeplot=0)),pch=16,col=colours[i])
    }
}
dev.off()


pdf(paste0(out.prefix,"qqplots.number_over_threshold.pdf"),height=5,width=10)
par(mar=c(2,2,1,1),mfrow=c(1,2),cex=1,oma=c(2,2,1,0))
for(i in 1:length(cutoffs)){
    p=cutoffs[i]
    maternal.counts=maternal[maternal$informative,paste("nrec.p.",p,sep="")]
    if(p==0){
        tmp <- qqr(maternal.counts,femalesize,colours[i],title=paste0(cohort, " maternal"),limits=range(20,76),xlowerlimit=20)
    } else{
        points(jitter(qqr(maternal.counts,femalesize,makeplot=0)),pch=16,col=colours[i])
    }
}
legend("topleft",paste0("Number with p > ",cutoffs),col=colours,pch=16,  bg="white",cex=0.8)
mtext("Expected recombinations genome-wide",1,outer=T,cex=1.3,padj=1.2)
mtext("Observed recombinations genome-wide",2,outer=T,cex=1.3,padj=-0.5)
for(i in 1:length(cutoffs)){
    p=cutoffs[i]
    paternal.counts=paternal[paternal$informative,paste("nrec.p.",p,sep="")]
    if(p==0){
        tmp <- qqr(paternal.counts,malesize,colours[i],title=paste0(cohort," paternal"),limits=range(9,45),xlowerlimit=10)
    } else {
        points(jitter(qqr(paternal.counts,malesize,makeplot=0)),pch=16,col=colours[i])
    }
}
dev.off()
}

#write.table(maternal,paste(out.prefix,".maternal.summing_probabilities.txt",sep=""),quote=F,sep="\t",row.names=F)
#write.table(paternal,paste(out.prefix,".paternal.summing_probabilities.txt",sep=""),quote=F,sep="\t",row.names=F)



