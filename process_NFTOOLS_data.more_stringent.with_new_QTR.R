cohorts=c("NTR","QTR370","QTR610","QTRCoreExome","FC","GPC","VB","ORCADES")
dirs=c("NTR","QTR","QTR","QTR","FC","GPC","VB","ORCADES")

events.files=c()
counts.files=c()

for(c in 1:length(cohorts)){
events.files=c(events.files,paste0(dirs[c],"/NFTOOLS/output_data/",cohorts[c],".more_stringent.all_chrs.inrecomb_k5.list.txt.events"))
counts.files=c(counts.files,paste0(dirs[c],"/NFTOOLS/final_data/",cohorts[c],".more_stringent.all_chrs.k_5.family_stats.txt"))
}

events.files[2:4] = c("QTR/NFTOOLS/output_data/all_QTR_370K_samples.all_chrs.no_errors.inrecomb_k5.list.txt.events",
"QTR/NFTOOLS/output_data/all_QTR_610K_samples.all_chrs.no_errors.inrecomb_k5.list.txt.events",
"QTR/NFTOOLS/output_data/new_samples.CoreExome.all_chrs.no_errors.inrecomb_k5.list.txt.events")

counts.files[2:4] = c("QTR/NFTOOLS/final_data/all_QTR_370K_samples.all_chrs.no_errors.k_5.family_stats.txt",
"QTR/NFTOOLS/final_data/all_QTR_610K_samples.all_chrs.no_errors.k_5.family_stats.txt",
"QTR/NFTOOLS/final_data/new_samples.CoreExome.all_chrs.no_errors.k_5.family_stats.txt")


fam.files=c("NTR/input_for_SHAPEIT2/NTR.1.more_stringent.generr_removed.fam",
"QTR/input_for_SHAPEIT2/all_QTR_370K_samples.1.post_MERLIN.more_stringent.generr_removed.fam",
"QTR/input_for_SHAPEIT2/all_QTR_610K_samples.1.post_MERLIN.more_stringent.generr_removed.fam",
"QTR/input_for_SHAPEIT2/new_samples.CoreExome.1.post_MERLIN.more_stringent.generr_removed.fam",
    "FC/input_for_SHAPEIT2/FC.1.more_stringent.generr_removed.fam",    
"GPC/input_for_SHAPEIT2/GPC.1.more_stringent.generr_removed.fam",    
"VB/input_for_SHAPEIT2/VB.1.more_stringent.generr_removed.fam",
    "ORCADES/input_for_SHAPEIT2/ORCADES.1.more_stringent.generr_removed.fam")
    
names(counts.files)=c("NTR","QTR370","QTR610","QTRCoreExome","FC","GPC","VB","ORCADES") 
names(events.files)=names(counts.files)
names(fam.files)=names(counts.files)
counts.by.parents=list()
all.events=list()
remove.outliers=function(r,fam,cohort){
        #remove families in FC cohort due to MZ twins/duplicates (76), and because of low maternal counts, reason unknown (52)
    if(cohort=="FC"){
        r=r[!r$child %in% fam[fam[,1] %in% c("52", "76"),2],]
    }
              #remove families FVG to due to pedigree erros
    if(cohort=="FVG"){
        r=r[!r$child %in% fam[fam[,1] %in% c("69", "148"),2],]
    }
                                        #remove families in VB: 143 and 64 each contain a cryptic pair of MZ twins or duplicate samples
    if(cohort=="VB"){
        r=r[!r$child %in% fam[fam[,1] %in% c("64", "143"),2],]
    }
     #remove families from NTR: 10689 seems to have three times the same sample; 11880 is strange (maybe too much missingness in mother)
    if(cohort=="NTR_v1"|cohort=="NTR_v2"){
        r=r[!r$child %in% fam[fam[,1] %in% c("10689", "11880"),2],]
    }
  #remove dodgy family 188 (mother APP5117398, father APP5117401 or APP5117402_dummyfather from GPC); think relationships are mis-specified
    if(cohort=="GPC"){
        r=r[!r$child %in% c("APP5117399","APP5119341"),]
    }
    return(r)
}
all.fams=list()

outliers=list()
for(c in 1:length(counts.files)){
    counts=read.delim(counts.files[c],header=T,stringsAsFactors=F,skip=2,colClasses=c("integer",rep("character",3),rep("numeric",8)))
    counts.by.parents[[names(counts.files)[c]]]=counts
    fam=read.delim(fam.files[c],header=F,colClasses=c(rep("character",4),rep("numeric",2)),sep="")
    all.fams[[names(counts.files)[c]]]=fam
    events=read.delim(events.files[c],header=T,stringsAsFactors=F,
colClasses=c("character","integer","integer","character","character","integer","character","integer"))
    events=remove.outliers(events,fam,names(counts.files)[c])
    events.by.child=as.data.frame(table(events[,c("child","sex")]))
    events.by.child$family=events[match(events.by.child$child,events$child),"fam"]
    family.count=table(events.by.child$family)/2
    events.by.child$nkids = family.count[as.character( events.by.child$family)]
    outliers[[names(counts.files)[c]]]=events.by.child[(events.by.child$nkids>2)&((events.by.child$sex=="male" & (events.by.child$Freq<9 |events.by.child$Freq > 45))|
(events.by.child$sex=="female" &
       (events.by.child$Freq<10|events.by.child$Freq > 76)))|(events.by.child$nkids==1) &  ((events.by.child$sex=="male" &
      (events.by.child$Freq<2*9 |events.by.child$Freq > 2*45))|(events.by.child$sex=="female" &(events.by.child$Freq<2*10|events.by.child$Freq > 2*76))),]
    events.by.child=events.by.child[!(events.by.child$nkids>2 & ((events.by.child$sex=="male" & (events.by.child$Freq<9 |events.by.child$Freq > 45))|
(events.by.child$sex=="female" &
         (events.by.child$Freq<10 |events.by.child$Freq > 76)))|(events.by.child$nkids==1) &  
((events.by.child$sex=="male" & (events.by.child$Freq<2*9 |events.by.child$Freq > 2*45))|(events.by.child$sex=="female"
      &(events.by.child$Freq<2*10|events.by.child$Freq > 2*76)))),]
    events.by.child$cohort=names(counts.files)[c]
    all.events[[names(counts.files)[c]]]=events.by.child

}

qtr370 = all.events[["QTR370"]]
qtr370$duo = paste(qtr370$child,qtr370$sex,sep="-")

qtr610 = all.events[["QTR610"]]
qtr610$duo = paste(qtr610$child,qtr610$sex,sep="-")

qtrCoreExome = all.events[["QTRCoreExome"]]
qtrCoreExome$duo = paste(qtrCoreExome$child,qtrCoreExome$sex,sep="-")

qtr370.keep  = qtr370[!qtr370$family %in% qtr370[qtr370$duo %in% qtr610$duo,"family"],]
qtrCoreExome.keep = qtrCoreExome[!qtrCoreExome$family %in% qtrCoreExome[qtrCoreExome$duo %in% qtr610$duo|qtrCoreExome$duo %in% qtr370$duo ,"family"],]

qtr370.chuck  = qtr370[qtr370$family %in% qtr370[qtr370$duo %in% qtr610$duo,"family"],]
qtrCoreExome.chuck = qtrCoreExome[qtrCoreExome$family %in% qtrCoreExome[qtrCoreExome$duo %in% qtr610$duo|qtrCoreExome$duo %in% qtr370$duo ,"family"],]

all.events[["QTR370"]]=qtr370.keep[,1:6]
all.events[["QTRCoreExome"]]=qtrCoreExome.keep[,1:6]

informative = rbind(unlist(lapply(all.events,function(events){sum(events$sex=="female" & events$nkids>2)})),unlist(lapply(all.events,function(events){
sum(events$sex=="male" & events$nkids>2)})))
rowSums(informative[,-1])

save(all.events,file="NFTOOLS_results/NFTOOLS_counts.all_families_except_outliers.more_stringent.with_new_QTR.RData")

qqr <- function(k1,rate,col1,makeplot=TRUE,title,limits,xlowerlimit,...) {
      k1 <- sort(k1)
        n <- length(k1)
        e <- qpois(1:n/(n+1),rate)
        l <- c(min(k1),min(max(k1),rate+6*sqrt(rate)))
        etmp <- e[k1>=l[2]]
        ret <- cbind(jitter(etmp),l[2])

        if(makeplot) {
                 # plot(e,k1,type='n',xlab="",ylab="",ylim=l,main=title,...)
 #                  plot(e,k1,type='n',xlab="",ylab="",ylim=limits,xlim=c(xlowerlimit,max(limits)),main=title,...)
                              plot(e,k1,type='n',xlab="",ylab="",ylim=limits,xlim=limits,main=title,...)
                         abline(0,1,col=1,lty=2);grid()
                         points(jitter(e[k1<l[2]]),jitter((k1[k1<l[2]])),col=col1,pch=16)
#                         points(ret,col=col1,pch=17)
               }

        return(  cbind(e,k1))
  }
male.lengths <- c(1.9512,1.8955,1.6071,1.4654,1.512,1.3762,1.2835,1.0794,1.1725,1.3389,1.0936,1.3554,1.0131,0.9462,1.0257,1.081,1.0856,0.9862,0.9264,0.7472,0.4731,
0.4896)
female.lengths <- c(3.4541,3.2541,2.7564,2.5906,2.6019,2.4159,2.3033,2.0994,1.982,2.1813,1.9553,2.0664,1.5588,1.4236,1.5496,1.4962,1.6153,1.4257,1.2682,1.2197,0.764,
0.8276)
malesize <- sum(male.lengths)
femalesize <- sum(female.lengths)

for(i in 1:length(all.events)){
    rec=all.events[[i]]
    rec.2kids=rec[rec$nkids==1,]
    rec.2kids$Freq=rec.2kids$Freq/2
    rec=rec[rec$nkids>2,]
#    pdf(paste0("/home/hilary/maternal_age_recombination/qqplots_NFTOOLS_more_stringent/",names(all.events)[i],"-qq.pdf"),width=12,height=6)
 pdf(paste0("NFTOOLS_results/",names(all.events)[i],"-qq.pdf"),width=12,height=6)
    par(mar=c(2,2,1,1),mfrow=c(1,2),cex=1,oma=c(2,2,0,0))
    tmp <- qqr(subset(rec,sex=="male")$Freq,malesize,"red",title=paste(names(all.events)[i]," paternal"),
               limits=c(9,45),xlowerlimit=10)
    points(jitter(qqr(subset(rec.2kids,sex=="male")$Freq,malesize,makeplot=0)),pch=16,col="darkgreen")
    legend("topleft",c("phase known","phased unknown"),pch=16,col=c("red","darkgreen"))
    cat(names(all.events)[i],"\t",range(subset(rec,sex=="male")$Freq),"\n")
    if(sum(rec$sex=="male" & rec$Freq>50)>0){
       print(rec[rec$sex=="male" & rec$Freq>50,])
    }
    tmp <- qqr(subset(rec,sex=="female")$Freq,femalesize,"red",title=paste(names(all.events)[i]," maternal"),
               limits=c(10,76),xlowerlimit=20)
        points(jitter(qqr(subset(rec.2kids,sex=="female")$Freq,femalesize,makeplot=0)),pch=16,col="darkgreen")  
    cat(names(all.events)[i],"\t",range(subset(rec,sex=="female")$Freq),"\n")
    if(sum(rec$sex=="female" & rec$Freq>70)>0){
        print(rec[rec$sex=="female" & rec$Freq>70,])
    }
#    points(jitter(qqr(subset(rec.flt5,sex=="female")$Freq,femalesize,makeplot=0)),pch=16,col=cols[5])
    mtext("Expected recombinations genome-wide",1,outer=T,cex=1.6,padj=1.7)
    mtext("Observed recombinations genome-wide",2,outer=T,cex=1.6,padj=-1)
    dev.off()
}


#plotting histograms of distributions

for(i in 1:length(all.events)){
#    pdf(paste0("/home/hilary/maternal_age_recombination/qqplots_NFTOOLS_more_stringent/",names(all.events)[i],".histograms.pdf"),width=8,height=8)
    pdf(paste0("NFTOOLS_results/",names(all.events)[i],".histograms.pdf"),width=8,height=8)

    par(mfrow=c(2,2))
    rec=all.events[[i]]
    rec.2kids=rec[rec$nkids==1,]
    rec.2kids$Freq=rec.2kids$Freq/2
    rec=rec[rec$nkids>2,]
    hist(subset(rec,sex=="male")$Freq,main="male, phase known",xlab="Number of recombinations",col="lightblue")
    abline(v=25.9,col="darkgreen",lwd=2)
    abline(v=26.6,col="blue",lwd=2)
    abline(v=mean(subset(rec,sex=="male")$Freq),col="black",lwd=2)
     legend("topleft",c("deCODE","Adam","these data"),col=c("darkgreen","blue","black"),lwd=2,lty=1)

    hist(subset(rec.2kids,sex=="male")$Freq,main="male, phase unknown",xlab="Number of recombinations",col="lightblue")
    abline(v=25.9,col="darkgreen",lwd=2)
    abline(v=26.6,col="blue",lwd=2)
    abline(v=mean(subset(rec.2kids,sex=="male")$Freq),col="black",lwd=2)
    abline(v=mean(subset(rec,sex=="male")$Freq),col="black",lwd=2,lty=2)
    legend("topleft","phase known",col="black",lwd=2,lty=2)
    
    hist(subset(rec,sex=="female")$Freq,main="female, phase known",xlab="Number of recombinations",col="pink")
    abline(v=42.81,col="darkgreen",lwd=2)
    abline(v=41.6,col="blue",lwd=2)
    abline(v=mean(subset(rec,sex=="female")$Freq),col="black",lwd=2)

    hist(subset(rec.2kids,sex=="female")$Freq,main="female, phase unknown",xlab="Number of recombinations",col="pink")
    abline(v=42.81,col="darkgreen",lwd=2)
    abline(v=41.6,col="blue",lwd=2)
    abline(v=mean(subset(rec.2kids,sex=="female")$Freq),col="black",lwd=2)
    abline(v=mean(subset(rec,sex=="female")$Freq),col="black",lwd=2,lty=2)
    legend("topleft","phase known",col="black",lwd=2,lty=2)

    dev.off()
}
### compare NFTOOLS to duoHMM

duoHMM.files=c("NTR/duoHMM_results/NTR.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"QTR/duoHMM_results/all_QTR_370K_samples.all.post_MERLIN.more_stringent.generr_removed.recombination_count.all_duos.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"QTR/duoHMM_results/all_QTR_610K_samples.all.post_MERLIN.more_stringent.generr_removed.recombination_count.all_duos.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"QTR/duoHMM_results/new_samples.CoreExome.all.post_MERLIN.more_stringent.generr_removed.recombination_count.all_duos.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
    "FC/duoHMM_results/FC.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
    "GPC/duoHMM_results/GPC.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
    "VB/duoHMM_results/VB.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
    "ORCADES/duoHMM_results/ORCADES.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt")
    
names(duoHMM.files)=c("NTR","QTR370","QTR610","QTRCoreExome","FC","GPC","VB","ORCADES")
    
age.files=c("NTR/NTR_v2_minus_MZs.with_parental_age.txt","QTR/all_QTR_370K_samples.with_parental_age.txt","QTR/all_QTR_610K_samples.with_parental_age.txt",
"QTR/new_samples.CoreExome.with_parental_age.txt","FC/FC.with_parental_age.txt",    "GPC/GPC.with_parental_age.txt",    "VB/valborbera_b37-related.with_parental_age.txt",
"ORCADES/ORCADES.with_parental_age.txt")

names(age.files)=names(all.events)
library(scales)
NFTOOLS.inf.mat.data.with.age=NULL
NFTOOLS.uninf.mat.data.with.age=NULL
NFTOOLS.inf.pat.data.with.age=NULL
NFTOOLS.uninf.pat.data.with.age=NULL
for(i in 1:length(all.events)){
    events=all.events[[i]]
    events$cohort=names(all.events)[i]
    ages=read.delim(age.files[i],header=T,stringsAsFactors=F)
    ages$individual=gsub("-","_",ages$individual)
    ages$mother=gsub("-","_",ages$mother)
    ages$father=gsub("-","_",ages$father)
    mat.events.inf=events[events$sex=="female" & events$nkids>2,]
    pat.events.inf=events[events$sex=="male" & events$nkids>2,]
    mat.events.uninf=events[events$sex=="female" & events$nkids==1,]
    pat.events.uninf=events[events$sex=="male" & events$nkids==1,]
    mat.events.inf$age=ages[match(mat.events.inf$child,ages$individual),"maternal_age"]
    pat.events.inf$age=ages[match(pat.events.inf$child,ages$individual),"paternal_age"]
    mat.events.uninf$sib1.age=ages[match(mat.events.uninf$child,ages$individual),"maternal_age"]
    pat.events.uninf$sib1.age=ages[match(pat.events.uninf$child,ages$individual),"paternal_age"]
   
    mat.events.uninf$sib2.age=NA
    pat.events.uninf$sib2.age=NA
    mat.events.uninf$sib2=NA
    pat.events.uninf$sib2=NA
    for(j in 1:nrow(mat.events.uninf)){
        mat=ages[ages$individual==mat.events.uninf[j,"child"],"mother"]
        sibs=ages[ages$mother == mat,"individual"]
        sib=sibs[1]
        if(sibs[1]==mat.events.uninf[j,"child"]){            sib=sibs[2]        }
        mat.events.uninf$sib2[j]=sib
        mat.events.uninf$sib2.age[j]=ages[ages$individual==sib,"maternal_age"]
    }
    for(j in 1:nrow(pat.events.uninf)){
        pat=ages[ages$individual==pat.events.uninf[j,"child"],"father"]
        sibs=ages[ages$father == pat,"individual"]
        sib=sibs[1]
        if(sibs[1]==pat.events.uninf[j,"child"]){            sib=sibs[2]        }
        pat.events.uninf$sib2[j]=sib
        pat.events.uninf$sib2.age[j]=ages[ages$individual==sib,"paternal_age"]
        
    }
    mat.events.uninf$av.Freq=mat.events.uninf$Freq/2
    pat.events.uninf$av.Freq=pat.events.uninf$Freq/2
    NFTOOLS.inf.mat.data.with.age=rbind(NFTOOLS.inf.mat.data.with.age,mat.events.inf)
    NFTOOLS.inf.pat.data.with.age=rbind(NFTOOLS.inf.pat.data.with.age,pat.events.inf)
    NFTOOLS.uninf.mat.data.with.age=rbind(NFTOOLS.uninf.mat.data.with.age,mat.events.uninf)
    NFTOOLS.uninf.pat.data.with.age=rbind(NFTOOLS.uninf.pat.data.with.age,pat.events.uninf)
    

    duoHMM=read.delim(duoHMM.files[names(all.events)[i]],header=T,stringsAsFactors=F)
    duoHMM.mat=duoHMM[duoHMM$sex=="Female",]
    duoHMM.pat=duoHMM[duoHMM$sex=="Male",]
    mat.events.inf$duoHMM.Freq=duoHMM.mat[match(mat.events.inf$child,duoHMM.mat$CHILD),"nrec"]
    pat.events.inf$duoHMM.Freq=duoHMM.pat[match(pat.events.inf$child,duoHMM.pat$CHILD),"nrec"]

    mat.events.uninf$duoHMM.sib1.Freq=duoHMM.mat[match(mat.events.uninf$child,duoHMM.mat$CHILD),"nrec"]
    pat.events.uninf$duoHMM.sib1.Freq=duoHMM.pat[match(pat.events.uninf$child,duoHMM.pat$CHILD),"nrec"]

    mat.events.uninf$duoHMM.sib2.Freq=duoHMM.mat[match(mat.events.uninf$sib2,duoHMM.mat$CHILD),"nrec"]
    pat.events.uninf$duoHMM.sib2.Freq=duoHMM.pat[match(pat.events.uninf$sib2,duoHMM.pat$CHILD),"nrec"]

    mat.events.uninf$duoHMM.total.Freq=mat.events.uninf$duoHMM.sib1.Freq + mat.events.uninf$duoHMM.sib2.Freq
    pat.events.uninf$duoHMM.total.Freq=pat.events.uninf$duoHMM.sib1.Freq + pat.events.uninf$duoHMM.sib2.Freq

#    pdf(paste0("/home/hilary/maternal_age_recombination/qqplots_NFTOOLS_more_stringent/",names(all.events)[i],".comparing_NFTOOLS_to_duoHMM.informative_vs_2kid_families.pdf"),width=8,height=8)
    pdf(paste0("NFTOOLS_results/",names(all.events)[i],".comparing_NFTOOLS_to_duoHMM.informative_vs_2kid_families.pdf"),width=8,height=8)
    par(mfrow=c(2,2))
    plot(jitter(mat.events.inf$Freq),jitter(mat.events.inf$duoHMM.Freq),xlab="NFTOOLS",ylab="duoHMM",main="Maternal, phase known",col=alpha("red",0.5),pch=19)
    abline(a=0,b=1)
    plot(jitter(pat.events.inf$Freq),jitter(pat.events.inf$duoHMM.Freq),xlab="NFTOOLS",ylab="duoHMM",main="Paternal, phase known",col=alpha("blue",0.5),pch=19)
    abline(a=0,b=1)
    plot(jitter(mat.events.uninf$Freq),jitter(mat.events.uninf$duoHMM.total.Freq),xlab="NFTOOLS total for 2 kids",ylab="duoHMM total for 2 kids",main="Maternal, phase unknown",
col=alpha("red",0.5),pch=19)
    abline(a=0,b=1)
    plot(jitter(pat.events.uninf$Freq),jitter(pat.events.uninf$duoHMM.total.Freq),xlab="NFTOOLS total for 2 kids",ylab="duoHMM total for 2 kids",main="Paternal, phase unknown",
col=alpha("blue",0.5),pch=19)
    abline(a=0,b=1)
    dev.off()

}

#summaries of distributions
cohort.distributions=NULL
summarise.distribution=function(x){
    return(c(length(x),mean(x),quantile(x,c(0,0.25,0.5,0.75,1)),sd(x)))
}

for(i in 1:length(all.events)){
    rec=all.events[[i]]
    rec.2kids=rec[rec$nkids==1,]
    rec.2kids$Freq=rec.2kids$Freq/2
    rec=rec[rec$nkids>2,]
    cohort.distributions=rbind(cohort.distributions,c(summarise.distribution(subset(rec,sex=="male")$Freq),summarise.distribution( subset(rec.2kids,sex=="male")$Freq),
summarise.distribution( subset(rec,sex=="female")$Freq),        summarise.distribution( subset(rec.2kids,sex=="female")$Freq)))
}

rownames(cohort.distributions)=names(all.events)
colnames(cohort.distributions)=sapply(c("male.phase.known","male.phase.unknown","female.phase.known","female.phase.unknown"),function(x){paste0(x,c(".n",".mean",".min",".25th",".median",".75th",".max",".sd"))})

#write.table(cohort.distributions,"/home/hilary/maternal_age_recombination/qqplots_NFTOOLS_more_stringent/summary_of_distributions.NFTOOLS.txt",quote=F,sep="\t")
write.table(cohort.distributions,"NFTOOLS_results/summary_of_distributions.NFTOOLS.with_new_QTR.txt",quote=F,sep="\t")
