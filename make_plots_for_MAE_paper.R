outdir="/home/hilary/maternal_age_recombination/paper/figures/"
load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.more_stringent.RData")
if(FALSE){
data1.mat=data0.mat
data1.pat=data0.pat
}

if(FALSE){#run for NFTOOLS
    rm(list=ls())
    load("/well/donnelly/hilary/maternal_age_and_recombination/NFTOOLS_data_for_RSTAN.more_stringent.RData")
    colnames(data1.mat)[3]="nrec"
    colnames(data1.pat)[3]="nrec"
}

#make box-and-whisker plots of duoHMM counts by binned age, for informative meioses
data1.mat$binned.age = 20
data1.mat$binned.age[data1.mat$age.at.birth >20 & data1.mat$age.at.birth<=25] = 25
data1.mat$binned.age[data1.mat$age.at.birth >25 & data1.mat$age.at.birth<=30] = 30
data1.mat$binned.age[data1.mat$age.at.birth >30 & data1.mat$age.at.birth<=35] = 35
data1.mat$binned.age[data1.mat$age.at.birth >35 & data1.mat$age.at.birth<=40] = 40
data1.mat$binned.age[data1.mat$age.at.birth >40] = 45

data1.pat$binned.age = 20
data1.pat$binned.age[data1.pat$age.at.birth >20 & data1.pat$age.at.birth<=25] = 25
data1.pat$binned.age[data1.pat$age.at.birth >25 & data1.pat$age.at.birth<=30] = 30
data1.pat$binned.age[data1.pat$age.at.birth >30 & data1.pat$age.at.birth<=35] = 35
data1.pat$binned.age[data1.pat$age.at.birth >35 & data1.pat$age.at.birth<=40] = 40
data1.pat$binned.age[data1.pat$age.at.birth >40] = 45

data1.mat$nrec.normalised=data1.mat$nrec-mean(data1.mat$nrec[data1.mat$binned.age==25])
data1.pat$nrec.normalised=data1.pat$nrec-mean(data1.pat$nrec[data1.pat$binned.age==25])

data1.mat.by.bin=split(data1.mat,data1.mat$binned.age)
data1.pat.by.bin=split(data1.pat,data1.pat$binned.age)

mat.binned.rec = as.data.frame(do.call("rbind",lapply(data1.mat.by.bin,function(rec){
mean.rec=    mean(rec$nrec.normalised)
lower.rec=mean.rec-qnorm(0.975)*sd(rec$nrec.normalised)/sqrt(length(rec$nrec))
upper.rec=mean.rec+qnorm(0.975)*sd(rec$nrec.normalised)/sqrt(length(rec$nrec))
return(c(mean.rec,lower.rec,upper.rec))
})))
colnames(mat.binned.rec)=c("mean","lower.95CI","upper.95CI")
pat.binned.rec = as.data.frame(do.call("rbind",lapply(data1.pat.by.bin,function(rec){
mean.rec=    mean(rec$nrec.normalised)
lower.rec=mean.rec-qnorm(0.975)*sd(rec$nrec.normalised)/sqrt(length(rec$nrec))
upper.rec=mean.rec+qnorm(0.975)*sd(rec$nrec.normalised)/sqrt(length(rec$nrec))
return(c(mean.rec,lower.rec,upper.rec))
})))
colnames(pat.binned.rec)=c("mean","lower.95CI","upper.95CI")
mat.binned.rec$x=c(17.5,22.5,27.5,32.5,37.5,42.5)
pat.binned.rec$x=c(17.5,22.5,27.5,32.5,37.5,42.5)+0.5

#pdf(paste0(outdir,"/crossover_count_by_binned_age.duoHMM_informative_duos_only.pdf"),height=5,width=5)
pdf(paste0(outdir,"/crossover_count_by_binned_age.duoHMM_informative_nuc_fams_only.pdf"),height=5,width=5)
#pdf(paste0(outdir,"/crossover_count_by_binned_age.NFTOOLS_informative_nuc_fams_only.pdf"),height=5,width=5)
plot(mat.binned.rec$x,mat.binned.rec$mean,xlab="Parental age",ylab="Difference from 20-25 age group",col="red",ylim=c(-4,4),pch=19)
points(pat.binned.rec$x,pat.binned.rec$mean,col="blue",pch=19)
apply(mat.binned.rec,1,function(bounds){segments(x0=bounds[4],y0=bounds[2],x1=bounds[4],y1=bounds[3],col="red")})
apply(pat.binned.rec,1,function(bounds){segments(x0=bounds[4],y0=bounds[2],x1=bounds[4],y1=bounds[3],col="blue")})
legend("topright",c("maternal","paternal"),col=c("red","blue"),pch=19,lty=1)
abline(h=0)
dev.off()

##### Meta-analysis
library(nlme)
library(metafor)
load("frequentist_analysis/more_stringent/frequentist_analysis_of_MAE.RData")

mycols=c("black","blue","red","green","orange","red4","purple","darkgreen","violetred2")
names(mycols)=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB","ORCADES")

#fixed effects meta-analysis of betas
#informative_nuclear_families_duoHMM
meta.fe=rma.uni(yi=summary.duohmm2.mat.results.by.cohort[1,-ncol(summary.duohmm2.mat.results.by.cohort)],sei=summary.duohmm2.mat.results.by.cohort[2,-ncol(summary.duohmm2.mat.results.by.cohort)],method="FE",weighted=T)
pdf("/home/hilary/maternal_age_recombination/paper/figures/box_plot_of_beta_Age_from_LMM_on_informative_nuclear_families_duoHMM.pdf",height=4,width=6)
plot(1:8,c(summary.duohmm2.mat.results.by.cohort[1,-ncol(summary.duohmm2.mat.results.by.cohort)],meta.fe$b),main="beta_age from LMM",ylim=c(-0.7,1),col=c(rep("white",7),"black"),pch=19,xaxt="n",xlab="Cohort",ylab="beta_age",cex=0.5)
axis(1,at=1:8,labels=FALSE,tick=F,cex.axis=0.5)
mylabels=c(gsub("duohmm.in.nftools.","",colnames(summary.duohmm2.mat.results.by.cohort)[-ncol(summary.duohmm2.mat.results.by.cohort)]),"meta-analysis")
text(1:8, par("usr")[3] - 0.2, labels = mylabels, srt = 45, pos = 1, xpd = TRUE,cex=0.8)
abline(h=0)
sample.sizes= table(data1.mat.duohmm2$cohort)
for(i in 1:7){
my.x=i
my.y=summary.duohmm2.mat.results.by.cohort[1,i]
scaling=0.05
offset=0.02*sample.sizes[i]/sample.sizes[2]
segments(x0=i,y0=summary.duohmm2.mat.results.by.cohort[1,i]-1.96*summary.duohmm2.mat.results.by.cohort[2,i],x1=i,y1=summary.duohmm2.mat.results.by.cohort[1,i]+1.96*summary.duohmm2.mat.results.by.cohort[2,i])
polygon(c(my.x-scaling,my.x-scaling,my.x+scaling,my.x+scaling),c(my.y-offset,my.y+offset,my.y+offset,my.y-offset),col=mycols[names(sample.sizes)[i]],border=mycols[names(sample.sizes)[i]])
}
segments(x0=8,y0=meta.fe$ci.lb,x1=8,y1=meta.fe$ci.ub)
dev.off()

#informative_nuclear_families_duoHMM
meta.fe=rma.uni(yi=summary.nftools.mat.results.by.cohort[1,-ncol(summary.nftools.mat.results.by.cohort)],sei=summary.nftools.mat.results.by.cohort[2,-ncol(summary.nftools.mat.results.by.cohort)],method="FE",weighted=T)
pdf("/home/hilary/maternal_age_recombination/paper/figures/box_plot_of_beta_Age_from_LMM_on_informative_nuclear_families_NFTOOLS.pdf",height=4,width=6)
plot(1:8,c(summary.nftools.mat.results.by.cohort[1,-ncol(summary.nftools.mat.results.by.cohort)],meta.fe$b),main="beta_age from LMM",ylim=c(-0.7,1),col=c(rep("white",7),"black"),pch=19,xaxt="n",xlab="Cohort",ylab="beta_age",cex=0.5)
axis(1,at=1:8,labels=FALSE,tick=F,cex.axis=0.5)
mylabels=c(gsub("nftools.","",colnames(summary.nftools.mat.results.by.cohort)[-ncol(summary.nftools.mat.results.by.cohort)]),"meta-analysis")
text(1:8, par("usr")[3] - 0.2, labels = mylabels, srt = 45, pos = 1, xpd = TRUE,cex=0.8)
abline(h=0)
sample.sizes= table(data1.mat.nftools$cohort)
colnames(summary.nftools.mat.results.by.cohort)=gsub("nftools.","",colnames(summary.nftools.mat.results.by.cohort))
summary.nftools.mat.results.by.cohort=summary.nftools.mat.results.by.cohort[,c("FC","GPC","NTR","ORCADES","QTR370","QTR610","VB","all")]
samples.sizes=sample.sizes[colnames(summary.nftools.mat.results.by.cohort)[1:7]]

for(i in 1:7){
my.x=i
my.y=summary.nftools.mat.results.by.cohort[1,i]
scaling=0.05
offset=0.02*sample.sizes[i]/sample.sizes[2]
segments(x0=i,y0=summary.nftools.mat.results.by.cohort[1,i]-1.96*summary.nftools.mat.results.by.cohort[2,i],x1=i,y1=summary.nftools.mat.results.by.cohort[1,i]+1.96*summary.nftools.mat.results.by.cohort[2,i])
polygon(c(my.x-scaling,my.x-scaling,my.x+scaling,my.x+scaling),c(my.y-offset,my.y+offset,my.y+offset,my.y-offset),col=mycols[names(sample.sizes)[i]],border=mycols[names(sample.sizes)[i]])
}
segments(x0=8,y0=meta.fe$ci.lb,x1=8,y1=meta.fe$ci.ub)
dev.off()


####for all informative duos
meta.fe.all=rma.uni(yi=summary.duohmm.mat.results.by.cohort[1,1:9],sei=summary.duohmm.mat.results.by.cohort[2,1:9],method="FE",weighted=T)
pdf("/home/hilary/maternal_age_recombination/paper/figures//box_plot_of_beta_Age_from_LMM_on_all_informative_duos_duoHMM.pdf",height=4,width=6)
plot(1:10,c(summary.duohmm.mat.results.by.cohort[1,1:9],meta.fe.all$b),main="beta_age from LMM",ylim=c(-0.7,1),col=c(rep("white",9),"black"),pch=19,xaxt="n",xlab="Cohort",ylab="beta_age",cex=0.5)
axis(1,at=1:10,labels=FALSE,tick=F,cex.axis=0.5,las=2)
mylabels=c(gsub("duohmm.","",colnames(summary.duohmm.mat.results.by.cohort)[1:9]),"meta-analysis")
text(1:10, par("usr")[3] - 0.2, labels = mylabels, srt = 45, pos = 1, xpd = TRUE,cex=0.8)
abline(h=0)
sample.sizes= table(data1.mat.duohmm$cohort)
for(i in 1:9){
my.x=i
my.y=summary.duohmm.mat.results.by.cohort[1,i]
scaling=0.05
offset=0.02*sample.sizes[i]/sample.sizes[2]
segments(x0=i,y0=summary.duohmm.mat.results.by.cohort[1,i]-1.96*summary.duohmm.mat.results.by.cohort[2,i],x1=i,y1=summary.duohmm.mat.results.by.cohort[1,i]+1.96*summary.duohmm.mat.results.by.cohort[2,i])
polygon(c(my.x-scaling,my.x-scaling,my.x+scaling,my.x+scaling),c(my.y-offset,my.y+offset,my.y+offset,my.y-offset),col=mycols[names(sample.sizes)[i]],border=mycols[names(sample.sizes)[i]])
}
segments(x0=10,y0=meta.fe.all$ci.lb,x1=10,y1=meta.fe.all$ci.ub)

dev.off()

##### Scatterplots
load("duoHMM_data_for_RSTAN.more_stringent.RData")
library(nlme)
library(scales)

pdf("/home/hilary/maternal_age_recombination/paper/figures/scatterplots_of_duoHMM_crossovers_vs_age.informative_duos.pdf",height=9,width=9)
cohorts=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB","ORCADES")
par(mfrow=c(3,3))
for(c in 1:length(cohorts)){
    my.lme=summary(lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data1.mat[data1.mat$cohort==cohorts[c],]))
    plot(data1.mat$age.at.birth[data1.mat$cohort!=cohorts[c]],data1.mat$nrec[data1.mat$cohort!=cohorts[c]],main=cohorts[c],xlab="maternal age at birth",ylab="Number of crossovers",col=alpha("grey",0.3),pch=19,
     xlim=range(data1.mat$age.at.birth),ylim=range(data1.mat$nrec),cex=0.5)
    points(data1.mat$age.at.birth[data1.mat$cohort==cohorts[c]],data1.mat$nrec[data1.mat$cohort==cohorts[c]],col=alpha(mycols[cohorts[c]],0.5),pch=19,cex=0.5)
    abline(a=my.lme$coefficients$fixed[1],b=my.lme$coefficients$fixed[2])
    legend("topleft",paste0("slope = ",round(my.lme$coefficients$fixed[2],2)))
    legend("topright",paste0("p = ",round(my.lme$tTable[2,5],3)))
}
dev.off()

cohorts=c("FC", "GPC","NTR","QTR370", "QTR610", "VB" ,"ORCADES")
pdf("/home/hilary/maternal_age_recombination/paper/figures/scatterplots_of_duoHMM_crossovers_vs_age.informative_nuc_fams.pdf",height=9,width=9)
par(mfrow=c(3,3))
for(c in 1:length(cohorts)){
    my.lme=summary(lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data0.mat[data0.mat$cohort==cohorts[c],]))
    plot(data0.mat$age.at.birth[data0.mat$cohort!=cohorts[c]],data0.mat$nrec[data0.mat$cohort!=cohorts[c]],main=cohorts[c],xlab="maternal age at birth",ylab="Number of crossovers",col=alpha("grey",0.4),pch=19,
         xlim=range(data0.mat$age.at.birth),ylim=range(data0.mat$nrec),cex=0.5)
    points(data0.mat$age.at.birth[data0.mat$cohort==cohorts[c]],data0.mat$nrec[data0.mat$cohort==cohorts[c]],col=alpha(mycols[cohorts[c]],0.5),pch=19,cex=0.5)
    abline(a=my.lme$coefficients$fixed[1],b=my.lme$coefficients$fixed[2])
    legend("topleft",paste0("slope = ",round(my.lme$coefficients$fixed[2],2)))
    legend("topright",paste0("p = ",round(my.lme$tTable[2,5],3)))
}
dev.off()


data1.mat$inf.nuc.fam=data1.mat$duo %in% data0.mat$duo
data1.mat$point.type=19
data1.mat$point.type[!data1.mat$inf.nuc.fam]=3
data1.mat$col="black"
data1.mat$col[!data1.mat$inf.nuc.fam]="red"

pdf("/home/hilary/maternal_age_recombination/paper/figures/scatterplots_of_duoHMM_crossovers_vs_age.informative_duos_vs_informative_nuc_fams.pdf",height=12,width=10)
cohorts=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB","ORCADES")
par(mfrow=c(3,3))
for(c in 1:length(cohorts)){
    my.lme=summary(lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data1.mat[data1.mat$cohort==cohorts[c],]))
    if(sum(data1.mat$cohort==cohorts[c] & data1.mat$inf.nuc.fam)>0){
        my.lme2=summary(lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data1.mat[data1.mat$cohort==cohorts[c] & data1.mat$inf.nuc.fam,]))
    }
    
    plot(data1.mat$age.at.birth[data1.mat$cohort!=cohorts[c]],data1.mat$nrec[data1.mat$cohort!=cohorts[c]],main=cohorts[c],xlab="maternal age at birth",ylab="Number of crossovers",col=alpha("grey",0.2),pch=19,
     xlim=range(data1.mat$age.at.birth),ylim=c(20,80),cex=0.5,cex.main=2,cex.axis=1.5,cex.lab=1.5)
    points(data1.mat$age.at.birth[data1.mat$cohort==cohorts[c]],data1.mat$nrec[data1.mat$cohort==cohorts[c]],pch=data1.mat$point.type[data1.mat$cohort==cohorts[c]],cex=0.5,col=data1.mat$col[data1.mat$cohort==cohorts[c]])
    abline(a=my.lme$coefficients$fixed[1],b=my.lme$coefficients$fixed[2])
    if(sum(data1.mat$cohort==cohorts[c] & data1.mat$inf.nuc.fam)>0){
        abline(a=my.lme2$coefficients$fixed[1],b=my.lme2$coefficients$fixed[2],lty=2)
    }
    legend("topleft",c(paste0("slope_all = ",round(my.lme$coefficients$fixed[2],2)),paste0("p_all = ",round(my.lme$tTable[2,5],3))))
    if(sum(data1.mat$cohort==cohorts[c] & data1.mat$inf.nuc.fam)>0){
        legend("topright",c(paste0("slope_nuc = ",round(my.lme2$coefficients$fixed[2],2)),paste0("p_nuc = ",round(my.lme2$tTable[2,5],3))))
    }

}
dev.off()



#### Compare NFTOOLS to duoHMM counts
rm(list=ls())
load("duoHMM_data_for_RSTAN.more_stringent.RData")
duohmm.mat=data0.mat

load("NFTOOLS_data_for_RSTAN.more_stringent.RData")
nftools.mat = data1.mat

both=cbind(nftools.mat,duohmm.mat[match(as.character(nftools.mat$child),as.character(duohmm.mat$CHILD)),])
both=both[!is.na(both$CHILD),]

both$diff=both$Freq-both$nrec
both$col=NA
both$col=mycols[both$cohort]

pdf("/home/hilary/maternal_age_recombination/paper/figures/NFTOOLS_vs_duoHMM_counts.by_cohort.inform_nuc_fams_only.pdf",height=9,width=9)
par(mfrow=c(3,3))
cohorts=unique(both$cohort)
for(i in 1:length(cohorts)){
    plot(x=jitter(both$Freq[both$cohort==cohorts[i]]),y=jitter(both$nrec[both$cohort==cohorts[i]]),col=alpha(both$col[both$cohort==cohorts[i]],0.5),xlab="NFTOOLS",ylab="duoHMM",pch=19,main=cohorts[i])
    abline(a=0,b=1)
legend("topleft",paste0("r = ",round(cor(both$Freq[both$cohort==cohorts[i]],both$nrec[both$cohort==cohorts[i]]),digits=2)))
}
dev.off()

###make qq plots of duoHMM,

qqr <- function(k1,rate,col1,makeplot=TRUE,title,limits,xlowerlimit,...) {
      k1 <- sort(k1)
        n <- length(k1)
        e <- qpois(1:n/(n+1),rate)
        l <- c(min(k1),min(max(k1),rate+6*sqrt(rate)))
        etmp <- e[k1>=l[2]]
        ret <- cbind(jitter(etmp),l[2])

        if(makeplot) {
            #       plot(e,k1,type='n',xlab="",ylab="",ylim=limits,xlim=c(xlowerlimit,max(limits)),main=title,...)
                         plot(e,k1,type='n',xlab="",ylab="",ylim=limits,xlim=limits,main=title,cex=0.5,...)
                         abline(0,1,col=1,lty=2);grid()
                               points(jitter(e[k1<l[2]]),jitter((k1[k1<l[2]])),col=col1,pch=16,cex=0.5)
                         #      points(ret,col=col1,pch=17)
                     }

        return(  cbind(e,k1))
  }


male.lengths <- c(1.9512,1.8955,1.6071,1.4654,1.512,1.3762,1.2835,1.0794,1.1725,1.3389,1.0936,1.3554,1.0131,0.9462,1.0257,1.081,1.0856,0.9862,0.9264,0.7472,0.4731,0.4896)
female.lengths <- c(3.4541,3.2541,2.7564,2.5906,2.6019,2.4159,2.3033,2.0994,1.982,2.1813,1.9553,2.0664,1.5588,1.4236,1.5496,1.4962,1.6153,1.4257,1.2682,1.2197,0.764,0.8276)
malesize <- sum(male.lengths)
femalesize <- sum(female.lengths)
cols=c("red","darkgreen","purple","blue","pink")

cohort.list=read.delim("cohort_list.txt",header=F,stringsAsFactors=F)
pdf("/home/hilary/maternal_age_recombination/paper/figures/qqplots_duoHMM.informative_duos_only.showing_filtering.pdf",height=16,width=12)
par(mfrow=c(6,4),mar=c(2,2,1,1),cex=1,oma=c(2,2,0,0))
plot(1,1,col="white",xaxt="n",yaxt="n", bty='n')
            legend("topleft",c("Unfiltered","Filter 1: double xovers < X SNPs apart filtered","Filter 1 + clustered duplicate positions removed","Filter 3: Duplicate positions removed",
                               "Filter4: Duplicate positions within families removed"),col=cols[c(1,2,5,3,4)],pch=16,  bg="white",cex=0.5)
plot(1,1,col="white",xaxt="n",yaxt="n", bty='n')
for(i in 1:nrow(cohort.list)){
    mydir=cohort.list[i,1]
    cohort=cohort.list[i,2]
    prefix=cohort.list[i,3]
    argv=c(paste0(mydir,"/duoHMM_results/",prefix,".more_stringent.generr_removed.recombinations.txt"),"","","",cohort)

    rec=read.delim(gsub("recombinations","recombination_count.informative_duos_only.p_0.5.annotated",gsub("1.","all.",argv[1],fixed=T)),header=T)
    rec.flt1=read.delim(gsub("recombinations","recombination_count.informative_duos_only.p_0.5.annotated.double_xovers_within_X_SNPs_removed",gsub("1.","all.",argv[1],fixed=T)),header=T)
    rec.flt2=read.delim(gsub("recombinations","recombination_count.informative_duos_only.p_0.5.annotated.duplicate_xovers_removed",gsub("1.","all.",argv[1],fixed=T)),header=T)
    rec.flt3=read.delim(gsub("recombinations","recombination_count.informative_duos_only.p_0.5.annotated.duplicate_xovers_within_families_removed",gsub("1.","all.",argv[1],fixed=T)),header=T)
    rec.flt5=read.delim(gsub("recombinations","recombination_count.informative_duos_only.p_0.5.annotated.double_xovers_within_X_SNPs_and_clustered_duplicate_xovers_removed",gsub("1.","all.",argv[1],fixed=T)),header=T)
    
    if(nrow(rec)>0){
        print("QQ-PLOT")
#        par(mar=c(2,2,1,1),mfrow=c(1,2),cex=1,oma=c(2,2,0,0))
        tmp <- qqr(subset(rec,sex=="Male")$nrec,malesize,cols[1],title=paste(argv[5]," paternal"),
                                        #           limits=range(c(subset(rec,sex=="Male")$nrec,subset(rec.flt,sex=="Male")$nrec,subset(rec.flt1,sex=="Male")$nrec,subset(rec.flt2,sex=="Male")$nrec)),xlowerlimit=10)
                   limits=range(9,45),xlowerlimit=10)
        points(jitter(qqr(subset(rec.flt2,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[3],cex=0.5)
        points(jitter(qqr(subset(rec.flt3,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[4],cex=0.5)
        points(jitter(qqr(subset(rec.flt5,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[5],cex=0.5)
        points(jitter(qqr(subset(rec.flt1,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[2],cex=0.5)
        if(i==0){
            legend("topleft",c("Unfiltered","Filter 1: double xovers < X SNPs apart filtered","Filter 1 + clustered duplicate positions removed","Filter 3: Duplicate positions removed",
                               "Filter4: Duplicate positions within families removed"),col=cols[c(1,2,5,3,4)],pch=16,  bg="white",cex=0.5)
        }
        tmp <- qqr(subset(rec,sex=="Female")$nrec,femalesize,cols[1],title=paste(argv[5]," maternal"),
                   limits=range(10,76),xlowerlimit=20)
        points(jitter(qqr(subset(rec.flt2,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[3],cex=0.5)
        points(jitter(qqr(subset(rec.flt3,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[4],cex=0.5)
        points(jitter(qqr(subset(rec.flt5,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[5],cex=0.5)
        points(jitter(qqr(subset(rec.flt1,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[2],cex=0.5)
    }
}
    mtext("Expected recombinations genome-wide",1,outer=T,cex=1.6,padj=0.8)
    mtext("Observed recombinations genome-wide",2,outer=T,cex=1.6,padj=-0.5)

dev.off()

    
pdf("/home/hilary/maternal_age_recombination/paper/figures/qqplots_duoHMM.double_xovers_within_X_SNPs_removed.showing_effect_of_family_type.pdf",height=16,width=12)
par(mfrow=c(6,4),mar=c(2,2,1,1),cex=1,oma=c(2,2,0,0))
mycols=c("darkgreen","green","red3","red","purple","purple3","orange","orange2")
mypch=c(21,24,21,24,21,24,21,24)
plot(1,1,col="white",xaxt="n",yaxt="n", bty='n')
legend("topleft",c("3 kids, both parents","3 kids, one parent","2 kids, both parents","2 kids, one parent","1 kid, both parents","1 kid, one parent","3 generations, both parents","3 generations, one parent"),
                  pch=mypch,col=mycols,cex=0.8,pt.bg=mycols)
plot(1,1,col="white",xaxt="n",yaxt="n", bty='n')
for(i in 1:nrow(cohort.list)){
    mydir=cohort.list[i,1]
    cohort=cohort.list[i,2]
    prefix=cohort.list[i,3]
    print(prefix)
    argv=c(paste0(mydir,"/duoHMM_results/",prefix,".more_stringent.generr_removed.recombinations.txt"),"",paste0("qqplots_duoHMM_v3/",cohort,".more_stringent"),"",cohort)

    load(paste0(argv[3],".multiple_filters.RData"))
    tmp <- qqr(subset(rec.flt1,sex=="Male")$nrec,malesize,"white",title=paste(argv[5]," paternal"),limits=range(0,45),xlowerlimit=10)
    if(nrow(subset(rec.flt1.informative.2gen.2parents,sex="Male"))>0){points(jitter(qqr(subset(rec.flt1.informative.2gen.2parents,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[1],col=mycols[1],bg=mycols[1],cex=0.5)}
    if(nrow(subset(rec.flt1.informative.2gen.1parent,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.informative.2gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[2],col=mycols[2],bg=mycols[2],cex=0.5)}
    if(nrow(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[3],col=mycols[3],bg=mycols[3],cex=0.5)}
    if(nrow(subset(rec.flt1.noninf.2kids.2gen.1parent,sex="Male"))>0){points(jitter(qqr(subset(rec.flt1.noninf.2kids.2gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[4],col=mycols[4],bg=mycols[4],cex=0.5)}
    if(nrow(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[5],col=mycols[5],bg=mycols[5],cex=0.5)}
    if(nrow(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[6],col=mycols[6],bg=mycols[6],cex=0.5)}
    if(nrow(subset(rec.flt1.informative.3gen.2parents,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.2parents,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[7],col=mycols[7],bg=mycols[7],cex=0.5)}
    if(nrow(subset(rec.flt1.informative.3gen.1parent,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[8],col=mycols[8],bg=mycols[8],cex=0.5)}
                                        #legend("topleft",c("Unfiltered",paste0("<",thr/1000,"kb filtered")," < X SNPs apart filtered","Duplicate positions removed"),col=cols[c(1:2,4,3)],pch=16,bg="white")
    #tmp <- qqr(subset(rec.flt1.all,sex=="Female")$nrec,femalesize,"white",title=paste(argv[5]," maternal"),limits=range(subset(rec.flt1.all,sex=="Female")$nrec),xlowerlimit=20)
    tmp <- qqr(subset(rec.flt1,sex=="Female")$nrec,femalesize,"white",title=paste(argv[5]," maternal"),limits=range(0,76),xlowerlimit=20)
    if(nrow(subset(rec.flt1.informative.2gen.2parents,sex="Female"))>0){points(jitter(qqr(subset(rec.flt1.informative.2gen.2parents,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[1],col=mycols[1],bg=mycols[1],cex=0.5)}
    if(nrow(subset(rec.flt1.informative.2gen.1parent,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.informative.2gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[2],col=mycols[2],bg=mycols[2],cex=0.5)}
    if(nrow(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[3],col=mycols[3],bg=mycols[3],cex=0.5)}
    if(nrow(subset(rec.flt1.noninf.2kids.2gen.1parent,sex="Female"))>0){points(jitter(qqr(subset(rec.flt1.noninf.2kids.2gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[4],col=mycols[4],bg=mycols[4],cex=0.5)}
    if(nrow(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[5],col=mycols[5],bg=mycols[5],cex=0.5)}
    if(nrow(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[6],col=mycols[6],bg=mycols[6],cex=0.5)}
    if(nrow(subset(rec.flt1.informative.3gen.2parents,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.2parents,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[7],col=mycols[7],bg=mycols[7],cex=0.5)}
    if(nrow(subset(rec.flt1.informative.3gen.1parent,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[8],col=mycols[8],bg=mycols[8],cex=0.5)}
}
       mtext("Expected recombinations genome-wide",1,outer=T,cex=1.6,padj=0.8)
    mtext("Observed recombinations genome-wide",2,outer=T,cex=1.6,padj=-0.5)


dev.off()

#### QQplots for NFTOOLS

load("NFTOOLS_results/NFTOOLS_counts.all_families_except_outliers.more_stringent.RData")
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
male.lengths <- c(1.9512,1.8955,1.6071,1.4654,1.512,1.3762,1.2835,1.0794,1.1725,1.3389,1.0936,1.3554,1.0131,0.9462,1.0257,1.081,1.0856,0.9862,0.9264,0.7472,0.4731,0.4896)
female.lengths <- c(3.4541,3.2541,2.7564,2.5906,2.6019,2.4159,2.3033,2.0994,1.982,2.1813,1.9553,2.0664,1.5588,1.4236,1.5496,1.4962,1.6153,1.4257,1.2682,1.2197,0.764,0.8276)
malesize <- sum(male.lengths)
femalesize <- sum(female.lengths)
pdf(paste0("/home/hilary/maternal_age_recombination/paper/figures/NFTOOLS-qq.pdf"),width=12,height=12)
par(mfrow=c(4,4),mar=c(2,2,1,1),cex=1,oma=c(2,2,0,0))
plot(1,1,col="white",xaxt="n",yaxt="n", bty='n')
    legend("topleft",c("phase known","phase unknown, averaged"),pch=16,col=c("red","darkgreen"),cex=0.9)
plot(1,1,col="white",xaxt="n",yaxt="n", bty='n')

for(i in 1:length(all.events)){
    rec=all.events[[i]]
    rec.2kids=rec[rec$nkids==1,]
    rec.2kids$Freq=rec.2kids$Freq/2
    rec=rec[rec$nkids>2,]
    tmp <- qqr(subset(rec,sex=="male")$Freq,malesize,"red",title=paste(names(all.events)[i]," paternal"),
               limits=c(9,45),xlowerlimit=10)
    points(jitter(qqr(subset(rec.2kids,sex=="male")$Freq,malesize,makeplot=0)),pch=16,col="darkgreen")

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
}
    mtext("Expected recombinations genome-wide",1,outer=T,cex=1.6,padj=0.8)
    mtext("Observed recombinations genome-wide",2,outer=T,cex=1.6,padj=-0.5)

dev.off()
