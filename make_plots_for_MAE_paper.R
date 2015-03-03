if(FALSE){
  outdir="/home/hilary/maternal_age_recombination/paper/figures/"
#load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.more_stringent.RData")
  load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.more_stringent.no_GPC.RData")
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
mat.table = table(data1.mat$binned.age)
pat.table = table(data1.pat$binned.age)

#pdf(paste0(outdir,"/crossover_count_by_binned_age.duoHMM_informative_duos_only.pdf"),height=5,width=5)
#pdf(paste0(outdir,"/crossover_count_by_binned_age.duoHMM_informative_nuc_fams_only.pdf"),height=5,width=5,useDingbats = FALSE )
  pdf(paste0(outdir,"/crossover_count_by_binned_age.duoHMM_informative_nuc_fams_only.no_GPC.pdf"),height=5,width=5,useDingbats = FALSE )
#pdf(paste0(outdir,"/crossover_count_by_binned_age.NFTOOLS_informative_nuc_fams_only.pdf"),height=5,width=5)
par(oma=c(0,0,1,0),mar=c(5,4,2,1))
plot(mat.binned.rec$x,mat.binned.rec$mean,xlab="Parental age, binned",ylab="Difference from 20-25 age group",col="red",ylim=c(-5,6),pch=19)
points(pat.binned.rec$x,pat.binned.rec$mean,col="blue",pch=19)
apply(mat.binned.rec,1,function(bounds){segments(x0=bounds[4],y0=bounds[2],x1=bounds[4],y1=bounds[3],col="red")})
apply(pat.binned.rec,1,function(bounds){segments(x0=bounds[4],y0=bounds[2],x1=bounds[4],y1=bounds[3],col="blue")})
text(mat.binned.rec$x,rep(6,6),labels = mat.table,col="red")
text(mat.binned.rec$x,rep(5.3,6),labels = pat.table,col="blue")

legend("bottom",c("maternal","paternal"),col=c("red","blue"),pch=19,lty=1)
    mtext("Figure 1",3,outer=T,cex=1.2,padj=0.3,at=0.1,font=2)
abline(h=0)
dev.off()

##### Meta-analysis
library(nlme)
library(metafor)
#load("frequentist_analysis/more_stringent/frequentist_analysis_of_MAE.RData")
#load("frequentist_analysis/more_stringent/frequentist_analysis_of_MAE.no_nkids_term.RData")
#load("frequentist_analysis/frequentist_analysis_of_MAE.no_nkids_term.with_new_QTR.RData")

mycols=c("black","blue","red","green","orange","red4","purple","darkgreen","violetred2")
names(mycols)=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB","ORCADES")

colnames(summary.duohmm.mat.results.by.cohort) = gsub("duohmm.","",colnames(summary.duohmm.mat.results.by.cohort))
colnames(summary.duohmm2.mat.results.by.cohort) = gsub("duohmm.in.nftools.","",colnames(summary.duohmm2.mat.results.by.cohort) )
colnames(summary.nftools.mat.results.by.cohort) = gsub("nftools.","",colnames(summary.nftools.mat.results.by.cohort))

colnames(summary.duohmm.pat.results.by.cohort) = gsub("duohmm.","",colnames(summary.duohmm.pat.results.by.cohort))
colnames(summary.duohmm2.pat.results.by.cohort) = gsub("duohmm.in.nftools.","",colnames(summary.duohmm2.pat.results.by.cohort) )
colnames(summary.nftools.pat.results.by.cohort) = gsub("nftools.","",colnames(summary.nftools.pat.results.by.cohort))


nftools=data.frame(t(summary.nftools.mat.results.by.cohort[c(1,2,5),]))
duohmm.all =data.frame(t(summary.duohmm.mat.results.by.cohort[c(1,2,5),]))
duohmm.nuc =data.frame(t(summary.duohmm2.mat.results.by.cohort[c(1,2,5),]))


nftools.pat=data.frame(t(summary.nftools.pat.results.by.cohort[c(1,2,5),]))
duohmm.all.pat =data.frame(t(summary.duohmm.pat.results.by.cohort[c(1,2,5),]))
duohmm.nuc.pat =data.frame(t(summary.duohmm2.pat.results.by.cohort[c(1,2,5),]))

#need to do this for paternal to, and redo meta-analysis excluding GPC

combined = cbind(duohmm.all,duohmm.nuc[rownames(duohmm.all),],nftools[rownames(duohmm.all),],
  duohmm.all.pat[rownames(duohmm.all),],duohmm.nuc.pat[rownames(duohmm.all),],nftools.pat[rownames(duohmm.all),])


colnames(combined) = c(paste0("mat.",c("beta_Age_duohmm_all","se_duohmm_all","p_duohmm_all", "beta_Age_duohmm_nuc","se_duohmm_nuc","p_duohmm_nuc","beta_Age_nftools","se_nftools","p_nftools")),
          paste0("pat.",c("beta_Age_duohmm_all","se_duohmm_all","p_duohmm_all", "beta_Age_duohmm_nuc","se_duohmm_nuc","p_duohmm_nuc","beta_Age_nftools","se_nftools","p_nftools")))

write.table(combined,"/home/hilary/maternal_age_recombination/paper/summary_of_LMM_results.txt",quote=F,sep="\t")



summary.duohmm2.mat.results.by.cohort.no.gpc   = summary.duohmm2.mat.results.by.cohort[,-2]
summary.nftools.mat.results.by.cohort.no.gpc= summary.nftools.mat.results.by.cohort[,-5]
summary.duohmm.mat.results.by.cohort.no.gpc = summary.duohmm.mat.results.by.cohort[,-4]
  
meta.fe.duohmm.i.n.f=rma.uni(yi=summary.duohmm2.mat.results.by.cohort.no.gpc[1,-ncol(summary.duohmm2.mat.results.by.cohort.no.gpc)],sei=summary.duohmm2.mat.results.by.cohort.no.gpc[2,-ncol(summary.duohmm2.mat.results.by.cohort.no.gpc)],method="FE",weighted=T)
meta.fe.nftools=rma.uni(yi=summary.nftools.mat.results.by.cohort.no.gpc[1,-ncol(summary.nftools.mat.results.by.cohort.no.gpc)],sei=summary.nftools.mat.results.by.cohort.no.gpc[2,-ncol(summary.nftools.mat.results.by.cohort.no.gpc)],method="FE",weighted=T)
meta.fe.duohmm=rma.uni(yi=summary.duohmm.mat.results.by.cohort.no.gpc[1,1:8],sei=summary.duohmm.mat.results.by.cohort.no.gpc[2,1:8],method="FE",weighted=T)

meta.fe.duohmm.i.n.f.with.gpc=rma.uni(yi=summary.duohmm2.mat.results.by.cohort[1,-ncol(summary.duohmm2.mat.results.by.cohort)],sei=summary.duohmm2.mat.results.by.cohort[2,-ncol(summary.duohmm2.mat.results.by.cohort)],method="FE",weighted=T)
meta.fe.nftools.with.gpc=rma.uni(yi=summary.nftools.mat.results.by.cohort[1,-ncol(summary.nftools.mat.results.by.cohort)],sei=summary.nftools.mat.results.by.cohort[2,-ncol(summary.nftools.mat.results.by.cohort)],method="FE",weighted=T)
meta.fe.duohmm.with.gpc=rma.uni(yi=summary.duohmm.mat.results.by.cohort[1,1:9],sei=summary.duohmm.mat.results.by.cohort[2,1:9],method="FE",weighted=T)



summary.duohmm2.pat.results.by.cohort.no.gpc   = summary.duohmm2.pat.results.by.cohort[,-2]
summary.nftools.pat.results.by.cohort.no.gpc= summary.nftools.pat.results.by.cohort[,-5]
summary.duohmm.pat.results.by.cohort.no.gpc = summary.duohmm.pat.results.by.cohort[,-3]
  
pat.meta.fe.duohmm.i.n.f=rma.uni(yi=summary.duohmm2.pat.results.by.cohort.no.gpc[1,-ncol(summary.duohmm2.pat.results.by.cohort.no.gpc)],sei=summary.duohmm2.pat.results.by.cohort.no.gpc[2,-ncol(summary.duohmm2.pat.results.by.cohort.no.gpc)],method="FE",weighted=T)
pat.meta.fe.nftools=rma.uni(yi=summary.nftools.pat.results.by.cohort.no.gpc[1,-ncol(summary.nftools.pat.results.by.cohort.no.gpc)],sei=summary.nftools.pat.results.by.cohort.no.gpc[2,-ncol(summary.nftools.pat.results.by.cohort.no.gpc)],method="FE",weighted=T)
pat.meta.fe.duohmm=rma.uni(yi=summary.duohmm.pat.results.by.cohort.no.gpc[1,1:8],sei=summary.duohmm.pat.results.by.cohort.no.gpc[2,1:8],method="FE",weighted=T)


pat.meta.fe.duohmm.i.n.f.with.gpc=rma.uni(yi=summary.duohmm2.pat.results.by.cohort[1,-ncol(summary.duohmm2.pat.results.by.cohort)],sei=summary.duohmm2.pat.results.by.cohort[2,-ncol(summary.duohmm2.pat.results.by.cohort)],method="FE",weighted=T)
pat.meta.fe.nftools.with.gpc=rma.uni(yi=summary.nftools.pat.results.by.cohort[1,-ncol(summary.nftools.pat.results.by.cohort)],sei=summary.nftools.pat.results.by.cohort[2,-ncol(summary.nftools.pat.results.by.cohort)],method="FE",weighted=T)
pat.meta.fe.duohmm.with.gpc=rma.uni(yi=summary.duohmm.pat.results.by.cohort[1,1:9],sei=summary.duohmm.pat.results.by.cohort[2,1:9],method="FE",weighted=T)

keep=  c("b","se","zval","pval","ci.lb","ci.ub","QE","QEp")

meta.summary = cbind(unlist(meta.fe.duohmm)[keep],unlist(meta.fe.duohmm.i.n.f)[keep],unlist(meta.fe.nftools)[keep],
      unlist(meta.fe.duohmm.with.gpc)[keep],unlist(meta.fe.duohmm.i.n.f.with.gpc)[keep],unlist(meta.fe.nftools.with.gpc)[keep],
      unlist(pat.meta.fe.duohmm)[keep],unlist(pat.meta.fe.duohmm.i.n.f)[keep],unlist(pat.meta.fe.nftools)[keep],
      unlist(pat.meta.fe.duohmm.with.gpc)[keep],unlist(pat.meta.fe.duohmm.i.n.f.with.gpc)[keep],unlist(pat.meta.fe.nftools.with.gpc)[keep])
colnames(meta.summary) = c(paste0("mat.no.gpc.",c("duohmm.all","duohmm.nuc","nftools")),paste0("mat.with.gpc.",c("duohmm.all","duohmm.nuc","nftools")),
          paste0("pat.no.gpc.",c("duohmm.all","duohmm.nuc","nftools")),paste0("pat.with.gpc.",c("duohmm.all","duohmm.nuc","nftools")))


write.table(meta.summary,"/home/hilary/maternal_age_recombination/paper/summary_of_meta_analysis_results.txt",quote=F,sep="\t")

pdf("/home/hilary/maternal_age_recombination/paper/figures/box_plot_of_beta_Age_from_LMM_on_different_family_types.no_nkids_term.no_GPC_in_metaanalysis.pdf",height=4,width=6,useDingbats=F)
par(oma=c(0,0,1,0),mar=c(5,4,2,1))
plot(1:10,rep(0,10),main="",ylim=c(-0.7,1.2),col=c(rep("white",10)),pch=19,xaxt="n",xlab="",ylab="",cex=0.5,xlim=c(1,10.5))
mtext(expression(beta["age"]),2,line=2.5)
mtext("Cohort",1,line=3.2)

axis(1,at=1:10,labels=FALSE,tick=F,cex.axis=0.5,las=2)
mylabels=c(gsub("duohmm.","",colnames(summary.duohmm.mat.results.by.cohort)[1:9]),"meta-analysis")
text(1:10, par("usr")[3] - 0.17, labels = mylabels, srt = 45, pos = 1, xpd = TRUE,cex=0.8)
legend("topleft",c("NFTOOLS, informative nuclear families","duoHMM, informative nuclear families","duoHMM, all informative meioses"),lty=c(1,2,4),cex=0.8)
    mtext("Figure 2",3,outer=T,cex=1.2,padj=0.3,at=0.1,font=2)
abline(h=0)
xs = 1:9
  names(xs) = colnames(summary.duohmm.mat.results.by.cohort)[1:9]
  ref = 142
  summary.nftools.mat.results.by.cohort=summary.nftools.mat.results.by.cohort[,c("FC","GPC","NTR","ORCADES","QTR370","QTR610","VB","all")]
  sample.sizes=table(data1.mat.nftools$cohort)
print(sample.sizes)
#sample.sizes[colnames(summary.nftools.mat.results.by.cohort)[1:7]]
for(i in 1:7){
my.x=xs[colnames(summary.nftools.mat.results.by.cohort)[i]]
my.y=summary.nftools.mat.results.by.cohort[1,i]
scaling=0.05
offset=0.02*sample.sizes[colnames(summary.nftools.mat.results.by.cohort)[i]]/ref
segments(x0=my.x,y0=summary.nftools.mat.results.by.cohort[1,i]-1.96*summary.nftools.mat.results.by.cohort[2,i],x1=my.x,y1=summary.nftools.mat.results.by.cohort[1,i]+1.96*summary.nftools.mat.results.by.cohort[2,i],lty=1)
polygon(c(my.x-scaling,my.x-scaling,my.x+scaling,my.x+scaling),c(my.y-offset,my.y+offset,my.y+offset,my.y-offset),col=mycols[colnames(summary.nftools.mat.results.by.cohort)[i]],border=mycols[colnames(summary.nftools.mat.results.by.cohort)[i]])
}
points(10,meta.fe.nftools$b,pch=19)
segments(x0=10,y0=meta.fe.nftools$ci.lb,x1=10,y1=meta.fe.nftools$ci.ub)

sample.sizes= table(data1.mat.duohmm2$cohort)
print(sample.sizes)
for(i in 1:7){
my.x=xs[colnames(summary.duohmm2.mat.results.by.cohort)[i]]+0.2
my.y=summary.duohmm2.mat.results.by.cohort[1,i]
scaling=0.05
offset=0.02*sample.sizes[colnames(summary.duohmm2.mat.results.by.cohort)[i]]/ref
segments(x0=my.x,y0=summary.duohmm2.mat.results.by.cohort[1,i]-1.96*summary.duohmm2.mat.results.by.cohort[2,i],x1=my.x,y1=summary.duohmm2.mat.results.by.cohort[1,i]+1.96*summary.duohmm2.mat.results.by.cohort[2,i],lty=2)
polygon(c(my.x-scaling,my.x-scaling,my.x+scaling,my.x+scaling),c(my.y-offset,my.y+offset,my.y+offset,my.y-offset),col=mycols[colnames(summary.duohmm2.mat.results.by.cohort)[i]],border=mycols[colnames(summary.duohmm2.mat.results.by.cohort)[i]])
}
  points(10.2,meta.fe.duohmm.i.n.f$b,pch=19)
segments(x0=10.2,y0=meta.fe.duohmm.i.n.f$ci.lb,x1=10.2,y1=meta.fe.duohmm.i.n.f$ci.ub,lty=2)


sample.sizes= table(data1.mat.duohmm$cohort)
print(sample.sizes)
for(i in 1:9){
my.x=i+0.4
my.y=summary.duohmm.mat.results.by.cohort[1,i]
scaling=0.05
offset=0.02*sample.sizes[colnames(summary.duohmm.mat.results.by.cohort)[i]]/ref
segments(x0=my.x,y0=summary.duohmm.mat.results.by.cohort[1,i]-1.96*summary.duohmm.mat.results.by.cohort[2,i],x1=my.x,y1=summary.duohmm.mat.results.by.cohort[1,i]+1.96*summary.duohmm.mat.results.by.cohort[2,i],lty=4)
polygon(c(my.x-scaling,my.x-scaling,my.x+scaling,my.x+scaling),c(my.y-offset,my.y+offset,my.y+offset,my.y-offset),col=mycols[colnames(summary.duohmm.mat.results.by.cohort)[i]],border=mycols[colnames(summary.duohmm.mat.results.by.cohort)[i]])
}
points(10.4,meta.fe.duohmm$b,pch=19)
segments(x0=10.4,y0=meta.fe.duohmm$ci.lb,x1=10.4,y1=meta.fe.duohmm$ci.ub,lty=4)


dev.off()



#fixed effects meta-analysis of betas
#informative_nuclear_families_duoHMM
meta.fe=rma.uni(yi=summary.duohmm2.mat.results.by.cohort[1,-ncol(summary.duohmm2.mat.results.by.cohort)],sei=summary.duohmm2.mat.results.by.cohort[2,-ncol(summary.duohmm2.mat.results.by.cohort)],method="FE",weighted=T)
pdf("/home/hilary/maternal_age_recombination/paper/figures/box_plot_of_beta_Age_from_LMM_on_informative_nuclear_families_duoHMM.no_nkids_term.pdf",height=4,width=6)
plot(1:8,c(summary.duohmm2.mat.results.by.cohort[1,-ncol(summary.duohmm2.mat.results.by.cohort)],meta.fe$b),main="",ylim=c(-0.7,1),col=c(rep("white",7),"black"),pch=19,xaxt="n",xlab="Cohort",
     ylab="beta_age",cex=0.5)
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

#informative_nuclear_families_NFTOOLS
#meta.fe=rma.uni(yi=summary.duohmm2.mat.results.by.cohort[1,-ncol(summary.duohmm2.mat.results.by.cohort)],sei=summary.duohmm2.mat.results.by.cohort[2,-ncol(summary.duohmm2.mat.results.by.cohort)],method="FE",weighted=T)
meta.fe=rma.uni(yi=summary.nftools.mat.results.by.cohort[1,-ncol(summary.nftools.mat.results.by.cohort)],sei=summary.nftools.mat.results.by.cohort[2,-ncol(summary.nftools.mat.results.by.cohort)],method="FE",weighted=T)
pdf("/home/hilary/maternal_age_recombination/paper/figures/box_plot_of_beta_Age_from_LMM_on_informative_nuclear_families_NFTOOLS.no_nkids_term.pdf",height=4,width=6)
plot(1:8,c(summary.nftools.mat.results.by.cohort[1,-ncol(summary.nftools.mat.results.by.cohort)],meta.fe$b),main="",ylim=c(-0.7,1),col=c(rep("white",7),"black"),pch=19,xaxt="n",xlab="Cohort",ylab="beta_age",cex=0.5)
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
pdf("/home/hilary/maternal_age_recombination/paper/figures//box_plot_of_beta_Age_from_LMM_on_all_informative_duos_duoHMM.no_nkids_term.pdf",height=4,width=6)
plot(1:10,c(summary.duohmm.mat.results.by.cohort[1,1:9],meta.fe.all$b),main="",ylim=c(-0.7,1),col=c(rep("white",9),"black"),pch=19,xaxt="n",xlab="Cohort",ylab="beta_age",cex=0.5)
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
  mycols=c("black","blue","red","green","orange","red4","purple","darkgreen","violetred2")
names(mycols)=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB","ORCADES")

both$diff=both$Freq-both$nrec
both$col=NA
both$col=mycols[both$cohort]

#check outliers
  tail(both[order(abs(both$diff)),])
fc = read.delim("FC/duoHMM_results/FC.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",header=T)
#table(fc[fc$CHILD=="65_556_65" & fc$sex=="Female","chr"])
# 1  2  3  4  5  6  7  8  9 10 11 13 14 15 17 18 19 22
#   1  4  1  3  2  1  1 21  1  3  2  3  1  2  1  2  1  1 
#excess on chromosome 8; no excess of mendel errors or missingness - maybe just due to phasing errors? most are called with 0.5<p<0.8, vs p>0.9 for the other chromosomes
#need to remove this individual
gpc = read.delim("GPC/duoHMM_results/GPC.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",header=T)
table(gpc[gpc$CHILD=="APP5139271" & gpc$sex=="Female","chr"])   
#nothing unusual here
#in NFTOOLS data, high numbers seem to be spread across the chromosome
#not obvious missingness or number of mendelian errors is any different

both = both[both$CHILD != "65_556_65",]

pdf("/home/hilary/maternal_age_recombination/paper/figures/NFTOOLS_vs_duoHMM_counts.by_cohort.inform_nuc_fams_only.pdf",height=9,width=9)
par(mfrow=c(3,3),mar=c(2,2,2,2),cex=1,oma=c(2,2,3,2))
cohorts=unique(both$cohort)
for(i in 1:length(cohorts)){
    plot(x=jitter(both$Freq[both$cohort==cohorts[i]]),y=jitter(both$nrec[both$cohort==cohorts[i]]),col=alpha(both$col[both$cohort==cohorts[i]],0.5),xlab="",ylab="",pch=19,main=cohorts[i],cex.axis=1,cex.main=1.5,cex.lab=1.5)
    abline(a=0,b=1)
legend("topleft",paste0("r = ",round(cor(both$Freq[both$cohort==cohorts[i]],both$nrec[both$cohort==cohorts[i]]),digits=2)),cex=1)
}
    mtext("Supplementary Figure 3",3,outer=T,cex=2,padj=0.2,at=0.2,font=2,line=1.5)
    mtext("number of crossovers called by NFTOOLS",1,outer=T,cex=1.6,padj=0.8)
    mtext("number of crossovers called by duoHMM",2,outer=T,cex=1.6,padj=0.8,line=1.5)
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
pdf("/home/hilary/maternal_age_recombination/paper/figures/FigureS2_qqplots_duoHMM.informative_duos_only.showing_filtering.pdf",height=16,width=12,useDingbats=F)


par(mfrow=c(6,4),mar=c(2,2,1,1),cex=1,oma=c(2,2,3,0))
#plot(1,1,col="white",xaxt="n",yaxt="n", bty='n')
#            legend("topleft",c("Unfiltered","Filter 1: double xovers < X SNPs apart filtered","Filter 1 + clustered duplicate positions removed","Filter 3: Duplicate positions removed",
#                               "Filter4: Duplicate positions within families removed"),col=cols[c(1,2,5,3,4)],pch=16,  bg="white",cex=0.5)
            legend("topleft",c("Unfiltered","double xovers < X SNPs apart filtered"),col=cols[c(1,2)],pch=16,  bg="white",cex=0.5)
#plot(1,1,col="white",xaxt="n",yaxt="n", bty='n')
for(i in 1:nrow(cohort.list)){
    mydir=cohort.list[i,1]
    cohort=cohort.list[i,2]
    prefix=cohort.list[i,3]
    argv=c(paste0(mydir,"/duoHMM_results/",prefix,".more_stringent.generr_removed.recombinations.txt"),"","","",cohort)

    rec=read.delim(gsub("recombinations","recombination_count.informative_duos_only.p_0.5.annotated",gsub("1.","all.",argv[1],fixed=T)),header=T)
    rec.flt1=read.delim(gsub("recombinations","recombination_count.informative_duos_only.p_0.5.annotated.double_xovers_within_X_SNPs_removed",gsub("1.","all.",argv[1],fixed=T)),header=T)
#    rec.flt2=read.delim(gsub("recombinations","recombination_count.informative_duos_only.p_0.5.annotated.duplicate_xovers_removed",gsub("1.","all.",argv[1],fixed=T)),header=T)
#    rec.flt3=read.delim(gsub("recombinations","recombination_count.informative_duos_only.p_0.5.annotated.duplicate_xovers_within_families_removed",gsub("1.","all.",argv[1],fixed=T)),header=T)
#    rec.flt5=read.delim(gsub("recombinations","recombination_count.informative_duos_only.p_0.5.annotated.double_xovers_within_X_SNPs_and_clustered_duplicate_xovers_removed",gsub("1.","all.",argv[1],fixed=T)),header=T)
    
    if(nrow(rec)>0){
        print("QQ-PLOT")
#        par(mar=c(2,2,1,1),mfrow=c(1,2),cex=1,oma=c(2,2,0,0))
        tmp <- qqr(subset(rec,sex=="Male")$nrec,malesize,cols[1],title=paste(argv[5]," paternal"),
                                        #           limits=range(c(subset(rec,sex=="Male")$nrec,subset(rec.flt,sex=="Male")$nrec,subset(rec.flt1,sex=="Male")$nrec,subset(rec.flt2,sex=="Male")$nrec)),xlowerlimit=10)
                   limits=range(9,45),xlowerlimit=10)
  #      points(jitter(qqr(subset(rec.flt2,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[3],cex=0.5)
  #      points(jitter(qqr(subset(rec.flt3,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[4],cex=0.5)
  #      points(jitter(qqr(subset(rec.flt5,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[5],cex=0.5)
        points(jitter(qqr(subset(rec.flt1,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[2],cex=0.5)
        if(i==1){
  #          legend("topleft",c("Unfiltered","Filter 1: double xovers < X SNPs apart filtered","Filter 1 + clustered duplicate positions removed","Filter 3: Duplicate positions removed",
  #                             "Filter4: Duplicate positions within families removed"),col=cols[c(1,2,5,3,4)],pch=16,  bg="white",cex=0.5)
            legend("topleft",c("unfiltered","filtered"),col=cols[c(1,2)],pch=16,  bg="white",cex=0.8)
        }
        tmp <- qqr(subset(rec,sex=="Female")$nrec,femalesize,cols[1],title=paste(argv[5]," maternal"),
                   limits=range(10,76),xlowerlimit=20)
 #       points(jitter(qqr(subset(rec.flt2,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[3],cex=0.5)
 #       points(jitter(qqr(subset(rec.flt3,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[4],cex=0.5)
 #       points(jitter(qqr(subset(rec.flt5,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[5],cex=0.5)
        points(jitter(qqr(subset(rec.flt1,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[2],cex=0.5)
    }
}
    mtext("Expected recombinations genome-wide",1,outer=T,cex=1.6,padj=0.8)
    mtext("Observed recombinations genome-wide",2,outer=T,cex=1.6,padj=-0.5)
    mtext("Supplementary Figure 2",3,outer=T,cex=2,padj=0.3,at=0.2,font=2,line=1.5)
dev.off()

    
pdf("/home/hilary/maternal_age_recombination/paper/figures/qqplots_duoHMM.double_xovers_within_X_SNPs_removed.showing_effect_of_family_type.pdf",height=16,width=12,useDingbats=F)
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

    mtext("Supplementary Figure 8",3,outer=T,cex=2,padj=0.3,at=0.2,font=2,line=1.5)

dev.off()
#par(mfrow=c(6,4),mar=c(2,2,1,1),cex=1,oma=c(2,2,0,0))
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
                         points(jitter(e[k1<l[2]]),jitter((k1[k1<l[2]])),col=alpha(col1,0.5),pch=16)
#                         points(ret,col=col1,pch=17)
               }
        return(  cbind(e,k1))
  }
male.lengths <- c(1.9512,1.8955,1.6071,1.4654,1.512,1.3762,1.2835,1.0794,1.1725,1.3389,1.0936,1.3554,1.0131,0.9462,1.0257,1.081,1.0856,0.9862,0.9264,0.7472,0.4731,0.4896)
female.lengths <- c(3.4541,3.2541,2.7564,2.5906,2.6019,2.4159,2.3033,2.0994,1.982,2.1813,1.9553,2.0664,1.5588,1.4236,1.5496,1.4962,1.6153,1.4257,1.2682,1.2197,0.764,0.8276)
malesize <- sum(male.lengths)
femalesize <- sum(female.lengths)
pdf(paste0("/home/hilary/maternal_age_recombination/paper/figures/NFTOOLS-qq.pdf"),width=12,height=12,useDingbats=F)
par(mfrow=c(4,4),mar=c(2,2,1,1),cex=1,oma=c(2,2,3,0))

#plot(1,1,col="white",xaxt="n",yaxt="n", bty='n')
#    legend("topleft",c("phase known","phase unknown, averaged"),pch=16,col=c("red","darkgreen"),cex=0.9)
#plot(1,1,col="white",xaxt="n",yaxt="n", bty='n')
for(i in 1:length(all.events)){
    rec=all.events[[i]]
    rec.2kids=rec[rec$nkids==1,]
    rec.2kids$Freq=rec.2kids$Freq/2
    rec=rec[rec$nkids>2,]
    tmp <- qqr(subset(rec,sex=="male")$Freq,malesize,"red",title=paste(names(all.events)[i]," paternal"),
               limits=c(9,45),xlowerlimit=10)
#    points(jitter(qqr(subset(rec.2kids,sex=="male")$Freq,malesize,makeplot=0)),pch=16,col="darkgreen")

    cat(names(all.events)[i],"\t",range(subset(rec,sex=="male")$Freq),"\n")
    if(sum(rec$sex=="male" & rec$Freq>50)>0){
       print(rec[rec$sex=="male" & rec$Freq>50,])
    }
    tmp <- qqr(subset(rec,sex=="female")$Freq,femalesize,"red",title=paste(names(all.events)[i]," maternal"),
               limits=c(10,76),xlowerlimit=20)
 #       points(jitter(qqr(subset(rec.2kids,sex=="female")$Freq,femalesize,makeplot=0)),pch=16,col="darkgreen")  
    cat(names(all.events)[i],"\t",range(subset(rec,sex=="female")$Freq),"\n")
    if(sum(rec$sex=="female" & rec$Freq>70)>0){
        print(rec[rec$sex=="female" & rec$Freq>70,])
    }
}
    mtext("Expected recombinations genome-wide",1,outer=T,cex=1.6,padj=0.8)
    mtext("Observed recombinations genome-wide",2,outer=T,cex=1.6,padj=-0.5)
#par(mfrow=c(4,4),mar=c(2,2,1,1),cex=1,oma=c(2,2,3,0))
    mtext("Supplementary Figure 1",3,outer=T,cex=2,padj=0.3,at=0.15,font=2,line=1.5)
dev.off()



#### plot age distribution
load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.more_stringent.RData")
outdir="/home/hilary/maternal_age_recombination/paper/figures/"
mycols=c("black","blue","red","green","orange","red4","purple","darkgreen","violetred2","yellow","skyblue")
names(mycols)=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB","ORCADES","VIS","KORCULA")

pdf(paste0(outdir,"distribution_of_parental_age.pdf"),height=5,width=10)
par(mfrow=c(1,2),mar=c(2,2,2,1),oma=c(2,2,1,0))

for(i in 1:length(mycols)){
  if(i==1){
    plot(density(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"]),col=mycols[i],main="Maternal",xlab="age at birth",ylim=c(0,0.12),lwd=2,xlim=
         range(all.data[all.data$sex=="Female" ,"age.at.birth"]))
} else {
  lines(density(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"]),col=mycols[i],lwd=2)
}
#curve(dnorm(x,mean=mean(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"]),sd=sd(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"])),lty=2,col=mycols[i],
#            add=T )
}
legend("topright",names(mycols),col=mycols,lwd=2,cex=0.8)
for(i in 1:length(mycols)){
  if(i==1){
    plot(density(all.data[all.data$sex=="Male" & all.data$cohort==names(mycols)[i],"age.at.birth"]),col=mycols[i],main="Paternal",xlab="age at birth" ,ylim=c(0,0.12),lwd=2,
         xlim=range(all.data[all.data$sex=="Male" ,"age.at.birth"]))
} else {
  lines(density(all.data[all.data$sex=="Male" & all.data$cohort==names(mycols)[i],"age.at.birth"]),col=mycols[i],lwd=2)
}
}
mtext("Figure S16",3,outer=T,cex=1.5,padj=0.3,at=0.05,font=2)
dev.off()



pdf(paste0(outdir,"distribution_of_parental_age.pdf"),height=15,width=15)
par(mfrow=c(6,4),mar=c(2,2,2,1),oma=c(2,2,1,0))

for(i in 1:length(mycols)){
#x=      hist(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"],main=paste0(names(mycols)[i]," maternal"),xlab="age at birth",lwd=2,freq=F ,
#             ylim=c(0,0.12),xlim=range(all.data[all.data$sex=="Female","age.at.birth"]),breaks=seq(from=12,to=48,by=1),plot=F)

print(x$breaks)
#    plot(density(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"]),col=mycols[i],main="Maternal",xlab="age at birth",ylim=c(0,0.1),lwd=2)
      hist(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"],main=paste0(names(mycols)[i]," maternal"),xlab="age at birth",lwd=2,freq=F ,
           xlim=range(all.data[all.data$sex=="Female","age.at.birth"]),breaks=20,cex.main=2,cex.lab=2,cex.axis=2)
#,           ylim=c(0,0.12)
curve(dnorm(x,mean=mean(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"]),sd=sd(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"])),lty=2,col=mycols[i],
            add=T ,lwd=3)

#    plot(density(all.data[all.data$sex=="Male" & all.data$cohort==names(mycols)[i],"age.at.birth"]),col=mycols[i],main="Paternal",xlab="age at birth" ,ylim=c(0,0.1),lwd=2)
      hist(all.data[all.data$sex=="Male" & all.data$cohort==names(mycols)[i],"age.at.birth"],main=paste0(names(mycols)[i]," paternal"),xlab="age at birth",lwd=2,freq=F ,
           xlim=range(all.data[all.data$sex=="Male","age.at.birth"]),breaks=20,cex.main=2,cex.lab=2,cex.axis=2)
#ylim=c(0,0.14),
#seq(from=12,to=62,by=1))
curve(dnorm(x,mean=mean(all.data[all.data$sex=="Male" & all.data$cohort==names(mycols)[i],"age.at.birth"]),sd=sd(all.data[all.data$sex=="Male" & all.data$cohort==names(mycols)[i],"age.at.birth"])),lty=2,col=mycols[i],
            add=T ,lwd=3)

#  lines(density(all.data[all.data$sex=="Male" & all.data$cohort==names(mycols)[i],"age.at.birth"]),col=mycols[i],lwd=2)

}
dev.off()


pdf(paste0(outdir,"distribution_of_parental_age.excluding_one_of_each_twin_pair.pdf"),height=15,width=15)
par(mfrow=c(6,4),mar=c(2,2,2,1),oma=c(2,2,1,0))
all.data.wo.twins=unique(all.data[,c("PARENT","sex","age.at.birth","cohort")])
all.data.original=all.data
all.data=all.data.wo.twins
for(i in 1:length(mycols)){

  hist(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"],main=paste0(names(mycols)[i]," maternal"),xlab="age at birth",lwd=2,freq=F ,
           xlim=range(all.data[all.data$sex=="Female","age.at.birth"]),breaks=20,cex.main=2,cex.lab=2,cex.axis=2)
curve(dnorm(x,mean=mean(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"]),sd=sd(all.data[all.data$sex=="Female" & all.data$cohort==names(mycols)[i],"age.at.birth"])),lty=2,col=mycols[i],
            add=T ,lwd=3)
      hist(all.data[all.data$sex=="Male" & all.data$cohort==names(mycols)[i],"age.at.birth"],main=paste0(names(mycols)[i]," paternal"),xlab="age at birth",lwd=2,freq=F ,
           xlim=range(all.data[all.data$sex=="Male","age.at.birth"]),breaks=20,cex.main=2,cex.lab=2,cex.axis=2)
curve(dnorm(x,mean=mean(all.data[all.data$sex=="Male" & all.data$cohort==names(mycols)[i],"age.at.birth"]),sd=sd(all.data[all.data$sex=="Male" & all.data$cohort==names(mycols)[i],"age.at.birth"])),lty=2,col=mycols[i],
            add=T ,lwd=3)
}
dev.off()


##### qq plots for informative meioses using different p cutoffs in duoHMM
}
qqr <- function(k1,rate,col1,makeplot=TRUE,title,limits,xlowerlimit,mypch=16,...) {
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
                               points(jitter(e[k1<l[2]]),jitter((k1[k1<l[2]])),col=col1,pch=mypch)
                         #      points(ret,col=col1,pch=17)
                     }

        return(  cbind(e,k1))
  }


colours=c("black","blue","red","green","darkgreen","skyblue","red4","purple","orange","pink3")
male.lengths <- c(1.9512,1.8955,1.6071,1.4654,1.512,1.3762,1.2835,1.0794,1.1725,1.3389,1.0936,1.3554,1.0131,0.9462,1.0257,1.081,1.0856,0.9862,0.9264,0.7472,0.4731,0.4896)
female.lengths <- c(3.4541,3.2541,2.7564,2.5906,2.6019,2.4159,2.3033,2.0994,1.982,2.1813,1.9553,2.0664,1.5588,1.4236,1.5496,1.4962,1.6153,1.4257,1.2682,1.2197,0.764,0.8276)
malesize <- sum(male.lengths)
femalesize <- sum(female.lengths)

cohort.list=read.delim("cohort_list.txt",header=F,stringsAsFactors=F)
pdf("qqplots_duoHMM_v3/all_cohorts..more_stringent.qqplots_showing_different_p_thresholds.pdf",height=30,width=15)
par(mar=c(2,2,2,1),mfrow=c(11,4),cex=1,oma=c(2,2,2,0))

for(c in 1:nrow(cohort.list)){
  cohort=cohort.list[c,2]
  maternal=read.delim(paste0("qqplots_duoHMM_v3/",cohort,".more_stringent.informative_meioses..maternal_recombination_counts_different_p_thresholds.txt"),header=T,stringsAsFactors=F)
paternal=read.delim(paste0("qqplots_duoHMM_v3/",cohort,".more_stringent.informative_meioses..paternal_recombination_counts_different_p_thresholds.txt"),header=T,stringsAsFactors=F)
cutoffs=c(0,0.3,0.5,0.7,0.9)
for(i in 1:length(cutoffs)){
    p=cutoffs[i]
    maternal.counts=maternal[maternal$informative,paste("E.nrec.p.",p,sep="")]
    if(p==0){
        tmp <- qqr(maternal.counts,femalesize,colours[i],title=paste0(cohort, " maternal, expected counts"),limits=range(20,76),xlowerlimit=20,mypch=16)
    } else{
        points(jitter(qqr(maternal.counts,femalesize,makeplot=0)),pch=16,col=colours[i])
    }
}
  if(c==1){
    legend("topleft",paste0("Expected number for p > ",cutoffs),col=colours,pch=16,  bg="white",cex=0.7)
}
mtext("Figure S17",3,outer=T,cex=2,padj=0.3,at=0.05)
                                        #mtext("Expected recombinations genome-wide",1,outer=T,cex=1.3,padj=1.2)
#mtext("Observed recombinations genome-wide",2,outer=T,cex=1.3,padj=-0.5)
  for(i in 1:length(cutoffs)){
    p=cutoffs[i]
    maternal.counts=maternal[maternal$informative,paste("nrec.p.",p,sep="")]
    if(p==0){
        tmp <- qqr(maternal.counts,femalesize,colours[i],title=paste0(cohort, " maternal, hard cutoff"),limits=range(20,76),xlowerlimit=20,mypch=17)
    } else{
        points(jitter(qqr(maternal.counts,femalesize,makeplot=0)),col=colours[i],pch=17)
    }
}
if(c==1){
  legend("topleft",paste0("Number with p > ",cutoffs),col=colours,pch=17,  bg="white",cex=0.7)
}
  for(i in 1:length(cutoffs)){
    p=cutoffs[i]
    paternal.counts=paternal[paternal$informative,paste("E.nrec.p.",p,sep="")]
    if(p==0){
        tmp <- qqr(paternal.counts,malesize,colours[i],title=paste0(cohort," paternal, expected counts"),limits=range(9,45),xlowerlimit=10,mypch=16)
    } else {
        points(jitter(qqr(paternal.counts,malesize,makeplot=0)),pch=16,col=colours[i])
    }
}

for(i in 1:length(cutoffs)){
    p=cutoffs[i]
    paternal.counts=paternal[paternal$informative,paste("nrec.p.",p,sep="")]
    if(p==0){
        tmp <- qqr(paternal.counts,malesize,colours[i],title=paste0(cohort," paternal, hard cutoff"),limits=range(9,45),xlowerlimit=10,mypch=17)
    } else {
        points(jitter(qqr(paternal.counts,malesize,makeplot=0)),col=colours[i],pch=17)
    }
}
}
mtext("Expected recombinations genome-wide",1,outer=T,cex=1.3,padj=1.2)
mtext("Observed recombinations genome-wide",2,outer=T,cex=1.3,padj=-0.5)

dev.off()




