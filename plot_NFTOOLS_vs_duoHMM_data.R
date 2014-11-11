 load("duoHMM_data_for_RSTAN.NTR_v2.RData")
duohmm.mat=data0.mat

 load("NFTOOLS_data_for_RSTAN.NTR_v2.RData")
nftools.mat = data1.mat


both$diff=both$Freq-both$nrec

both=cbind(nftools.mat,duohmm.mat[match(as.character(nftools.mat$child),as.character(duohmm.mat$CHILD)),])
both=both[!is.na(both$CHILD),]
library(scales)
both$col=NA
mycols=c("black","blue","red","green","orange","red4","purple","darkgreen")
names(mycols)=c("CARL","FC","FVG","GPC","NTR_v2","QTR_370K","QTR_610K","VB")
both$col=mycols[both$cohort]

pdf("/home/hilary/maternal_age_recombination/NFTOOLS_vs_duoHMM_counts.inform_nuc_fams_only.including_NTR_v2.pdf",height=5,width=5)
plot(x=jitter(both$Freq),y=jitter(both$nrec),col=alpha(both$col,0.5),xlab="NFTOOLS",ylab="duoHMM",pch=19)
abline(a=0,b=1)
legend("topleft",paste0("r = ",round(cor(both$Freq,both$nrec),digits=2))
dev.off()

pdf("/home/hilary/maternal_age_recombination/NFTOOLS_vs_duoHMM_counts.by_cohort.inform_nuc_fams_only.including_NTR_v2.pdf",height=5,width=8)
par(mfrow=c(2,3))
cohorts=unique(both$cohort)
for(i in 1:length(cohorts)){
    plot(x=jitter(both$Freq[both$cohort==cohorts[i]]),y=jitter(both$nrec[both$cohort==cohorts[i]]),col=alpha(both$col[both$cohort==cohorts[i]],0.5),xlab="NFTOOLS",ylab="duoHMM",pch=19,main=cohorts[i])
    abline(a=0,b=1)
legend("topleft",paste0("r = ",round(cor(both$Freq[both$cohort==cohorts[i]],both$nrec[both$cohort==cohorts[i]]),digits=2)))
}
dev.off()



       ###outliers both[abs(both$diff)>5,
 #NTR 10869   999009188 was an outlier --> 44 in NFTOOLS, 55 in duoHMM -- noted 2% missingness in data
#       999604629 and 999677284 - child of 999604628 - both outliers; looks like the crossovers may have been incorrectly assigned in duoHMM?
