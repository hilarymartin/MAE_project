
setwd( "/gpfs1/well/donnelly/hilary/maternal_age_and_recombination/")
p.cutoff=c(0,0.5,0.9,0.95)

#files=c("QTR/duoHMM_results/QTR.370K.generr_removed.total_recombination_counts.paternal.p_0.txt","QTR/duoHMM_results/QTR.610K.generr_removed.total_recombination_counts.paternal.p_0.txt",
#"NTR/duoHMM_results/NTR.generr_removed.total_recombination_counts.paternal.p_0.txt","FC/duoHMM_results/FC.generr_removed.total_recombination_counts.paternal.p_0.txt",
#"ORCADES/duoHMM_results/ORCADES_raw_merged.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.txt","CARL/duoHMM_results/CARL.generr_removed.total_recombination_counts.paternal.p_0.txt",
#"VB/duoHMM_results/VB.generr_removed.total_recombination_counts.paternal.p_0.txt","FVG/duoHMM_results/FVG.generr_removed.total_recombination_counts.paternal.p_0.txt",
#"VIS_KORCULA/duoHMM_results/VIS.generr_removed.total_recombination_counts.paternal.p_0.txt","VIS_KORCULA/duoHMM_results/KORCULA.generr_removed.total_recombination_counts.paternal.p_0.txt")

files=c("QTR/duoHMM_results/QTR.370K.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.txt","QTR/duoHMM_results/QTR.610K.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.txt","NTR/duoHMM_results/NTR.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.txt","FC/duoHMM_results/FC.post_MERLIN.generr_removed.total_recombination_counts.paternal.p_0.txt","ORCADES/duoHMM_results/ORCADES_raw_merged.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.txt","CARL/duoHMM_results/CARL.generr_removed.total_recombination_counts.paternal.p_0.txt","VB/duoHMM_results/VB.generr_removed.total_recombination_counts.paternal.p_0.txt","FVG/duoHMM_results/FVG.generr_removed.total_recombination_counts.paternal.p_0.txt","VIS_KORCULA/duoHMM_results/VIS.post_MERLIN.generr_removed.total_recombination_counts.paternal.p_0.txt","VIS_KORCULA/duoHMM_results/KORCULA.post_MERLIN.generr_removed.total_recombination_counts.paternal.p_0.txt")


#files=c("QTR/duoHMM_results/QTR.370K.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","QTR/duoHMM_results/QTR.610K.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","NTR/duoHMM_results/NTR.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","FC/duoHMM_results/FC.post_MERLIN.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","ORCADES/duoHMM_results/ORCADES_raw_merged.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","CARL/duoHMM_results/CARL.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","VB/duoHMM_results/VB.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","FVG/duoHMM_results/FVG.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","VIS_KORCULA/duoHMM_results/VIS.post_MERLIN.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","VIS_KORCULA/duoHMM_results/KORCULA.post_MERLIN.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt")
names(files)=c("QTR.370","QTR.610K","NTR","FC","ORCADES.intersection","CARL","VB","FVG","VIS","KORCULA")

#files=c("QTR/duoHMM_results/QTR.370K.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","QTR/duoHMM_results/QTR.610K.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","NTR/duoHMM_results/NTR.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","FC/duoHMM_results/FC.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","ORCADES/duoHMM_results/ORCADES_raw_merged.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","CARL/duoHMM_results/CARL.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","VB/duoHMM_results/VB.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","FVG/duoHMM_results/FVG.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","VIS_KORCULA/duoHMM_results/VIS.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt","VIS_KORCULA/duoHMM_results/KORCULA.generr_removed.total_recombination_counts.paternal.p_0.no_double_recombinants_within_1000000bp.txt")
#names(files)=c("QTR.370","QTR.610K","NTR","FC","ORCADES.intersection","CARL","VB","FVG","VIS","KORCULA")

#mydir="/home/hilary/maternal_age_recombination/exploring_duoHMM_calls"
mydir="/home/hilary/maternal_age_recombination/exploring_duoHMM_calls_post_MERLIN"
#mydir="/home/hilary/maternal_age_recombination/exploring_duoHMM_calls_post_MERLIN_no_double_recombinants_within_1Mb/"
#mydir="/home/hilary/maternal_age_recombination/exploring_duoHMM_calls_no_double_recombinants_within_1Mb/"


remove.outliers=function(mat.counts,i){
#remove families in FC cohort due to MZ twins/duplicates (76), and because of low paternal counts, reason unknown (52)
if(i==4){
mat.counts=mat.counts[!mat.counts[,1] %in% c("52", "76"),]
}
#remove families FVG to due to pedigree erros
if(i==8){
mat.counts=mat.counts[!mat.counts[,1] %in% c("69","148"),]
}
#remove families in VB: 143 and 64 each contain a cryptic pair of MZ twins or duplicate samples
if(i==7){
mat.counts=mat.counts[!mat.counts[,1] %in% c("64","143"),]
}
#remove families from NTR: 10689 seems to have three times the same sample; 11880 is strange (maybe too much missingness in father)
if(i==3){
mat.counts=mat.counts[!mat.counts[,1] %in% c("10689","11880"),]
}
return(mat.counts)
}

for(p in p.cutoff){
matfiles=gsub("p_0",paste("p_",p,sep=""),files)
#pdf(paste(mydir,"/","total_paternal_recombination_counts.p_",p,".pdf",sep=""),height=5*3,width=2*5)
pdf(paste(mydir,"/","total_paternal_recombination_counts.excluding_problem_families.p_",p,".pdf",sep=""),height=5*3,width=2*5)
par(mfrow=c(5,2))
for(i in 1:length(matfiles)){
mat.counts=read.delim(matfiles[i],header=T)
mat.counts=remove.outliers(mat.counts,i)
hist1=hist(mat.counts$total,plot=F)
hist(mat.counts$total,main=paste("paternal counts ",names(files)[i],sep=""),xlab=paste("crossover count, p=",p,sep=""),xlim=range(hist1$breaks),ylim=c(0,max(hist1$counts)))
#if(i==1){
legend("topright",c("all duos","duos from families with >2 kids"),fill=c("white","grey"),cex=0.8)
legend("topleft",c("mean","median","all duos","duos from families with >2 kids"),lty=c(1,1,1,3),col=c("red","blue","black","black"),cex=0.8)
#}
if(sum(mat.counts$total_same_father>2)>0){
mat.counts2=mat.counts[mat.counts$total_same_father>2,]
hist(mat.counts2$total,main="",xlab="",col="grey",add=T)
abline(v=median(mat.counts2$total),col="red",lty=3)
abline(v=mean(mat.counts2$total),col="blue",lty=3)
}
abline(v=median(mat.counts$total),col="red")
abline(v=mean(mat.counts$total),col="blue")
}
dev.off()
}

#print out quantiles of distribution
for(p in p.cutoff){
matfiles=gsub("p_0",paste("p_",p,sep=""),files)
stats=NULL

for(i in 1:length(matfiles)){
cohort.stats=c()
mat.counts=read.delim(matfiles[i],header=T)
mat.counts=remove.outliers(mat.counts,i)
cohort.stats=c(cohort.stats,nrow(mat.counts),mean(mat.counts$total),quantile(mat.counts$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts$total))

#at least 3 kids with same father; only two generations 
if(sum(mat.counts$total_same_father>2 & mat.counts$n_genotyped_paternal_grandparents==0)>0){
mat.counts2=mat.counts[mat.counts$total_same_father>2 & mat.counts$n_genotyped_paternal_grandparents==0,]
cohort.stats=c(cohort.stats,nrow(mat.counts2),mean(mat.counts2$total),quantile(mat.counts2$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts2$total))
}else {
cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)
}

#at least 3 kids with same father and mother; only two generations 
if(sum(mat.counts$total_same_father>2 & mat.counts$mother !="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0)>0){
mat.counts2=mat.counts[mat.counts$total_same_father>2 & mat.counts$mother !="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0,]
cohort.stats=c(cohort.stats,nrow(mat.counts2),mean(mat.counts2$total),quantile(mat.counts2$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts2$total))
}else {
cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)
}


#at least 3 kids with same father but mother missing for all; only two generations 
if(sum(mat.counts$total_same_father>2 & mat.counts$mother =="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0)>0){
mat.counts2=mat.counts[mat.counts$total_same_father>2 & mat.counts$mother =="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0,]
cohort.stats=c(cohort.stats,nrow(mat.counts2),mean(mat.counts2$total),quantile(mat.counts2$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts2$total))
}else {
cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)
}

#2 kids with same father and mother; only two generations 
if(sum(mat.counts$total_same_father==2 & mat.counts$mother !="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0)>0){
mat.counts2=mat.counts[mat.counts$total_same_father==2 & mat.counts$mother !="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0,]
cohort.stats=c(cohort.stats,nrow(mat.counts2),mean(mat.counts2$total),quantile(mat.counts2$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts2$total))
}else {
cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)
}

#2 kids with same father but mother missing for all; only two generations 
if(sum(mat.counts$total_same_father==2 & mat.counts$mother =="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0)>0){
mat.counts2=mat.counts[mat.counts$total_same_father==2 & mat.counts$mother =="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0,]
cohort.stats=c(cohort.stats,nrow(mat.counts2),mean(mat.counts2$total),quantile(mat.counts2$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts2$total))
}else {
cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)
}

#1 kid with same father and mother; only two generations 
if(sum(mat.counts$total_same_father==1 & mat.counts$mother !="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0)>0){
mat.counts2=mat.counts[mat.counts$total_same_father==1 & mat.counts$mother !="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0,]
cohort.stats=c(cohort.stats,nrow(mat.counts2),mean(mat.counts2$total),quantile(mat.counts2$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts2$total))
}else {
cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)
}

#1 kid with same father but mother missing for all; only two generations 
if(sum(mat.counts$total_same_father==1 & mat.counts$mother =="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0)>0){
mat.counts2=mat.counts[mat.counts$total_same_father==1 & mat.counts$mother =="0" & mat.counts$total_same_parents==mat.counts$total_same_father & mat.counts$n_genotyped_paternal_grandparents==0,]
cohort.stats=c(cohort.stats,nrow(mat.counts2),mean(mat.counts2$total),quantile(mat.counts2$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts2$total))
}else {
cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)
}

#1 or both paternal grandparents
if(sum(mat.counts$n_genotyped_paternal_grandparents>0)>0){
mat.counts2=mat.counts[mat.counts$n_genotyped_paternal_grandparents>0,]
cohort.stats=c(cohort.stats,nrow(mat.counts2),mean(mat.counts2$total),quantile(mat.counts2$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts2$total))
}else {
cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)
}

#1 or both paternal grandparents, both parents same and genotyped
if(sum(mat.counts$n_genotyped_paternal_grandparents>0 & mat.counts$mother !="0" & mat.counts$total_same_parents==mat.counts$total_same_father )>0){
mat.counts2=mat.counts[mat.counts$n_genotyped_paternal_grandparents>0 & mat.counts$mother !="0" & mat.counts$total_same_parents==mat.counts$total_same_father,]
cohort.stats=c(cohort.stats,nrow(mat.counts2),mean(mat.counts2$total),quantile(mat.counts2$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts2$total))
}else {
cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)
}


#1 or both paternal grandparents, mother missing
if(sum(mat.counts$n_genotyped_paternal_grandparents>0 & mat.counts$mother =="0" & mat.counts$total_same_parents==mat.counts$total_same_father )>0){
mat.counts2=mat.counts[mat.counts$n_genotyped_paternal_grandparents>0 & mat.counts$mother =="0" & mat.counts$total_same_parents==mat.counts$total_same_father,]
cohort.stats=c(cohort.stats,nrow(mat.counts2),mean(mat.counts2$total),quantile(mat.counts2$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts2$total))
}else {
cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)
}


#either all 3 kids with same father or 1 or both paternal grandparents
if(sum(mat.counts$n_genotyped_paternal_grandparents>0|mat.counts$total_same_father>2)>0){
mat.counts2=mat.counts[mat.counts$n_genotyped_paternal_grandparents>0|mat.counts$total_same_father>2,]
cohort.stats=c(cohort.stats,nrow(mat.counts2),mean(mat.counts2$total),quantile(mat.counts2$total,c(0,0.25,0.5,0.75,1)),sd(mat.counts2$total))
}else {
cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)
}
stats=rbind(stats,cohort.stats)
}
rownames(stats)=names(files)
colnames(stats)=c(sapply(c("all_families","min_3_kids","min_3_kids.both_parent_same","min_3_kids.mother_missing","2_kids.both_parent_same","2_kids.mother_missing","1_kid.both_parents","1_kid.mother_missing","at_least_1_mat_grandparent",
"at_least_1_mat_grandparent.both_parents_same","at_least_1_mat_grandparent.mother_missing","informative"),function(x){paste(x,c("n_duos","mean","min","25th_pc","median","75th_pc","max","stdev"),sep=".")}))
#write.table(stats,paste(mydir,"/","summary_of_total_paternal_recombination_counts.p_",p,".txt",sep=""),sep="\t",quote=F)
write.table(stats,paste(mydir,"/","summary_of_total_paternal_recombination_counts.p_",p,".excluding_problem_families.txt",sep=""),sep="\t",quote=F)

}

for(p in p.cutoff){
matfiles=gsub("p_0",paste("p_",p,sep=""),files)
pdf(paste(mydir,"/","total_paternal_recombination_counts.excluding_problem_families.p_",p,".all_vs_three_generation_pedigrees.pdf",sep=""),height=5*3,width=2*5)
#pdf(paste(mydir,"/","total_paternal_recombination_counts.p_",p,".all_vs_three_generation_pedigrees.pdf",sep=""),height=5*3,width=2*5)
par(mfrow=c(5,2))
for(i in 1:length(matfiles)){
mat.counts=read.delim(matfiles[i],header=T)
mat.counts=remove.outliers(mat.counts,i)

hist1=hist(mat.counts$total,plot=F)
hist(mat.counts$total,main=paste("paternal counts ",names(files)[i],sep=""),xlab=paste("crossover count, p=",p,sep=""),xlim=range(hist1$breaks),ylim=c(0,max(hist1$counts)))
#if(i==1){
legend("topright",c("all duos","duos from three-generation families"),fill=c("white","grey"),cex=0.8)
legend("topleft",c("mean","median","all duos","duos from three-generation families"),lty=c(1,1,1,3),col=c("red","blue","black","black"),cex=0.8)
#}
if(sum(mat.counts$n_genotyped_paternal_grandparents>0)>0){
mat.counts2=mat.counts[mat.counts$n_genotyped_paternal_grandparents>0,]
hist(mat.counts2$total,main="",xlab="",col="grey",add=T)
abline(v=median(mat.counts2$total),col="red",lty=3)
abline(v=mean(mat.counts2$total),col="blue",lty=3)
}
abline(v=median(mat.counts$total),col="red")
abline(v=mean(mat.counts$total),col="blue")
}
dev.off()
}



for(p in p.cutoff){
matfiles=gsub("p_0",paste("p_",p,sep=""),files)
pdf(paste(mydir,"/","total_paternal_recombination_counts.excluding_problem_families.p_",p,".all_vs_informative_meioses.pdf",sep=""),height=5*3,width=2*5)
#pdf(paste(mydir,"/","total_paternal_recombination_counts.p_",p,".all_vs_informative_meioses.pdf",sep=""),height=5*3,width=2*5)
par(mfrow=c(5,2))
for(i in 1:length(matfiles)){
mat.counts=read.delim(matfiles[i],header=T)
mat.counts=remove.outliers(mat.counts,i)

hist1=hist(mat.counts$total,plot=F)
hist(mat.counts$total,main=paste("paternal counts ",names(files)[i],sep=""),xlab=paste("crossover count, p=",p,sep=""),xlim=range(hist1$breaks),ylim=c(0,max(hist1$counts)))
#if(i==1){
legend("topright",c("all duos","duos from informative meioses"),fill=c("white","grey"),cex=0.8)
legend("topleft",c("mean","median","all duos","duos from informative meioses"),lty=c(1,1,1,3),col=c("red","blue","black","black"),cex=0.8)
#}
if(sum(mat.counts$n_genotyped_paternal_grandparents>0|mat.counts$total_same_father>2)>0){
mat.counts2=mat.counts[mat.counts$n_genotyped_paternal_grandparents>0|mat.counts$total_same_father>2,]
hist(mat.counts2$total,main="",xlab="",col="grey",add=T)
abline(v=median(mat.counts2$total),col="red",lty=3)
abline(v=mean(mat.counts2$total),col="blue",lty=3)
}
abline(v=median(mat.counts$total),col="red")
abline(v=mean(mat.counts$total),col="blue")
}
dev.off()
}




