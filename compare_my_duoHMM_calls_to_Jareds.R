setwd("/well/donnelly/hilary/maternal_age_and_recombination/")
jared.files=c("orkney889_b37-shapeit.nrec","clean-carl-flip_b37-shapeit.nrec","valborbera_b37-shapeit.nrec","clean-fvg-flip_b37-shapeit.nrec","VIS_for1000G-shapeit.nrec","KORCULA_for1000G-shapeit.nrec")
names(jared.files)=c("ORCADES","CARL","VB","FVG","VIS","KORCULA")

files=c("ORCADES/duoHMM_results/ORCADES_raw_merged.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.5.txt","CARL/duoHMM_results/CARL.generr_removed.total_recombination_counts.paternal.p_0.5.txt",
    "VB/duoHMM_results/VB.generr_removed.total_recombination_counts.paternal.p_0.5.txt","FVG/duoHMM_results/FVG.generr_removed.total_recombination_counts.paternal.p_0.5.txt",
    "VIS_KORCULA/duoHMM_results/VIS.generr_removed.total_recombination_counts.paternal.p_0.5.txt","VIS_KORCULA/duoHMM_results/KORCULA.generr_removed.total_recombination_counts.paternal.p_0.5.txt")
names(files)=c("ORCADES","CARL","VB","FVG","VIS","KORCULA")

remove.outliers=function(mat.counts,i){
    #remove families FVG to due to pedigree erros
    if(i==4){
        mat.counts=mat.counts[!mat.counts[,1] %in% c("69","148"),]
    }
    #remove families in VB: 143 and 64 each contain a cryptic pair of MZ twins or duplicate samples
    if(i==3){
        mat.counts=mat.counts[!mat.counts[,1] %in% c("64","143"),]
    }
    return(mat.counts)
}


library(scales)


myxlabs=c("my counts (with duoHMM error correction)","my counts (with MERLIN+duoHMM error correction)")[c(1,2,2,2,1,1)]

for(i in 1:length(jared.files)){
jared.counts=read.delim(paste0("Jareds_duoHMM_calls/",jared.files[i]),sep="",header=F,stringsAsFactors=F)

colnames(jared.counts)=c("child_parent","chr","n_crossovers","parental_sex")
jared.counts$child = sapply(jared.counts[,1],function(x){z=strsplit(x,"-"); return(z[[1]][1])})
jared.counts$parent = sapply(jared.counts[,1],function(x){z=strsplit(x,"-"); return(z[[1]][2])})
jared.kids=unique(jared.counts$child)

my.pat.counts=read.delim(files[i],header=T,stringsAsFactors=F)
my.mat.counts=read.delim(gsub("paternal","maternal",files[i]),header=T,stringsAsFactors=F)
my.pat.counts=remove.outliers(my.pat.counts,i)
my.mat.counts=remove.outliers(my.mat.counts,i)

my.pat.counts$child=as.character(my.pat.counts$child)
my.mat.counts$child=as.character(my.mat.counts$child)
my.pat.counts$mother=as.character(my.pat.counts$mother)
my.mat.counts$mother=as.character(my.mat.counts$mother)
my.pat.counts$father=as.character(my.pat.counts$father)
my.mat.counts$father=as.character(my.mat.counts$father)

if(i==1){
my.pat.counts$child=gsub("[A-Z]","",my.pat.counts$child,perl=T)
my.mat.counts$child=gsub("[A-Z]","",my.mat.counts$child,perl=T)
my.pat.counts$mother=gsub("[A-Z]","",my.pat.counts$mother,perl=T)
my.mat.counts$mother=gsub("[A-Z]","",my.mat.counts$mother,perl=T)
my.pat.counts$father=gsub("[A-Z]","",my.pat.counts$father,perl=T)
my.mat.counts$father=gsub("[A-Z]","",my.mat.counts$father,perl=T)
}

my.kids=unique(c(my.mat.counts$child,my.pat.counts$child))

jared.mat.counts.by.chr=jared.counts[jared.counts$parental_sex==2,]
jared.pat.counts.by.chr=jared.counts[jared.counts$parental_sex==1,]

my.pat.counts=my.pat.counts[my.pat.counts$child %in% jared.pat.counts.by.chr$child,]
my.mat.counts=my.mat.counts[my.mat.counts$child %in% jared.mat.counts.by.chr$child,]

jared.mat.kids=unique(jared.mat.counts.by.chr$child)
jared.pat.kids=unique(jared.pat.counts.by.chr$child)
jared.mat.counts=data.frame(sapply(1:22,function(chr){sapply(jared.mat.kids,function(kid){    if(sum(jared.mat.counts.by.chr$chr==chr & jared.mat.counts.by.chr$child== kid)>0){
    return(jared.mat.counts.by.chr[jared.mat.counts.by.chr$chr==chr & jared.mat.counts.by.chr$child== kid,"n_crossovers"])}else {    return(0)}})}))

jared.pat.counts=data.frame(sapply(1:22,function(chr){sapply(jared.pat.kids,function(kid){    if(sum(jared.pat.counts.by.chr$chr==chr & jared.pat.counts.by.chr$child== kid)>0){
    return(jared.pat.counts.by.chr[jared.pat.counts.by.chr$chr==chr & jared.pat.counts.by.chr$child== kid,"n_crossovers"])}else {    return(0)}})}))

colnames(jared.mat.counts)=paste0("chr",1:22)
colnames(jared.pat.counts)=paste0("chr",1:22)
rownames(jared.mat.counts)=jared.mat.kids
rownames(jared.pat.counts)=jared.pat.kids

jared.mat.counts$total=rowSums(jared.mat.counts)
jared.pat.counts$total=rowSums(jared.pat.counts)
pdf(paste0("/home/hilary/maternal_age_recombination/comparing_my_duoHMM_calls_to_Jareds/comparing_Jareds_duoHMM_calls_with_MERLIN_error_correction_to_mine_without.",names(files)[i],".p_0.05.pdf"),height=5,width=10)
par(mfrow=c(1,2))
plot(my.mat.counts$total,jared.mat.counts[my.mat.counts$child,"total"],xlab=myxlabs[i],ylab="Jared's counts with MERLIN error correction",main=paste0("maternal ",names(files)[i]),
     xlim=range(c(my.mat.counts$total,jared.mat.counts$total)), ylim=range(c(my.mat.counts$total,jared.mat.counts$total)),col=alpha("black",0.4),bg=alpha("black",0.4),pch=21)
abline(a=0,b=1,col="red")
legend("topleft",paste0("rho = ",round(cor(my.mat.counts$total,jared.mat.counts[my.mat.counts$child,"total"],method="pearson"),digits=2)))
plot(my.pat.counts$total,jared.pat.counts[my.pat.counts$child,"total"],xlab=myxlabs[i],ylab="Jared's counts with MERLIN error correction",main=paste0("paternal ",names(files)[i]),
     xlim=range(c(my.pat.counts$total,jared.pat.counts$total)), ylim=range(c(my.pat.counts$total,jared.pat.counts$total)),col=alpha("black",0.4),bg=alpha("black",0.4),pch=21)
abline(a=0,b=1,col="red")
legend("topleft",paste0("rho = ",round(cor(my.pat.counts$total,jared.pat.counts[my.pat.counts$child,"total"],method="pearson"),digits=2)))
dev.off()

}

