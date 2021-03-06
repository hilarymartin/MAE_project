
mendel=list()
imendel=list()
#whole genome for FC
mendel[[1]]=read.delim("FC/FC_raw.pedigrees.lmendel",header=T,sep="",stringsAsFactors=F)
imendel[[1]]=read.delim("FC/FC_raw.pedigrees.imendel",header=T,sep="",stringsAsFactors=F)

#by chromosome
mendel.cohorts=c("ORCADES/input_for_SHAPEIT2/ORCADES_raw_merged.1.lmendel","VIS_KORCULA/MERLIN/VIS.chr1.lmendel","VIS_KORCULA/MERLIN/KORCULA.chr1.lmendel","GPC/MERLIN/GPC.chr1.round3.lmendel")
prefixes=c("ORCADES/input_for_SHAPEIT2/ORCADES_raw_merged.","VIS_KORCULA/MERLIN/VIS.chr","VIS_KORCULA/MERLIN/KORCULA.chr","GPC/MERLIN/GPC.chr")
suffixes=c(".lmendel",".lmendel",".lmendel",".round3.lmendel")
suffixes2=c(".imendel",".imendel",".imendel",".round3.imendel")
names(prefixes)=c("ORCADES","VIS","KORCULA","GPC")

for(c in 1:length(mendel.cohorts)){
    for(i in 1:22){
        cohort.mendel=read.delim(paste0(prefixes[c],i,suffixes[c]),header=T,sep="",stringsAsFactors=F)
        if(i==1){
            mendel[[c+1]]=cohort.mendel
            imendel[[c+1]]=read.delim(paste0(prefixes[c],i,suffixes2[c]),header=T,sep="",stringsAsFactors=F)
        } else {
            mendel[[c+1]]=rbind( mendel[[c+1]],cohort.mendel)
        }
    }
}
names(mendel) = c("FC",names(prefixes))
names(imendel)=names(mendel)

###missingness
missing=list()
imissing=list()

missing.files=c("QTR/input_for_SHAPEIT2/QTR_minus_MZs.610K.1.post_MERLIN.max_0.05_missing.lmiss","QTR/input_for_SHAPEIT2/QTR_minus_MZs.370K.1.post_MERLIN.max_0.05_missing.lmiss",
    "NTR/MERLIN/NTR_v2_minus_MZs.1.no_errors.for_plink.lmiss","FC/MERLIN/FC.1.no_errors.genotyped_indivs_only.for_plink.lmiss","ORCADES/input_for_SHAPEIT2/ORCADES_raw_merged.1.no_mendel_errors.max_0.05_missing.lmiss",
    "FVG/input_for_SHAPEIT2/valborbera_b37-chr1.lmiss","VB/input_for_SHAPEIT2/clean-fvg-flip_b37-chr1.lmiss","CARL/input_for_SHAPEIT2/clean-carl-flip_b37-chr1.lmiss",
    "VIS_KORCULA/input_for_SHAPEIT2/VIS.chr1.no_mendel_errors.post_MERLIN.lmiss","VIS_KORCULA/input_for_SHAPEIT2/KORCULA.chr1.no_mendel_errors.post_MERLIN.lmiss",
    "GPC/input_for_SHAPEIT2/GPC.chr1.no_mendel_errors.post_MERLIN.lmiss")

prefixes=c("QTR/input_for_SHAPEIT2/QTR_minus_MZs.610K.","QTR/input_for_SHAPEIT2/QTR_minus_MZs.370K.","NTR/MERLIN/NTR_v2_minus_MZs.","FC/MERLIN/FC.","ORCADES/input_for_SHAPEIT2/ORCADES_raw_merged.",
    "FVG/input_for_SHAPEIT2/valborbera_b37-chr","VB/input_for_SHAPEIT2/clean-fvg-flip_b37-chr","CARL/input_for_SHAPEIT2/clean-carl-flip_b37-chr","VIS_KORCULA/input_for_SHAPEIT2/VIS.chr",
    "VIS_KORCULA/input_for_SHAPEIT2/KORCULA.chr","GPC/input_for_SHAPEIT2/GPC.chr")
suffixes=c(".post_MERLIN.max_0.05_missing.lmiss",".post_MERLIN.max_0.05_missing.lmiss",".no_errors.for_plink.lmiss",".no_errors.genotyped_indivs_only.for_plink.lmiss",".no_mendel_errors.max_0.05_missing.lmiss",".lmiss",
    ".lmiss",".lmiss",".no_mendel_errors.post_MERLIN.lmiss",".no_mendel_errors.post_MERLIN.lmiss",".no_mendel_errors.post_MERLIN.lmiss")

suffixes2=c(".post_MERLIN.max_0.05_missing.imiss",".post_MERLIN.max_0.05_missing.imiss",".no_errors.for_plink.imiss",".no_errors.genotyped_indivs_only.for_plink.imiss",".no_mendel_errors.max_0.05_missing.imiss",".imiss",
    ".imiss",".imiss",".no_mendel_errors.post_MERLIN.imiss",".no_mendel_errors.post_MERLIN.imiss",".no_mendel_errors.post_MERLIN.imiss")

names(missing.files)=c("QTR610","QTR370","NTR","FC","ORCADES","FVG","VB","CARL","VIS","KORCULA","GPC")

for(c in 1:length(missing.files)){
    for(i in 1:22){
        cohort.missing=read.delim(paste0(prefixes[c],i,suffixes[c]),header=T,sep="",stringsAsFactors=F)
        if(i==1){
            missing[[c]]=cohort.missing
            imissing[[c]]=read.delim(paste0(prefixes[c],i,suffixes2[c]),header=T,sep="",stringsAsFactors=F)
        } else {
            missing[[c]]=rbind( missing[[c]],cohort.missing)
        }
    }
}
names(missing) = names(missing.files)
names(imissing)=names(missing.files)

nindivs = sapply(imissing,nrow)

mendel.2=list()
for(m in 1:length(mendel)){
    d=mendel[[m]]
    e=missing[[names(mendel)[m]]]
    d=d[d$SNP %in% e$SNP,]
    d$prop.mendel=d$N/nrow(imendel[[m]])
    mendel.2[[m]]=d

}
names(mendel.2)=names(mendel)
pdf("/home/hilary/maternal_age_recombination/QC/proportion_of_Mendel_errors_by_SNP.final_SNP_set_used.pdf",height=6,width=9)
par(mfrow=c(2,3))
for(m in 1:length(mendel.2)){
    d=mendel.2[[m]]
    hist(d$prop.mendel,main=names(mendel.2)[m],xlab="Proportion with Mendel errors",ylab="Number of SNPs")
    cat(names(mendel.2)[m],"\n")
    print(nrow(d[d$prop.mendel>0.01,]))
}
dev.off()

mendel.summary= rbind(sapply(mendel.2,function(d){nrow(d[d$prop.mendel>0.0,])}),      sapply(mendel.2,function(d){nrow(d[d$prop.mendel>0.01,])}),sapply(mendel.2,function(d){nrow(d[d$prop.mendel>0.02,])}),
      sapply(mendel.2,function(d){nrow(d[d$prop.mendel>0.03,])}),sapply(mendel.2,function(d){nrow(d[d$prop.mendel>0.04,])}),sapply(mendel.2,function(d){nrow(d[d$prop.mendel>0.05,])}))
mendel.summary=rbind(mendel.summary,sapply(mendel.2,nrow))
rownames(mendel.summary)=c("ME>0","ME>0.01","ME>0.02","ME>0.03","ME>0.04","ME>0.05","total")
mendel.summary=rbind(mendel.summary,mendel.summary[1:6,]/mendel.summary[7,])
#rownames(mendel.summary)[8:nrow(mendel.summary)]=paste0("prop.SNPs.",rownames(mendel.summary)[8:nrow(mendel.summary)])
mendel.summary=rbind(mendel.summary,sapply(imendel,nrow))
rownames(mendel.summary)[nrow(mendel.summary)]="N.indivs"
write.table(mendel.summary,"/home/hilary/maternal_age_recombination/QC/proportion_of_SNPs_with_Mendel_errors.txt",quote=F,sep="\t")

pdf("/home/hilary/maternal_age_recombination/QC/proportion_of_missingness_SNP.final_SNP_set_used.pdf",height=5,width=5)
for(m in 1:length(missing)){
    d=missing[[m]]
    hist(d$F_MISS,main=names(missing)[m],xlab="Proportion missing",ylab="Number of SNPs")
}
dev.off()


pdf("/home/hilary/maternal_age_recombination/QC/proportion_of_missingness_SNP.final_SNP_set_used.pdf",height=9,width=12)
par(mfrow=c(3,4))
for(m in 1:length(missing)){
    d=missing[[m]]
    hist(d$F_MISS,main=names(missing)[m],xlab="Proportion missing",ylab="Number of SNPs",xlim=c(0,0.05))
}
dev.off()

missing.summary=rbind(sapply(missing,function(d){nrow(d[d$F_MISS>0,])}),sapply(missing,function(d){nrow(d[d$F_MISS>0.01,])}),sapply(missing,function(d){nrow(d[d$F_MISS>0.02,])}),
    sapply(missing,function(d){nrow(d[d$F_MISS>0.03,])}),sapply(missing,function(d){nrow(d[d$F_MISS>0.04,])}),sapply(missing,function(d){nrow(d[d$F_MISS>0.05,])}),sapply(missing,function(d){nrow(d)}))
rownames(missing.summary)=c("missing>0","missing>0.01","missing>0.02","missing>0.03","missing>0.04","missing>0.05","total")
write.table(missing.summary,"/home/hilary/maternal_age_recombination/QC/proportion_of_SNPs_with_missingness.txt",quote=F,sep="\t")
