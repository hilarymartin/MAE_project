#if(FALSE){
  argv <- commandArgs(TRUE)
i=as.numeric(argv[1])

load("NFTOOLS_results/NFTOOLS_counts.all_families_except_outliers.more_stringent.with_new_QTR.RData")

duohmm.files=c("FC/duoHMM_results/FC.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"GPC/duoHMM_results/GPC.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"NTR/duoHMM_results/NTR.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"ORCADES/duoHMM_results/ORCADES.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"QTR/duoHMM_results/all_QTR_370K_samples.all.post_MERLIN.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"QTR/duoHMM_results/all_QTR_610K_samples.all.post_MERLIN.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"QTR/duoHMM_results/new_samples.CoreExome.all.post_MERLIN.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
    "VB/duoHMM_results/VB.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt")
names(duohmm.files)=c("FC","GPC","NTR","ORCADES","QTR370","QTR610","QTRCoreExome","VB")

#"FVG/duoHMM_results/FVG.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt"
#"VIS_KORCULA/duoHMM_results/KORCULA.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt"
#"VIS_KORCULA/duoHMM_results/VIS.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt"
#"CARL/duoHMM_results/CARL.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt"

nftools.files=c("FC/NFTOOLS/output_data/FC.1.more_stringent.inrecomb_k5.list.txt.events","GPC/NFTOOLS/output_data/GPC.1.more_stringent.inrecomb_k5.list.txt.events",
    "NTR/NFTOOLS/output_data/NTR.1.more_stringent.inrecomb_k5.list.txt.events","ORCADES/NFTOOLS/output_data/ORCADES.1.more_stringent.inrecomb_k5.list.txt.events",
"QTR/NFTOOLS/output_data/all_QTR_370K_samples.1.post_MERLIN.more_stringent.inrecomb_k5.list.txt.events","QTR/NFTOOLS/output_data/all_QTR_610K_samples.1.post_MERLIN.more_stringent.inrecomb_k5.list.txt.events",
    "QTR/NFTOOLS/output_data/new_samples.CoreExome.1.post_MERLIN.more_stringent.inrecomb_k5.list.txt.events",    "VB/NFTOOLS/output_data/VB.1.more_stringent.inrecomb_k5.list.txt.events")
names(nftools.files)=names(duohmm.files)
library(GenomicRanges)


filtered.nftools=all.events[[names(duohmm.files)[i]]]
filtered.nftools.inf=filtered.nftools[filtered.nftools$nkids>1,]

duohmm=read.delim(duohmm.files[i],header=T,stringsAsFactors=F)
duohmm=duohmm[duohmm$informative.2gen.2parents,]

duohmm=duohmm[duohmm$CHILD %in% filtered.nftools.inf$child,]

nftools=NULL
for(chr in 1:22){
    nftools.chr=read.delim(gsub("1.",paste0(chr,"."),nftools.files[i],fixed=T),header=T,stringsAsFactors=F)
    nftools.chr$chr=chr
    nftools=rbind(nftools,nftools.chr)
}
nftools.mat=nftools[nftools$sex=="female",]
nftools.mat=nftools.mat[nftools.mat$child %in% filtered.nftools[filtered.nftools$sex=="female","child"],]
nftools.pat=nftools[nftools$sex=="male",]
nftools.pat=nftools.pat[nftools.pat$child %in% filtered.nftools[filtered.nftools$sex=="male","child"],]
duohmm.mat=duohmm[duohmm$sex=="Female" & duohmm$CHILD %in% nftools.mat$child,]
duohmm.pat=duohmm[duohmm$sex=="Male" & duohmm$CHILD %in% nftools.pat$child,]
nftools.mat=nftools.mat[ nftools.mat$child  %in% duohmm.mat$CHILD ,]
nftools.pat=nftools.pat[ nftools.pat$child  %in% duohmm.pat$CHILD ,]

n.overlaps.mat = sapply(unique(duohmm.mat$CHILD),function(child){
    child.duohmm.crossovers=GRanges(seqnames=duohmm.mat[duohmm.mat$CHILD==child,"chr"],ranges=IRanges(start=duohmm.mat[duohmm.mat$CHILD==child,"START"],
                                                                                           end=duohmm.mat[duohmm.mat$CHILD==child,"END"]))
    child.nftools.crossovers=(GRanges(seqnames=nftools.mat[nftools.mat$child ==child,"chr"],ranges=IRanges(start=nftools.mat[nftools.mat$child ==child,"left"],
                                                                                                  nftools.mat[nftools.mat$child ==child,"right"])))
    return(sum(countOverlaps(child.duohmm.crossovers,child.nftools.crossovers)))        
})


overlaps.mat = sapply(unique(duohmm.mat$CHILD),function(child){
    child.duohmm.crossovers=GRanges(seqnames=duohmm.mat[duohmm.mat$CHILD==child,"chr"],ranges=IRanges(start=duohmm.mat[duohmm.mat$CHILD==child,"START"],
                                                                                           end=duohmm.mat[duohmm.mat$CHILD==child,"END"]))
    child.nftools.crossovers=(GRanges(seqnames=nftools.mat[nftools.mat$child ==child,"chr"],ranges=IRanges(start=nftools.mat[nftools.mat$child ==child,"left"],
                                                                                                  nftools.mat[nftools.mat$child ==child,"right"])))
    return(subsetByOverlaps(child.duohmm.crossovers,child.nftools.crossovers))
})
if(FALSE){
overlapping = cbind(duohmm[queryHits(myoverlaps),],nftools[subjectHits(myoverlaps),])

#matching start and end
overlapping$type[overlapping$START==overlapping$left & overlapping$END==overlapping$right] = 1
#one contained within the other
overlapping$type[overlapping$START<overlapping$left & overlapping$END>overlapping$right] = 2
overlapping$type[overlapping$START>overlapping$left & overlapping$END<overlapping$right] = 3

#matching start, contained within
overlapping$type[overlapping$START==overlapping$left & overlapping$END<overlapping$right] = 4
#matching start, overlapping end
overlapping$type[overlapping$START==overlapping$left & overlapping$END>overlapping$right] = 4

overlapping$type[overlapping$START<overlapping$left & overlapping$END==overlapping$right] = 4
overlapping$type[overlapping$START>overlapping$left & overlapping$END==overlapping$right] = 4
}
n.overlaps.pat = sapply(unique(duohmm.pat$CHILD),function(child){
    child.duohmm.crossovers=GRanges(seqnames=duohmm.pat[duohmm.pat$CHILD==child,"chr"],ranges=IRanges(start=duohmm.pat[duohmm.pat$CHILD==child,"START"],
                                                                                           end=duohmm.pat[duohmm.pat$CHILD==child,"END"]))
    child.nftools.crossovers=(GRanges(seqnames=nftools.pat[nftools.pat$child ==child,"chr"],ranges=IRanges(start=nftools.pat[nftools.pat$child ==child,"left"],
                                                                                                nftools.pat[nftools.pat$child ==child,"right"])))
    return(sum(countOverlaps(child.duohmm.crossovers,child.nftools.crossovers)))        
})
names(n.overlaps.mat)=unique(duohmm.mat$CHILD)
names(n.overlaps.pat)=unique(duohmm.pat$CHILD)


n.duohmm.crossovers.mat = sapply(unique(duohmm.mat$CHILD),function(child){return(sum(duohmm.mat$CHILD==child))})
n.nftools.crossovers.mat = sapply(unique(duohmm.mat$CHILD),function(child){return(sum(nftools.mat$child==child))})

names(n.duohmm.crossovers.mat)=unique(duohmm.mat$CHILD)
names(n.nftools.crossovers.mat)=unique(duohmm.mat$CHILD)

n.duohmm.crossovers.pat = sapply(unique(duohmm.pat$CHILD),function(child){return(sum(duohmm.pat$CHILD==child))})
n.nftools.crossovers.pat = sapply(unique(duohmm.pat$CHILD),function(child){return(sum(nftools.pat$child==child))})

names(n.duohmm.crossovers.pat)=unique(duohmm.pat$CHILD)
names(n.nftools.crossovers.pat)=unique(duohmm.pat$CHILD)


mat.overlap.by.child = as.data.frame(cbind(n.overlaps.mat[match(names(n.overlaps.mat),names(n.duohmm.crossovers.mat))],n.duohmm.crossovers.mat,n.nftools.crossovers.mat))
pat.overlap.by.child = as.data.frame(cbind(n.overlaps.pat[match(names(n.overlaps.pat),names(n.duohmm.crossovers.pat))],n.duohmm.crossovers.pat,n.nftools.crossovers.pat))
colnames(mat.overlap.by.child) = c("overlap","duohmm","nftools")
colnames(pat.overlap.by.child) = c("overlap","duohmm","nftools")
mat.overlap.by.child$diff=mat.overlap.by.child$duohmm - mat.overlap.by.child$nftools
pat.overlap.by.child$diff=pat.overlap.by.child$duohmm - pat.overlap.by.child$nftools
if(FALSE){
write.table(mat.overlap.by.child,paste0("comparison_of_NFTOOLS_vs_duoHMM_calls/",names(duohmm.files)[i],"_maternal_crossovers_overlap_by_child_from_informative_nuclear_families.txt"),quote=F,sep="\t")
write.table(pat.overlap.by.child,paste0("comparison_of_NFTOOLS_vs_duoHMM_calls/",names(duohmm.files)[i],"_paternal_overlap_by_child_from_informative_nuclear_families.txt"),quote=F,sep="\t")
}
  overlaps=c(names(duohmm.files[i]),sum(n.overlaps.mat),nrow(duohmm.mat),nrow(nftools.mat),sum(n.overlaps.mat)/nrow(duohmm.mat),sum(n.overlaps.mat)/nrow(nftools.mat),  sum(n.overlaps.pat),nrow(duohmm.pat),nrow(nftools.pat),sum(n.overlaps.pat)/nrow(duohmm.pat),
  sum(n.overlaps.pat)/nrow(nftools.pat))

names(overlaps) = c("cohort","mat.n.overlaps","mat.n.duohmm","mat.n.nftools","mat.prop.overlap.duohmm","mat.prop.overlap.nftools","pat.n.overlaps","pat.n.duohmm","pat.n.nftools","pat.prop.overlap.duohmm","pat.prop.overlap.nftools")

write.table(overlaps,paste0("comparison_of_NFTOOLS_vs_duoHMM_calls/",names(duohmm.files)[i],"_overlapping_NFTOOLS_and_duoHMM_crossovers.informative_nuclear_families.txt"),quote=F,sep="\t")

#}
if(FALSE){
    
all.overlaps=NULL
for(i in 1:length(nftools.files)){
  overlaps=read.delim(paste0("comparison_of_NFTOOLS_vs_duoHMM_calls/",names(duohmm.files)[i],"_overlapping_NFTOOLS_and_duoHMM_crossovers.informative_nuclear_families.txt"),skip=2,header=F)
if(i==1){
  all.overlaps=overlaps
} else {
  all.overlaps=cbind(all.overlaps,overlaps)
}
}
all.overlaps=all.overlaps[,c(2,4,6,8,10,12,14)]
colnames(all.overlaps)=names(duohmm.files)
rownames(all.overlaps)=c("mat.n.overlaps","mat.n.duohmm","mat.n.nftools","mat.prop.overlap.duohmm","mat.prop.overlap.nftools","pat.n.overlaps","pat.n.duohmm","pat.n.nftools","pat.prop.overlap.duohmm","pat.prop.overlap.nftools")

all.overlaps$total=rowSums(all.overlaps)

all.overlaps["mat.prop.overlap.duohmm","total"]= all.overlaps["mat.n.overlaps","total"]/all.overlaps["mat.n.duohmm","total"]
all.overlaps["mat.prop.overlap.nftools","total"]= all.overlaps["mat.n.overlaps","total"]/all.overlaps["mat.n.nftools","total"]
all.overlaps["pat.prop.overlap.duohmm","total"]= all.overlaps["pat.n.overlaps","total"]/all.overlaps["pat.n.duohmm","total"]
all.overlaps["pat.prop.overlap.nftools","total"]= all.overlaps["pat.n.overlaps","total"]/all.overlaps["pat.n.nftools","total"]

write.table(all.overlaps,"comparison_of_NFTOOLS_vs_duoHMM_calls/all_cohorts.overlapping_NFTOOLS_and_duoHMM_crossovers.informative_nuclear_families.txt",quote=F,sep="\t")


all.nftools.mat=list()
all.nftools.pat=list()
all.duohmm.mat=list()
all.duohmm.pat=list()
for(i in 1:length(duohmm.files)){
filtered.nftools=all.events[[names(duohmm.files)[i]]]
filtered.nftools.inf=filtered.nftools[filtered.nftools$nkids>1,]

duohmm=read.delim(duohmm.files[i],header=T,stringsAsFactors=F)
duohmm=duohmm[duohmm$informative.2gen.2parents,]

duohmm=duohmm[duohmm$CHILD %in% filtered.nftools.inf$child,]

  nftools=NULL
  for(chr in 1:22){
    nftools.chr=read.delim(gsub("1.",paste0(chr,"."),nftools.files[i],fixed=T),header=T,stringsAsFactors=F)
        nftools.chr$chr=chr
    nftools=rbind(nftools,nftools.chr)
}
nftools.mat=nftools[nftools$sex=="female",]
nftools.mat=nftools.mat[nftools.mat$child %in% filtered.nftools[filtered.nftools$sex=="female","child"],]
nftools.pat=nftools[nftools$sex=="male",]
nftools.pat=nftools.pat[nftools.pat$child %in% filtered.nftools[filtered.nftools$sex=="male","child"],]
duohmm.mat=duohmm[duohmm$sex=="Female" & duohmm$CHILD %in% nftools.mat$child,]
duohmm.pat=duohmm[duohmm$sex=="Male" & duohmm$CHILD %in% nftools.pat$child,]
nftools.mat=nftools.mat[ nftools.mat$child  %in% duohmm.mat$CHILD ,]
nftools.pat=nftools.pat[ nftools.pat$child  %in% duohmm.pat$CHILD ,]
all.nftools.mat[[i]]=nftools.mat
all.nftools.pat[[i]]=nftools.pat

all.duohmm.mat[[i]]=duohmm.mat
all.duohmm.pat[[i]]=duohmm.pat

}

names(all.nftools.mat)=names(duohmm.files)
names(all.nftools.pat)=names(duohmm.files)
names(all.duohmm.mat)=names(duohmm.files)
names(all.duohmm.pat)=names(duohmm.files)


#read in fams and check if the individuals have other relatives in the cohort that could be affecting duoHMM calls
fams=c("FC/duoHMM_results/FC.more_stringent.sample_information.txt","GPC/duoHMM_results/GPC.more_stringent.sample_information.txt","NTR/duoHMM_results/NTR.more_stringent.sample_information.txt","ORCADES/duoHMM_results/ORCADES.more_stringent.sample_information.txt",
    "QTR/duoHMM_results/QTR370.more_stringent.sample_information.txt","QTR/duoHMM_results/QTR610.more_stringent.sample_information.txt","VB/duoHMM_results/VB.more_stringent.sample_information.txt")
names(fams) = names(duohmm.files)

all.fams=list()
for(i in 1:length(fams)){
    fam=read.delim(fams[i],header=T,sep="")
    fam[,2]=gsub("-","_",fam[,2])
    fam[,3]=gsub("-","_",fam[,3])
    fam[,4]=gsub("-","_",fam[,4])
    all.fams[[names(fams)[i]]] = fam
}

all.outliers=NULL
for(i in 1:length(nftools.files)){
    mat.overlap.by.child=read.delim(paste0("comparison_of_NFTOOLS_vs_duoHMM_calls/",names(duohmm.files)[i],"_maternal_crossovers_overlap_by_child_from_informative_nuclear_families.txt"),header=T)
    pat.overlap.by.child=read.delim(paste0("comparison_of_NFTOOLS_vs_duoHMM_calls/",names(duohmm.files)[i],"_paternal_overlap_by_child_from_informative_nuclear_families.txt"),header=T)
    cutoff=5
    mat.outliers=mat.overlap.by.child[abs(mat.overlap.by.child$diff) > 5,]
    if(nrow(mat.outliers)>0){
        mat.outliers$sex="female"
        mat.outliers$cohort=names(nftools.files)[i]
        mat.outliers$child=rownames(mat.outliers)
    }
    pat.outliers=pat.overlap.by.child[abs(pat.overlap.by.child$diff) > 5,]
    if(nrow(pat.outliers)>0){
        pat.outliers$sex="male"
        pat.outliers$cohort=names(nftools.files)[i]
        pat.outliers$child=rownames(pat.outliers)
    }
    
    if(nrow(mat.outliers)>0 &&nrow(pat.outliers)>0){
        all.outliers=rbind(all.outliers,mat.outliers,pat.outliers)
    } else if(nrow(mat.outliers)>0){
        all.outliers=rbind(all.outliers,mat.outliers)
    } else if(nrow(pat.outliers)>0){
        all.outliers=rbind(all.outliers,pat.outliers)
    }
}

all.outliers=cbind(all.outliers,matrix(0,ncol=44,nrow=nrow(all.outliers)))

colnames(all.outliers)[8:(8+21)] = paste0("duohmm.chr",1:22)
colnames(all.outliers)[(8+22):ncol(all.outliers)] = paste0("nftools.chr",1:22)

all.outliers$nkids=NA
all.outliers$total.family.size=NA

for(j in 1:nrow(all.outliers)){
    cohort.events=all.events[[all.outliers[j,"cohort"]]]
if(all.outliers[j,"sex"]=="female"){
    duohmm.counts=table(all.duohmm.mat[[all.outliers[j,"cohort"]]][ all.duohmm.mat[[all.outliers[j,"cohort"]]]$CHILD==all.outliers[j,"child"],"chr"])
    nftools.counts=        table(all.nftools.mat[[all.outliers[j,"cohort"]]][ all.nftools.mat[[all.outliers[j,"cohort"]]]$child==all.outliers[j,"child"],"chr"])
    all.outliers[j,"nkids"] = cohort.events[cohort.events$child == all.outliers[j,"child"] & cohort.events$sex=="female","nkids"]
} else {
    duohmm.counts= table(all.duohmm.pat[[all.outliers[j,"cohort"]]][ all.duohmm.pat[[all.outliers[j,"cohort"]]]$CHILD==all.outliers[j,"child"],"chr"])
    nftools.counts= table(all.nftools.pat[[all.outliers[j,"cohort"]]][ all.nftools.pat[[all.outliers[j,"cohort"]]]$child==all.outliers[j,"child"],"chr"])
    all.outliers[j,"nkids"] = cohort.events[cohort.events$child == all.outliers[j,"child"] & cohort.events$sex=="male","nkids"]
}
names(duohmm.counts)=paste0("duohmm.chr",names(duohmm.counts))
names(nftools.counts)=paste0("nftools.chr",names(nftools.counts))
all.outliers[j,names(duohmm.counts)] = duohmm.counts
all.outliers[j,names(nftools.counts)] = nftools.counts
   all.outliers[j,"total.family.size"] = all.fams[[all.outliers[j,"cohort"]]][all.fams[[all.outliers[j,"cohort"]]]$V2 == all.outliers[j,"child"],"total_fam_size"]


}
write.table(all.outliers,"comparison_of_NFTOOLS_vs_duoHMM_calls/individuals_with_min_5_crossovers_different_between_NFTOOLS_and_duoHMMM.txt",quote=F,sep="\t")
lapply(all.duohmm.mat,function(x){colSums(x[,grep("informative",colnames(x))])})

}
