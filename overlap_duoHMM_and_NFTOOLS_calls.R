#if(FALSE){
  argv <- commandArgs(TRUE)
i=as.numeric(argv[1])
load("NFTOOLS_results/NFTOOLS_counts.all_families_except_outliers.more_stringent.RData")
duohmm.files=c("FC/duoHMM_results/FC.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"GPC/duoHMM_results/GPC.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"NTR/duoHMM_results/NTR.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"ORCADES/duoHMM_results/ORCADES.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"QTR/duoHMM_results/QTR_minus_MZs.370K.all.post_MERLIN.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"QTR/duoHMM_results/QTR_minus_MZs.610K.all.post_MERLIN.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"VB/duoHMM_results/VB.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt")
names(duohmm.files)=c("FC","GPC","NTR","ORCADES","QTR370","QTR610","VB")

#"FVG/duoHMM_results/FVG.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt"
#"VIS_KORCULA/duoHMM_results/KORCULA.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt"
#"VIS_KORCULA/duoHMM_results/VIS.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt"
#"CARL/duoHMM_results/CARL.all.more_stringent.generr_removed.recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt"

nftools.files=c("FC/NFTOOLS/output_data/FC.1.more_stringent.inrecomb_k5.list.txt.events","GPC/NFTOOLS/output_data/GPC.1.more_stringent.inrecomb_k5.list.txt.events",
    "NTR/NFTOOLS/output_data/NTR.1.more_stringent.inrecomb_k5.list.txt.events","ORCADES/NFTOOLS/output_data/ORCADES.1.more_stringent.inrecomb_k5.list.txt.events",
    "QTR/NFTOOLS/output_data/QTR370.1.more_stringent.inrecomb_k5.list.txt.events","QTR/NFTOOLS/output_data/QTR610.1.more_stringent.inrecomb_k5.list.txt.events",
    "VB/NFTOOLS/output_data/VB.1.more_stringent.inrecomb_k5.list.txt.events")
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
    n.overlaps.pat = sapply(unique(duohmm.pat$CHILD),function(child){
        child.duohmm.crossovers=GRanges(seqnames=duohmm.pat[duohmm.pat$CHILD==child,"chr"],ranges=IRanges(start=duohmm.pat[duohmm.pat$CHILD==child,"START"],
                                                                                             end=duohmm.pat[duohmm.pat$CHILD==child,"END"]))
        child.nftools.crossovers=(GRanges(seqnames=nftools.pat[nftools.pat$child ==child,"chr"],ranges=IRanges(start=nftools.pat[nftools.pat$child ==child,"left"],
                                                                                                  nftools.pat[nftools.pat$child ==child,"right"])))
        return(sum(countOverlaps(child.duohmm.crossovers,child.nftools.crossovers)))        
    })

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


}
