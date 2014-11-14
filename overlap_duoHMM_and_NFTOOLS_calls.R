

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
for(i in 1:length(duohmm.files)){
    duohmm=read.delim(duohmm.files[i],header=T,stringsAsFactors=F)
    nftools=NULL
    for(chr in 1:22){
        nftools.chr=read.delim(gsub("1.",paste0(chr,"."),nftools.files[i],fixed=T),header=T,stringsAsFactors=F)
        nftools.chr$chr=chr
        nftools=rbind(nftools,nftools.chr)
    }
    nftools.mat=nftools[nftools$sex=="female",]
    nftools.pat=nftools[nftools$sex=="male",]
    duohmm.mat=duohmm[duohmm$sex=="Female" & duohmm$CHILD %in% nftools.mat$child,]
    duohmm.pat=duohmm[duohmm$sex=="Male" & duohmm$CHILD %in% nftools.pat$child,]
    nftools.mat=nftools.mat[ nftools.mat$child  %in% duohmm.mat$CHILD ,]
    nftools.pat=nftools.pat[ nftools.pat$child  %in% duohmm.pat$CHILD ,]
    
    n.overlaps = sapply(unique(duohmm.mat$CHILD),function(child){
        child.duohmm.crossovers=GRanges(seqnames=duohmm.mat[duohmm.mat$CHILD==child,"chr"],ranges=IRanges(start=duohmm.mat[duohmm.mat$CHILD==child,"START"],
                                                                                               end=duohmm.mat[duohmm.mat$CHILD==child,"END"]))
        child.nftools.crossovers=(GRanges(seqnames=nftools.mat[nftools.mat$child ==child,"chr"],ranges=IRanges(start=nftools.mat[nftools.mat$child ==child,"left"],
                                                                                                    nftools.mat[nftools.mat$child ==child,"right"])))
        return(sum(countOverlaps(child.duohmm.crossovers,child.nftools.crossovers)))        
    })
#why are there so many more of nftools crossovers not overlapping duoHMM crossovers than vice versa? some bug... duplications in phase uninformed families?
