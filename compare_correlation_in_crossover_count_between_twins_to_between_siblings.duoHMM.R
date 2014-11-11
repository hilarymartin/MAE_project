
compare.correlations.by.simulation=function(counts,method,n){
  twin.corr.test=cor.test(counts[,1],counts[,2],method=method)
  twin.corr=twin.corr.test$estimate
  twin.p=twin.corr.test$p.value
  sim.corr=c()
  for(i in 1:n){
    random.counts=apply(counts,1,function(x){return(x[sample(2,size=1)])})
    random.corr.test=cor.test(counts[,3],random.counts,method=method)
    random.corr=random.corr.test$estimate
    sim.corr=c(sim.corr,random.corr)
  }
  empirical.p=sum(sim.corr>=twin.corr)/n
  return(c(twin.corr,twin.p,empirical.p))
}

qtr.610.nftools=read.delim("/home/hilary/TwinRecombination/Brisbane/NFTOOLS_maternal_counts_for_QTR_610K_families.twins_and_sibling.txt",header=T,stringsAsFactors=F)
qtr.all.nftools=read.delim("/home/hilary/TwinRecombination/Brisbane/NFTOOLS_maternal_counts_for_all_QTR_families_610_or_370_separately.twins_and_sibling.txt",header=T,stringsAsFactors=F)
ntr.nftools=read.delim("/home/hilary/TwinRecombination/Amsterdam/NFTOOLS_maternal_counts_for_all_NTR_families.twins_and_sibling.txt",header=T,stringsAsFactors=F)

qtr.370.nftools=qtr.all.nftools[!qtr.all.nftools$mother %in% qtr.610.nftools$mother,]
qtr.370.nftools.mat.counts.by.child=qtr.370.nftools[,c("twin1_mat_count","twin2_mat_count","sib_mat_count")]

mat.comparison.pearson.qtr370.nftools.pearson = compare.correlations.by.simulation(qtr.370.nftools.mat.counts.by.child,"pearson",1000)
mat.comparison.spearman.qtr370.nftools.spearman = compare.correlations.by.simulation(qtr.370.nftools.mat.counts.by.child,"spearman",1000)

#qtr.370.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.370K.generr_removed.total_recombination_counts.maternal.p_0.9.txt",header=T,stringsAsFactors=F)
#qtr.370.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.370K.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.maternal.p_0.5.txt",header=T,stringsAsFactors=F)
#qtr.370.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.370K.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.maternal.p_0.5.no_double_recombinants_within_1000000bp.txt",header=T,stringsAsFactors=F)

#qtr.370.duohmm.all=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR_minus_MZs.370K.all.post_MERLIN.max_0.05_missing.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_and_clustered_duplicate_xovers_removed.txt",header=T,stringsAsFactors=F)
qtr.370.duohmm.all=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR_minus_MZs.370K.all.post_MERLIN.max_0.05_missing.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",header=T,stringsAsFactors=F)
#qtr.610.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.610K.generr_removed.total_recombination_counts.maternal.p_0.9.txt",header=T,stringsAsFactors=F)
#qtr.610.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.610K.generr_removed.total_recombination_counts.maternal.p_0.5.txt",header=T,stringsAsFactors=F)
#qtr.610.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.610K.generr_removed.total_recombination_counts.maternal.p_0.5.no_double_recombinants_within_1000000bp.txt",header=T,stringsAsFactors=F)
#qtr.610.duohmm.all=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR_minus_MZs.610K.all.post_MERLIN.max_0.05_missing.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_and_clustered_duplicate_xovers_removed.txt",header=T,stringsAsFactors=F)
qtr.610.duohmm.all=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR_minus_MZs.610K.all.post_MERLIN.max_0.05_missing.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",header=T,stringsAsFactors=F)
#ntr.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/NTR/duoHMM_results/NTR.generr_removed.total_recombination_counts.maternal.p_0.9.txt",header=T,stringsAsFactors=F)
#ntr.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/NTR/duoHMM_results/NTR.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.maternal.p_0.5.txt",header=T,stringsAsFactors=F)
#ntr.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/NTR/duoHMM_results/NTR.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.maternal.p_0.5.no_double_recombinants_within_1000000bp.txt",header=T,stringsAsFactors=F)

#ntr.duohmm.all=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/NTR/duoHMM_results/NTR_minus_MZs.all.post_MERLIN.max_0.05_missing.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_and_clustered_duplicate_xovers_removed.txt",header=T,stringsAsFactors=F)
#ntr.duohmm.all=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/NTR/duoHMM_results/NTR_minus_MZs.all.post_MERLIN.max_0.05_missing.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",header=T,stringsAsFactors=F)
ntr.duohmm.all=read.delim("NTR/duoHMM_results/NTR_v2_minus_MZs.all.post_MERLIN.max_0.05_missing.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",header=T,stringsAsFactors=F)


qtr.duohmm.all=rbind(qtr.370.duohmm.all,qtr.610.duohmm.all)


colnames(qtr.370.duohmm.all)[c(1,7)]=c("child","total")
colnames(qtr.610.duohmm.all)[c(1,7)]=c("child","total")
colnames(ntr.duohmm.all)[c(1,7)]=c("child","total")

qtr.370.duohmm=qtr.370.duohmm.all[qtr.370.duohmm.all$sex=="Female",]
qtr.610.duohmm=qtr.610.duohmm.all[qtr.610.duohmm.all$sex=="Female",]
ntr.duohmm=ntr.duohmm.all[ntr.duohmm.all$sex=="Female",]

qtr.duohmm=rbind(qtr.370.duohmm,qtr.610.duohmm)


rownames(ntr.duohmm)=ntr.duohmm$child
rownames(qtr.duohmm)=qtr.duohmm$child
rownames(qtr.610.duohmm)=qtr.610.duohmm$child

qtr.mat.counts.by.child=cbind(qtr.duohmm[qtr.all.nftools$twin1,"total"],qtr.duohmm[qtr.all.nftools$twin2,"total"],qtr.duohmm[qtr.all.nftools$sib,"total"])
colnames(qtr.mat.counts.by.child)=c("twin1","twin2","sib")
rownames(qtr.mat.counts.by.child)=qtr.all.nftools$mother

qtr.610.mat.counts.by.child=cbind(qtr.610.duohmm[qtr.610.nftools$twin1,"total"],qtr.610.duohmm[qtr.610.nftools$twin2,"total"],qtr.610.duohmm[qtr.610.nftools$sib,"total"])
colnames(qtr.610.mat.counts.by.child)=c("twin1","twin2","sib")
rownames(qtr.610.mat.counts.by.child)=qtr.610.nftools$mother

qtr.370.mat.counts.by.child=qtr.mat.counts.by.child[!rownames(qtr.mat.counts.by.child) %in% rownames(qtr.610.mat.counts.by.child),]

ntr.mat.counts.by.child=cbind(ntr.duohmm[as.character(ntr.nftools$twin1),"total"],ntr.duohmm[as.character(ntr.nftools$twin2),"total"],ntr.duohmm[as.character(ntr.nftools$sib),"total"])
colnames(ntr.mat.counts.by.child)=c("twin1","twin2","sib")
rownames(ntr.mat.counts.by.child)=ntr.nftools$mother

combined.mat.counts.by.child=rbind(qtr.mat.counts.by.child,ntr.mat.counts.by.child)
combined.mat.counts.by.child.same.chip=rbind(qtr.610.mat.counts.by.child,ntr.mat.counts.by.child)

mat.comparison.pearson.qtr=compare.correlations.by.simulation(qtr.mat.counts.by.child,"pearson",1000)
mat.comparison.spearman.qtr=compare.correlations.by.simulation(qtr.mat.counts.by.child,"spearman",1000)
mat.comparison.pearson.ntr=compare.correlations.by.simulation(ntr.mat.counts.by.child,"pearson",1000)
mat.comparison.spearman.ntr=compare.correlations.by.simulation(ntr.mat.counts.by.child,"spearman",1000)
mat.comparison.pearson.qtr.610=compare.correlations.by.simulation(qtr.610.mat.counts.by.child,"pearson",1000)
mat.comparison.spearman.qtr.610=compare.correlations.by.simulation(qtr.610.mat.counts.by.child,"spearman",1000)
mat.comparison.pearson.qtr.370=compare.correlations.by.simulation(qtr.370.mat.counts.by.child,"pearson",1000)
mat.comparison.spearman.qtr.370=compare.correlations.by.simulation(qtr.370.mat.counts.by.child,"spearman",1000)
mat.comparison.pearson.combined=compare.correlations.by.simulation(combined.mat.counts.by.child,"pearson",1000)
mat.comparison.spearman.combined=compare.correlations.by.simulation(combined.mat.counts.by.child,"spearman",1000)
mat.comparison.pearson.combined=compare.correlations.by.simulation(combined.mat.counts.by.child,"pearson",1000)
mat.comparison.spearman.combined=compare.correlations.by.simulation(combined.mat.counts.by.child,"spearman",1000)
mat.comparison.pearson.combined.same.chip=compare.correlations.by.simulation(combined.mat.counts.by.child.same.chip,"pearson",1000)
mat.comparison.spearman.combined.same.chip=compare.correlations.by.simulation(combined.mat.counts.by.child.same.chip,"spearman",1000)

all.mat.comparisons=rbind(mat.comparison.pearson.qtr,mat.comparison.spearman.qtr,mat.comparison.pearson.qtr.610,mat.comparison.spearman.qtr.610,mat.comparison.pearson.qtr.370,mat.comparison.spearman.qtr.370,mat.comparison.pearson.ntr,mat.comparison.spearman.ntr,
mat.comparison.pearson.combined,mat.comparison.spearman.combined,mat.comparison.pearson.combined.same.chip,mat.comparison.spearman.combined.same.chip)
colnames(all.mat.comparisons)=c("correlation","p.val.cor.not.0","empirical.pval.twin.cor.gt.sib.cor")
rownames(all.mat.comparisons)=c("QTR.all.Pearson","QTR.all.Spearman","QTR.610.Pearson","QTR.610.Spearman","QTR.370.Pearson","QTR.370.Spearman","NTR.Pearson","NTR.Spearman","combined.Pearson","combined.Spearman","combined.same.chip.Pearson",
"combined.same.chip.Spearman")


qtr.610.nftools=read.delim("/home/hilary/TwinRecombination/Brisbane/NFTOOLS_paternal_counts_for_QTR_610K_families.twins_and_sibling.txt",header=T,stringsAsFactors=F)
qtr.all.nftools=read.delim("/home/hilary/TwinRecombination/Brisbane/NFTOOLS_paternal_counts_for_all_QTR_families_610_or_370_separately.twins_and_sibling.txt",header=T,stringsAsFactors=F)
ntr.nftools=read.delim("/home/hilary/TwinRecombination/Amsterdam/NFTOOLS_paternal_counts_for_all_NTR_families.twins_and_sibling.txt",header=T,stringsAsFactors=F)

qtr.370.nftools=qtr.all.nftools[!qtr.all.nftools$father %in% qtr.610.nftools$father,]
qtr.370.nftools.pat.counts.by.child=qtr.370.nftools[,c("twin1_pat_count","twin2_pat_count","sib_pat_count")]

pat.comparison.pearson.qtr370.nftools.pearson = compare.correlations.by.simulation(qtr.370.nftools.pat.counts.by.child,"pearson",1000)
pat.comparison.spearman.qtr370.nftools.spearman = compare.correlations.by.simulation(qtr.370.nftools.pat.counts.by.child,"spearman",1000)

nftools.qtr370.comparisons=rbind(mat.comparison.pearson.qtr370.nftools.pearson,mat.comparison.spearman.qtr370.nftools.spearman,pat.comparison.pearson.qtr370.nftools.pearson,pat.comparison.spearman.qtr370.nftools.spearman)
rownames(nftools.qtr370.comparisons)=c("mat.qtr370.nftools.pearson","mat.qtr370.nftools.spearman","pat.qtr370.nftools.pearson","pat.qtr370.nftools.spearman")

####Paternal
#write.table(nftools.qtr370.comparisons,"/home/hilary/maternal_age_recombination/comparing_correlation_in_crossover_count_between_twins_to_between_siblings.NFTOOLS.QTR370.txt",quote=F,sep="\t")


#qtr.370.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.370K.generr_removed.total_recombination_counts.paternal.p_0.9.txt",header=T,stringsAsFactors=F)
#qtr.370.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.370K.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.5.txt",header=T,stringsAsFactors=F)
#qtr.370.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.370K.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.5.no_double_recombinants_within_1000000bp.txt",header=T,stringsAsFactors=F)
#qtr.610.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.610K.generr_removed.total_recombination_counts.paternal.p_0.9.txt",header=T,stringsAsFactors=F)
#qtr.610.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.610K.generr_removed.total_recombination_counts.paternal.p_0.5.txt",header=T,stringsAsFactors=F)
#qtr.610.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/duoHMM_results/QTR.610K.generr_removed.total_recombination_counts.paternal.p_0.5.no_double_recombinants_within_1000000bp.txt",header=T,stringsAsFactors=F)
#ntr.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/NTR/duoHMM_results/NTR.generr_removed.total_recombination_counts.paternal.p_0.9.txt",header=T,stringsAsFactors=F)
#ntr.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/NTR/duoHMM_results/NTR.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.5.txt",header=T,stringsAsFactors=F)
#ntr.duohmm=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/NTR/duoHMM_results/NTR.post_MERLIN.max_0.05_missing.generr_removed.total_recombination_counts.paternal.p_0.5.no_double_recombinants_within_1000000bp.txt",header=T,stringsAsFactors=F)


qtr.370.duohmm=qtr.370.duohmm.all[qtr.370.duohmm.all$sex=="Male",]
qtr.610.duohmm=qtr.610.duohmm.all[qtr.610.duohmm.all$sex=="Male",]
ntr.duohmm=ntr.duohmm.all[ntr.duohmm.all$sex=="Male",]


qtr.duohmm=rbind(qtr.370.duohmm,qtr.610.duohmm)
rownames(ntr.duohmm)=ntr.duohmm$child
rownames(qtr.duohmm)=qtr.duohmm$child
rownames(qtr.610.duohmm)=qtr.610.duohmm$child

qtr.pat.counts.by.child=cbind(qtr.duohmm[qtr.all.nftools$twin1,"total"],qtr.duohmm[qtr.all.nftools$twin2,"total"],qtr.duohmm[qtr.all.nftools$sib,"total"])
colnames(qtr.pat.counts.by.child)=c("twin1","twin2","sib")
rownames(qtr.pat.counts.by.child)=qtr.all.nftools$father

qtr.610.pat.counts.by.child=cbind(qtr.610.duohmm[qtr.610.nftools$twin1,"total"],qtr.610.duohmm[qtr.610.nftools$twin2,"total"],qtr.610.duohmm[qtr.610.nftools$sib,"total"])
colnames(qtr.610.pat.counts.by.child)=c("twin1","twin2","sib")
rownames(qtr.610.pat.counts.by.child)=qtr.610.nftools$father

qtr.370.pat.counts.by.child=qtr.pat.counts.by.child[!rownames(qtr.pat.counts.by.child) %in% rownames(qtr.610.pat.counts.by.child),]

ntr.pat.counts.by.child=cbind(ntr.duohmm[as.character(ntr.nftools$twin1),"total"],ntr.duohmm[as.character(ntr.nftools$twin2),"total"],ntr.duohmm[as.character(ntr.nftools$sib),"total"])
colnames(ntr.pat.counts.by.child)=c("twin1","twin2","sib")
rownames(ntr.pat.counts.by.child)=ntr.nftools$father


combined.pat.counts.by.child=rbind(qtr.pat.counts.by.child,ntr.pat.counts.by.child)
combined.pat.counts.by.child.same.chip=rbind(qtr.610.pat.counts.by.child,ntr.pat.counts.by.child)

pat.comparison.pearson.qtr=compare.correlations.by.simulation(qtr.pat.counts.by.child,"pearson",1000)
pat.comparison.spearman.qtr=compare.correlations.by.simulation(qtr.pat.counts.by.child,"spearman",1000)
pat.comparison.pearson.ntr=compare.correlations.by.simulation(ntr.pat.counts.by.child,"pearson",1000)
pat.comparison.spearman.ntr=compare.correlations.by.simulation(ntr.pat.counts.by.child,"spearman",1000)
pat.comparison.pearson.qtr.610=compare.correlations.by.simulation(qtr.610.pat.counts.by.child,"pearson",1000)
pat.comparison.spearman.qtr.610=compare.correlations.by.simulation(qtr.610.pat.counts.by.child,"spearman",1000)
pat.comparison.pearson.qtr.370=compare.correlations.by.simulation(qtr.370.pat.counts.by.child,"pearson",1000)
pat.comparison.spearman.qtr.370=compare.correlations.by.simulation(qtr.370.pat.counts.by.child,"spearman",1000)
pat.comparison.pearson.combined=compare.correlations.by.simulation(combined.pat.counts.by.child,"pearson",1000)
pat.comparison.spearman.combined=compare.correlations.by.simulation(combined.pat.counts.by.child,"spearman",1000)
pat.comparison.pearson.combined=compare.correlations.by.simulation(combined.pat.counts.by.child,"pearson",1000)
pat.comparison.spearman.combined=compare.correlations.by.simulation(combined.pat.counts.by.child,"spearman",1000)
pat.comparison.pearson.combined.same.chip=compare.correlations.by.simulation(combined.pat.counts.by.child.same.chip,"pearson",1000)
pat.comparison.spearman.combined.same.chip=compare.correlations.by.simulation(combined.pat.counts.by.child.same.chip,"spearman",1000)

all.pat.comparisons=rbind(pat.comparison.pearson.qtr,pat.comparison.spearman.qtr,pat.comparison.pearson.qtr.610,pat.comparison.spearman.qtr.610,pat.comparison.pearson.qtr.370,
    pat.comparison.spearman.qtr.370,pat.comparison.pearson.ntr,pat.comparison.spearman.ntr,
pat.comparison.pearson.combined,pat.comparison.spearman.combined,pat.comparison.pearson.combined.same.chip,pat.comparison.spearman.combined.same.chip)
colnames(all.pat.comparisons)=c("correlation","p.val.cor.not.0","empirical.pval.twin.cor.gt.sib.cor")
rownames(all.pat.comparisons)=c("QTR.all.Pearson","QTR.all.Spearman","QTR.610.Pearson","QTR.610.Spearman","QTR.370.Pearson","QTR.370.Spearman","NTR.Pearson","NTR.Spearman","combined.Pearson",
            "combined.Spearman","combined.same.chip.Pearson","combined.same.chip.Spearman")

colnames(all.mat.comparisons)=paste0("maternal_",colnames(all.mat.comparisons))
colnames(all.pat.comparisons)=paste0("paternal_",colnames(all.pat.comparisons))
all.comparisons=cbind(all.mat.comparisons,all.pat.comparisons)
all.comparisons=cbind(all.mat.comparisons,all.pat.comparisons)
#write.table(all.comparisons,"/home/hilary/maternal_age_recombination/comparing_correlation_in_crossover_count_between_twins_to_between_siblings.duoHMM_p_0.9.txt",quote=F,sep="\t")
#write.table(all.comparisons,"/home/hilary/maternal_age_recombination/comparing_correlation_in_crossover_count_between_twins_to_between_siblings.duoHMM_p_0.5.QTR370_and_NTR_MERLIN_error_corrected.txt",quote=F,sep="\t")
#write.table(all.comparisons,"/home/hilary/maternal_age_recombination/comparing_correlation_in_crossover_count_between_twins_to_between_siblings.duoHMM_p_0.5.no_double_recombinants_within_1000000bp.QTR370_and_NTR_MERLIN_error_corrected.txt",quote=F,sep="\t")

#write.table(all.comparisons,"/home/hilary/maternal_age_recombination/comparing_correlation_in_crossover_count_between_twins_to_between_siblings.duoHMM_p_0.5.duoHMM_new_version.post_MERLIN.generr_removed.double_xovers_within_X_SNPs_and_clustered_duplicate_xovers_removed.txt",quote=F,sep="\t")
#write.table(all.comparisons,"/home/hilary/maternal_age_recombination/comparing_correlation_in_crossover_count_between_twins_to_between_siblings.duoHMM_p_0.5.duoHMM_new_version.post_MERLIN.generr_removed.double_xovers_within_X_SNPs.same_X_both_sexes.txt",quote=F,sep="\t")
write.table(all.comparisons,"/home/hilary/maternal_age_recombination/comparing_correlation_in_crossover_count_between_twins_to_between_siblings.duoHMM_p_0.5.duoHMM_new_version.post_MERLIN.generr_removed.double_xovers_within_X_SNPs.same_X_both_sexes.NTR_v2.txt",quote=F,sep="\t")
