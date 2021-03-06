count.files=c("CARL/duoHMM_results/CARL.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"FC/duoHMM_results/FC.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"FVG/duoHMM_results/FVG.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"NTR/duoHMM_results/NTR.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"ORCADES/duoHMM_results/ORCADES.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"QTR/duoHMM_results/all_QTR_370K_samples.all.post_MERLIN.more_stringent.generr_removed.recombination_count.all_duos.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"QTR/duoHMM_results/all_QTR_610K_samples.all.post_MERLIN.more_stringent.generr_removed.recombination_count.all_duos.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"QTR/duoHMM_results/new_samples.CoreExome.all.post_MERLIN.more_stringent.generr_removed.recombination_count.all_duos.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"VB/duoHMM_results/VB.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"VIS_KORCULA/duoHMM_results/KORCULA.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt",
"VIS_KORCULA/duoHMM_results/VIS.all.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.txt")
count.files=gsub("double_xovers_within_X_SNPs_removed.","",count.files)

age.files=c("CARL/clean-carl-flip_b37-related.with_parental_age.txt",
    "FC/FC.with_parental_age.txt",
    "FVG/clean-fvg-flip_b37-related.with_parental_age.txt",
    "NTR/NTR_v2_minus_MZs.with_parental_age.txt",
    "ORCADES/ORCADES.with_parental_age.txt",
    "QTR/all_QTR_370K_samples.with_parental_age.txt",
    "QTR/all_QTR_610K_samples.with_parental_age.txt",
    "QTR/new_samples.CoreExome.with_parental_age.txt",
    "VB/valborbera_b37-related.with_parental_age.txt",
    "VIS_KORCULA/KORCULA.with_parental_age.txt",
    "VIS_KORCULA/VIS.with_parental_age.txt")


names(age.files)=c("CARL","FC","FVG","NTR","ORCADES","QTR370","QTR610","QTRCoreExome","VB","KORCULA","VIS")
names(count.files)=names(age.files)

removed = NULL

family.type.count=NULL
all.data0=list()

for(i in 1:length(count.files)){
    counts=read.delim(count.files[i],header=T,stringsAsFactors=F,colClasses=c("character","character","logical","character","character","integer","integer",rep("logical",11)))
    age=read.delim(age.files[i],header=T,stringsAsFactors=F)#,colClasses=c(rep("character",4),rep("numeric",4)))

    age$maternal_age=as.numeric(age$maternal_age)
    age$paternal_age=as.numeric(age$paternal_age)
    if(names(count.files)[i]=="FC"){
        age[,2] =gsub("-","_",age[,2])
    }
    counts=counts[!is.na(counts$sex),]	
    mat.counts=counts[counts$sex=="Female",]
    pat.counts=counts[counts$sex=="Male",]
    mat.counts.2=unique(mat.counts[,c("CHILD","PARENT","informative","informative.2gen.2parents","noninf.2kids.2gen.2parents")])
    pat.counts.2=unique(pat.counts[,c("CHILD","PARENT","informative","informative.2gen.2parents","noninf.2kids.2gen.2parents")])
    family.type.count=rbind(family.type.count,c(length(unique(mat.counts.2$PARENT[mat.counts.2$informative])),length(unique(mat.counts.2$CHILD[mat.counts.2$informative])),
        length(unique(mat.counts.2$PARENT[mat.counts.2$informative.2gen.2parents])),length(unique(mat.counts.2$CHILD[mat.counts.2$informative.2gen.2parents])),
        length(unique(mat.counts.2$PARENT[mat.counts.2$noninf.2kids.2gen.2parents])),length(unique(mat.counts.2$CHILD[mat.counts.2$noninf.2kids.2gen.2parents])),
        length(unique(pat.counts.2$PARENT[pat.counts.2$informative])),length(unique(pat.counts.2$CHILD[pat.counts.2$informative])),
        length(unique(pat.counts.2$PARENT[pat.counts.2$informative.2gen.2parents])),length(unique(pat.counts.2$CHILD[pat.counts.2$informative.2gen.2parents])),
        length(unique(pat.counts.2$PARENT[pat.counts.2$noninf.2kids.2gen.2parents])),length(unique(pat.counts.2$CHILD[pat.counts.2$noninf.2kids.2gen.2parents]))))    
    mat.counts$age.at.birth = age[match(mat.counts$CHILD,age$individual),"maternal_age"]
    pat.counts$age.at.birth = age[match(pat.counts$CHILD,age$individual),"paternal_age"]
    cat(names(count.files)[i],"\t",nrow(mat.counts),"\t",nrow(pat.counts),"\n")
#remove any with dodgy ages, or missing ages
	removed=rbind(removed,mat.counts[is.na(mat.counts$age.at.birth) |(!is.na(mat.counts$age.at.birth) & !( mat.counts$age.at.birth>12 & mat.counts$age.at.birth<50)),],
	pat.counts[is.na(pat.counts$age.at.birth) |(!is.na(pat.counts$age.at.birth) &!pat.counts$age.at.birth>12),]) 
    mat.counts=mat.counts[!is.na(mat.counts$age.at.birth) & mat.counts$age.at.birth>12 & mat.counts$age.at.birth<50,]
    pat.counts=pat.counts[!is.na(pat.counts$age.at.birth) & pat.counts$age.at.birth>12,]
    cat(names(count.files)[i],"\t",nrow(mat.counts),"\t",nrow(pat.counts),"\n")
    counts=rbind(mat.counts,pat.counts)
    counts$cohort=names(count.files)[i]
     all.data0[[i]]=counts
}
names(all.data0) = names(count.files)

qtr370 = all.data0[["QTR370"]]
qtr370$cohort = "QTR370"
qtr610 = all.data0[["QTR610"]]
qtr610$cohort = "QTR610"
qtrCoreExome = all.data0[["QTRCoreExome"]]
qtrCoreExome$cohort = "QTRCoreExome"

qtr370.keep  = qtr370[!qtr370$PARENT %in% qtr370[qtr370$duo %in% qtr610$duo,"PARENT"],]
qtrCoreExome.keep = qtrCoreExome[!qtrCoreExome$PARENT %in% qtrCoreExome[qtrCoreExome$duo %in% qtr610$duo|qtrCoreExome$duo %in% qtr370$duo ,"PARENT"],]

qtr370.chuck  = qtr370[qtr370$PARENT %in% qtr370[qtr370$duo %in% qtr610$duo,"PARENT"],]
qtrCoreExome.chuck = qtrCoreExome[qtrCoreExome$PARENT %in% qtrCoreExome[qtrCoreExome$duo %in% qtr610$duo|qtrCoreExome$duo %in% qtr370$duo ,"PARENT"],]

all.data0[["QTR370"]]=qtr370.keep
all.data0[["QTRCoreExome"]]=qtrCoreExome.keep

write.table(rbind(qtr370.chuck,qtrCoreExome.chuck),
"QTR/duoHMM_results/duos_removed_from_QTR370_and_QTRCoreExome_due_to_overlap_with_other_chips.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t")


write.table(removed,"samples_removed_due_to_missing_or_suspicious_age_data.more_stringent.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t")


rownames(family.type.count)=names(age.files)
colnames(family.type.count)= c("mother.inf","mother.inf.n.kids","mother.inf.min3kids","mother.inf.min3kids.n.kids","mother.noninf.2kids","mother.noninf.2kids.n.kids","father.inf",
"father.inf.n.kids","father.inf.min3kids",            "father.inf.min3kids.n.kids",            "father.noninf.2kids","father.noninf.2kids.n.kids")

#write.table(family.type.count,"number_of_family_types_by_cohort.more_stringent.txt",quote=F,sep="\t")


all.data = do.call("rbind",all.data0)


write.table(all.data,"all_cohorts.post_MERLIN.generr_removed.recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t",row.names=F)

# for meioses from informative families (both parents + >2 kids)
data0=all.data[all.data$informative.2gen.2parents,]

data0.mat=data0[data0$sex=="Female",]
data0.pat=data0[data0$sex=="Male",]
mat.cohort.count0 = table(data0.mat$cohort)
pat.cohort.count0 = table(data0.pat$cohort)

data0.unique.mat=unique(data0.mat[,c("PARENT","cohort")])
data0.unique.pat=unique(data0.pat[,c("PARENT","cohort")])

mat.cohort.count.family = table(data0.unique.mat$cohort)
pat.cohort.count.family = table(data0.unique.pat$cohort)
family.type.count.inf.nuc=rbind(mat.cohort.count0,mat.cohort.count.family,pat.cohort.count0,pat.cohort.count.family)

data0.mat=data0.mat[data0.mat$cohort %in% names(mat.cohort.count0[mat.cohort.count0>20]),]
data0.pat=data0.pat[data0.pat$cohort %in% names(pat.cohort.count0[pat.cohort.count0>20]),]


mat.key.model0=unique(data.frame(as.integer(as.factor(data0.mat$cohort)),as.factor(data0.mat$cohort)))
pat.key.model0=unique(data.frame(as.integer(as.factor(data0.pat$cohort)),as.factor(data0.pat$cohort)))
colnames(mat.key.model0)=c("Cohort.type","Cohort")
colnames(pat.key.model0)=c("Cohort.type","Cohort")
write.table(mat.key.model0,"RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.informative_nuclear_families_only.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t",row.names=F)
write.table(pat.key.model0,"RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.informative_nuclear_families_only.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t",row.names=F)

data0.mat2=list(y=data0.mat$nrec,cohort=as.integer(as.factor(data0.mat$cohort)),family=as.integer(as.factor(data0.mat$PARENT)),Age=as.numeric(data0.mat$age.at.birth),J=nrow(data0.mat),
    I=length(unique(data0.mat$PARENT)),    C=length(unique(data0.mat$cohort)))
data0.pat2=list(y=data0.pat$nrec,cohort=as.integer(as.factor(data0.pat$cohort)),family=as.integer(as.factor(data0.pat$PARENT)),Age=as.numeric(data0.pat$age.at.birth),J=nrow(data0.pat),
    I=length(unique(data0.pat$PARENT)),    C=length(unique(data0.pat$cohort)))
cohort_by_family=rep(NA,data0.mat2$I)
for(i in 1:data0.mat2$I){
    cohort_by_family[i] <- data0.mat2$cohort[data0.mat2$family ==i][1]
}
data0.mat2$cohort_by_family = cohort_by_family

cohort_by_family=rep(NA,data0.pat2$I)
for(i in 1:data0.pat2$I){
    cohort_by_family[i] <- data0.pat2$cohort[data0.pat2$family ==i][1]
}
data0.pat2$cohort_by_family = cohort_by_family

#for informative meioses only with nkid>1, all cohorts, fit cohort effect and family effect
data1=all.data[all.data$nkid>1 & all.data$informative,]

data1.mat=data1[data1$sex=="Female",]
data1.pat=data1[data1$sex=="Male",]

mat.cohort.count = table(data1.mat$cohort)
pat.cohort.count = table(data1.pat$cohort)

data1.unique.mat=unique(data1.mat[,c("PARENT","cohort")])
data1.unique.pat=unique(data1.pat[,c("PARENT","cohort")])

mat.cohort.count.family1 = table(data1.unique.mat$cohort)
pat.cohort.count.family1 = table(data1.unique.pat$cohort)

family.type.count.inf=rbind(mat.cohort.count,mat.cohort.count.family1,pat.cohort.count,pat.cohort.count.family1)


                                        #only keep cohorts with more than 20
data1.mat=data1.mat[data1.mat$cohort %in% names(mat.cohort.count[mat.cohort.count>20]),]
data1.pat=data1.pat[data1.pat$cohort %in% names(pat.cohort.count[pat.cohort.count>20]),]

mat.key.model1=unique(data.frame(as.integer(as.factor(data1.mat$cohort)),as.factor(data1.mat$cohort)))
pat.key.model1=unique(data.frame(as.integer(as.factor(data1.pat$cohort)),as.factor(data1.pat$cohort)))
colnames(mat.key.model1)=c("Cohort.type","Cohort")
colnames(pat.key.model1)=c("Cohort.type","Cohort")
write.table(mat.key.model1,"RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t",row.names=F)
write.table(pat.key.model1,"RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t",row.names=F)

data1.mat2=list(y=data1.mat$nrec,cohort=as.integer(as.factor(data1.mat$cohort)),family=as.integer(as.factor(data1.mat$PARENT)),Age=as.numeric(data1.mat$age.at.birth),J=nrow(data1.mat),I=length(unique(data1.mat$PARENT)),    C=length(unique(data1.mat$cohort)))
data1.pat2=list(y=data1.pat$nrec,cohort=as.integer(as.factor(data1.pat$cohort)),family=as.integer(as.factor(data1.pat$PARENT)),Age=as.numeric(data1.pat$age.at.birth),J=nrow(data1.pat),I=length(unique(data1.pat$PARENT)),    C=length(unique(data1.pat$cohort)))
cohort_by_family=rep(NA,data1.mat2$I)
for(i in 1:data1.mat2$I){
    cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
}
data1.mat2$cohort_by_family = cohort_by_family
cohort_by_family=rep(NA,data1.pat2$I)
for(i in 1:data1.pat2$I){
    cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
}
data1.pat2$cohort_by_family = cohort_by_family

###look at distribution of mean per mother for cohort

data1.mat.by.cohort=split(data1.mat,data1.mat$cohort)
mat.mean.by.parent = list()
for(i in 1:length(data1.mat.by.cohort)){
    counts=data1.mat.by.cohort[[i]] 
    counts$age.at.birth=as.numeric(counts$age.at.birth)
    parents=unique(counts$PARENT)
    mean.by.parents=sapply(parents,function(p){mean(counts[counts$PARENT==p,"nrec"])})
    mean.age.by.parents=sapply(parents,function(p){mean(as.numeric(counts[counts$PARENT==p,"age.at.birth"]))})
    mat.mean.by.parent[[i]]=cbind(mean.by.parents,mean.age.by.parents)
    counts$adjusted.nrec=counts$nrec-mean.by.parents[counts$PARENT]
    counts$adjusted.age=counts$age.at.birth-mean.age.by.parents[counts$PARENT]
    n.children.per.parent=sapply(parents,function(p){sum(counts$PARENT==p)})
    counts$nkid2=n.children.per.parent[counts$PARENT]
    data1.mat.by.cohort[[i]]=counts
}
names(mat.mean.by.parent)=names(data1.mat.by.cohort)

mat.variance.in.mean.by.family=unlist(lapply(mat.mean.by.parent,function(x){var(x[,1])}))

data1.pat.by.cohort=split(data1.pat,data1.pat$cohort)
pat.mean.by.parent = list()
for(i in 1:length(data1.pat.by.cohort)){
    counts=data1.pat.by.cohort[[i]] 
    counts$age.at.birth=as.numeric(counts$age.at.birth)
    parents=unique(counts$PARENT)
    mean.by.parents=sapply(parents,function(p){mean(counts[counts$PARENT==p,"nrec"])})
    mean.age.by.parents=sapply(parents,function(p){mean(as.numeric(counts[counts$PARENT==p,"age.at.birth"]))})
    pat.mean.by.parent[[i]]=cbind(mean.by.parents,mean.age.by.parents)
    counts$adjusted.nrec=counts$nrec-mean.by.parents[counts$PARENT]
    counts$adjusted.age=counts$age.at.birth-mean.age.by.parents[counts$PARENT]
    n.children.per.parent=sapply(parents,function(p){sum(counts$PARENT==p)})
    counts$nkid2=n.children.per.parent[counts$PARENT]
    data1.pat.by.cohort[[i]]=counts

}
names(pat.mean.by.parent)=names(data1.pat.by.cohort)
pat.variance.in.mean.by.family=unlist(lapply(pat.mean.by.parent,function(x){var(x[,1])}))

data1.adjusted.mat= do.call("rbind",data1.mat.by.cohort)
data1.adjusted.pat= do.call("rbind",data1.pat.by.cohort)

data1.adjusted.mat2=list(y=data1.adjusted.mat$adjusted.nrec,cohort=as.integer(as.factor(data1.adjusted.mat$cohort)),family=as.integer(as.factor(data1.adjusted.mat$PARENT)),Age=as.numeric(data1.adjusted.mat$adjusted.age),
    J=nrow(data1.adjusted.mat),I=length(unique(data1.adjusted.mat$PARENT)), C=length(unique(data1.adjusted.mat$cohort)),nkids=data1.adjusted.mat$nkid)

data1.adjusted.pat2=list(y=data1.adjusted.pat$adjusted.nrec,cohort=as.integer(as.factor(data1.adjusted.pat$cohort)),family=as.integer(as.factor(data1.adjusted.pat$PARENT)),Age=as.numeric(data1.adjusted.pat$adjusted.age),
    J=nrow(data1.adjusted.pat),I=length(unique(data1.adjusted.pat$PARENT)), C=length(unique(data1.adjusted.pat$cohort)),nkids=data1.adjusted.pat$nkid)



variance.in.mean.by.family=cbind(mat.variance.in.mean.by.family,pat.variance.in.mean.by.family[match(names(mat.variance.in.mean.by.family),names(pat.variance.in.mean.by.family))])
colnames(variance.in.mean.by.family)=c("maternal","paternal")

write.table(variance.in.mean.by.family,"variance_in_mean_duoHMM_count_by_parent.informative_meioses_only.min_2_kids.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t")
##### Model 2
##########for informative meioses only, all cohorts, fit cohort effect, BUT NO FAMILY EFFECT

data1.2=all.data[all.data$informative,]

data1.2.mat=data1.2[data1.2$sex=="Female",]
data1.2.pat=data1.2[data1.2$sex=="Male",]

mat.cohort.count = table(data1.2.mat$cohort)
pat.cohort.count = table(data1.2.pat$cohort)
#only keep cohorts with more than 20
data1.2.mat=data1.2.mat[data1.2.mat$cohort %in% names(mat.cohort.count[mat.cohort.count>20]),]
data1.2.pat=data1.2.pat[data1.2.pat$cohort %in% names(pat.cohort.count[pat.cohort.count>20]),]

mat.key.model2=unique(data.frame(as.integer(as.factor(data1.2.mat$cohort)),as.factor(data1.2.mat$cohort)))
pat.key.model2=unique(data.frame(as.integer(as.factor(data1.2.pat$cohort)),as.factor(data1.2.pat$cohort)))
colnames(mat.key.model2)=c("Cohort.type","Cohort")
colnames(pat.key.model2)=c("Cohort.type","Cohort")

write.table(mat.key.model2,"RSTAN_output/key_for_maternal_cohorts_to_include_in_model_2.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t",row.names=F)
write.table(pat.key.model2,"RSTAN_output/key_for_paternal_cohorts_to_include_in_model_2.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t",row.names=F)

data1.2.mat2=list(y=data1.2.mat$nrec,cohort=as.integer(as.factor(data1.2.mat$cohort)),family=as.integer(as.factor(data1.2.mat$PARENT)),Age=as.numeric(data1.2.mat$age.at.birth),J=nrow(data1.2.mat),
    I=length(unique(data1.2.mat$PARENT)),    C=length(unique(data1.2.mat$cohort)))
data1.2.pat2=list(y=data1.2.pat$nrec,cohort=as.integer(as.factor(data1.2.pat$cohort)),family=as.integer(as.factor(data1.2.pat$PARENT)),Age=as.numeric(data1.2.pat$age.at.birth),J=nrow(data1.2.pat),
    I=length(unique(data1.2.pat$PARENT)),    C=length(unique(data1.2.pat$cohort)))

#### Model 3
############for nkid=>2, family effect, and fit N_called_crossovers~Binom(Y,p=f(cohort*family_type)); only for those with 

data2=all.data[all.data$nkid>1 & (all.data$informative.2gen.2parents|all.data$informative.2gen.1parent|all.data$noninf.2kids.2gen.2parents|all.data$noninf.2kids.2gen.1parent
|all.data$informative.3gen.2parents|all.data$informative.3gen.1parent),]

data2$cohort.family.type = data2$cohort
data2$cohort.family.type[data2$informative.2gen.2parents]=paste0(data2$cohort[data2$informative.2gen.2parents],".infor.2gen.2parents")
data2$cohort.family.type[data2$informative.2gen.1parent]=paste0(data2$cohort[data2$informative.2gen.1parent],".infor.2gen.1parent")
data2$cohort.family.type[data2$informative.3gen.2parents]=paste0(data2$cohort[data2$informative.3gen.2parents],".infor.3gen.2parents")
data2$cohort.family.type[data2$informative.3gen.1parent]=paste0(data2$cohort[data2$informative.3gen.1parent],".infor.3gen.1parent")
data2$cohort.family.type[data2$noninf.2kids.2gen.2parents]=paste0(data2$cohort[data2$noninf.2kids.2gen.2parents],".noninfor.2kids.2gen.2parents")
data2$cohort.family.type[data2$noninf.2kids.2gen.1parent]=paste0(data2$cohort[data2$noninf.2kids.2gen.1parent],".noninfor.2kids.2gen.1parent")

data2.mat=data2[data2$sex=="Female",]
data2.pat=data2[data2$sex=="Male",]

mat.cohort.count = table(data2.mat$cohort.family.type)
pat.cohort.count = table(data2.pat$cohort.family.type)


data2.unique.mat=unique(data2.mat[,c("PARENT","cohort")])
data2.unique.pat=unique(data2.pat[,c("PARENT","cohort")])



mat.cohort.count2 = table(data2.mat$cohort)
pat.cohort.count2 = table(data2.pat$cohort)
mat.cohort.count.family2 = table(data2.unique.mat$cohort)
pat.cohort.count.family2 = table(data2.unique.pat$cohort)

family.type.count.multikid=data.frame(rbind(mat.cohort.count2,mat.cohort.count.family2,pat.cohort.count2,pat.cohort.count.family2))
                                        #only keep cohorts with more than 20
data2.mat=data2.mat[data2.mat$cohort.family.type %in% names(mat.cohort.count[mat.cohort.count>20]),]
data2.pat=data2.pat[data2.pat$cohort.family.type %in% names(pat.cohort.count[pat.cohort.count>20]),]

data2.unique.mat=unique(data2.mat[,c("PARENT","cohort")])
data2.unique.pat=unique(data2.pat[,c("PARENT","cohort")])


mat.cohort.count2.no.small = table(data2.mat$cohort)
pat.cohort.count2.no.small  = table(data2.pat$cohort)

mat.cohort.count.family2.no.small  = table(data2.unique.mat$cohort)
pat.cohort.count.family2.no.small  = table(data2.unique.pat$cohort)

mat.family.type.count.multikid.no.small =rbind(mat.cohort.count2.no.small ,mat.cohort.count.family2.no.small)
pat.family.type.count.multikid.no.small =rbind(pat.cohort.count2.no.small ,pat.cohort.count.family2.no.small )
rowSums(mat.family.type.count.multikid.no.small)
rowSums(pat.family.type.count.multikid.no.small)

rownames(family.type.count.multikid) = c("mat.n.kids","mat.n.fams","pat.n.kids","pat.n.fams")
rownames(family.type.count.inf) = c("mat.n.kids","mat.n.fams","pat.n.kids","pat.n.fams")
rownames(family.type.count.inf.nuc) = c("mat.n.kids","mat.n.fams","pat.n.kids","pat.n.fams")

family.type.count.multikid$total=rowSums(family.type.count.multikid)
family.type.count.multikid$total.without.gpc = rowSums(family.type.count.multikid[,c("CARL" ,"FC", "FVG","KORCULA","NTR","ORCADES","QTR370","QTR610","VB","VIS")])

family.type.count.multikid$total.no.small=c(rowSums(mat.family.type.count.multikid.no.small),rowSums(pat.family.type.count.multikid.no.small))
family.type.count.multikid$total.no.small.without.gpc = c(rowSums(mat.family.type.count.multikid.no.small[,c("CARL" ,"FC", "FVG","NTR","ORCADES","QTR370","QTR610","VB")]),
  rowSums(pat.family.type.count.multikid.no.small[,c("FC", "FVG","NTR","ORCADES","QTR370","QTR610","VB")]))

write.table(family.type.count.multikid,"number_of_duos_by_cohorts.families_with_multiple_children.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t")
write.table(family.type.count.inf,"number_of_duos_by_cohorts.informative_duos.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t")
write.table(family.type.count.inf.nuc,"number_of_duos_by_cohorts.informative_nuc_fam_duos.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t")


mat.key.model3=unique(data.frame(as.integer(as.factor(data2.mat$cohort.family.type)),as.factor(data2.mat$cohort.family.type)))
pat.key.model3=unique(data.frame(as.integer(as.factor(data2.pat$cohort.family.type)),as.factor(data2.pat$cohort.family.type)))
colnames(mat.key.model3)=c("Cohort.type","Cohort")
colnames(pat.key.model3)=c("Cohort.type","Cohort")

write.table(mat.key.model3,"RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t",row.names=F)
write.table(pat.key.model3,"RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",quote=F,sep="\t",row.names=F)


mat.key2.model3=unique(data.frame(as.integer(as.factor(data2.mat$cohort)),as.factor(data2.mat$cohort)))
pat.key2.model3=unique(data.frame(as.integer(as.factor(data2.pat$cohort)),as.factor(data2.pat$cohort)))
colnames(mat.key2.model3)=c("Cohort.type","Cohort")
colnames(pat.key2.model3)=c("Cohort.type","Cohort")

write.table(mat.key2.model3,"RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.more_stringent.with_new_QTR.no_GPC.raw_counts.original_cohorts.txt",quote=F,sep="\t",row.names=F)
write.table(pat.key2.model3,"RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.more_stringent.with_new_QTR.no_GPC.raw_counts.original_cohorts.txt",quote=F,sep="\t",row.names=F)



##need to alter this to include ORCADES
mat.beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.more_stringent.with_new_QTR.txt",header=T,stringsAsFactors=F)[,1:5]
pat.beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.more_stringent.with_new_QTR.txt",header=T,stringsAsFactors=F)[,1:5]

mat.key.model3=unique(data.frame(as.integer(as.factor(data2.mat$cohort.family.type)),as.factor(data2.mat$cohort.family.type)))
pat.key.model3=unique(data.frame(as.integer(as.factor(data2.pat$cohort.family.type)),as.factor(data2.pat$cohort.family.type)))
colnames(mat.key.model3)=c("Cohort.type","Cohort")
colnames(pat.key.model3)=c("Cohort.type","Cohort")

mat.beta.parameters$Cohort.type=mat.key.model3[match(mat.beta.parameters$Cohort,mat.key.model3$Cohort,),"Cohort.type"]
pat.beta.parameters$Cohort.type=pat.key.model3[match(pat.beta.parameters$Cohort,pat.key.model3$Cohort),"Cohort.type"]
mat.beta.parameters=mat.beta.parameters[!is.na(mat.beta.parameters$Cohort.type),]
pat.beta.parameters=pat.beta.parameters[!is.na(pat.beta.parameters$Cohort.type),]

write.table(mat.beta.parameters,"parameters_for_beta_distribution.model3_maternal.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",col.names=F,row.names=F,sep="\t",quote=F)
write.table(pat.beta.parameters,"parameters_for_beta_distribution.model3_paternal.more_stringent.with_new_QTR.no_GPC.raw_counts.txt",col.names=F,row.names=F,sep="\t",quote=F)

data2.mat2=list(y=data2.mat$nrec,cohort=as.integer(as.factor(data2.mat$cohort.family.type)),family=as.integer(as.factor(data2.mat$PARENT)),Age=as.numeric(data2.mat$age.at.birth),J=nrow(data2.mat),
    I=length(unique(data2.mat$PARENT)),  C=length(unique(data2.mat$cohort.family.type)),alpha_parameters=mat.beta.parameters[order(mat.beta.parameters$Cohort.type),"alpha"],
    beta_parameters=c(mat.beta.parameters[order(mat.beta.parameters$Cohort.type),"beta"]),cohort2=as.integer(as.factor(data2.mat$cohort)),C2=length(unique(data2.mat$cohort)))
data2.pat2=list(y=data2.pat$nrec,cohort=as.integer(as.factor(data2.pat$cohort.family.type)),family=as.integer(as.factor(data2.pat$PARENT)),Age=as.numeric(data2.pat$age.at.birth),J=nrow(data2.pat),
    I=length(unique(data2.pat$PARENT)),    C=length(unique(data2.pat$cohort.family.type)),alpha_parameters=pat.beta.parameters[order(pat.beta.parameters$Cohort.type),"alpha"],
    beta_parameters=pat.beta.parameters[order(pat.beta.parameters$Cohort.type),"beta"],cohort2=as.integer(as.factor(data2.pat$cohort)),C2=length(unique(data2.pat$cohort)))


save.image("duoHMM_data_for_RSTAN.more_stringent.with_new_QTR.no_GPC.raw_counts.RData")
