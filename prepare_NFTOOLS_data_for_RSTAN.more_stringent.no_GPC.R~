load("NFTOOLS_results/NFTOOLS_counts.all_families_except_outliers.more_stringent.RData")

age.files=c("NTR/NTR_v2_minus_MZs.with_parental_age.txt",
    "QTR/QTR_minus_MZs.370K.with_parental_age.txt",
    "QTR/QTR_minus_MZs.610K.with_parental_age.txt",
    "FC/FC.with_parental_age.txt",
    "GPC/GPC.with_parental_age.txt",
    "VB/valborbera_b37-related.with_parental_age.txt","ORCADES/ORCADES.with_parental_age.txt"            )

names(age.files)=names(all.events)[1:7]

data1.mat=NULL
data1.pat=NULL
mat.mean.by.parent=list()
pat.mean.by.parent=list()
for(i in 1:length(all.events)){
    events=all.events[[i]]
    events$cohort=names(all.events)[i]
    events$family=paste0(events$cohort,"_",events$family)
    mat.events=events[events$sex=="female" & events$nkids>2,]
    pat.events=events[events$sex=="male" & events$nkids>2,]
    age=read.delim(age.files[names(all.events)[i]],header=T,stringsAsFactors=F)    
    age$individual=gsub("-","_",age$individual)
    mat.events$age.at.birth = age[match(mat.events$child,age$individual),"maternal_age"]
    pat.events$age.at.birth = age[match(pat.events$child,age$individual),"paternal_age"]

#remove any children with missing parental ages
    mat.events=mat.events[!is.na(mat.events$age.at.birth),]
    pat.events=pat.events[!is.na(pat.events$age.at.birth),]
    
#remove any families with fewer than 2 children remaining
    mat.family.count=table(mat.events$family)
    pat.family.count=table(pat.events$family)
    print(nrow(mat.events))
    print(nrow(pat.events))
    mat.events=mat.events[mat.events$family %in% names(mat.family.count[mat.family.count>1]),]
    pat.events=pat.events[pat.events$family %in% names(pat.family.count[pat.family.count>1]),]
    print(nrow(mat.events))
    print(nrow(pat.events))

        
    mat.parents=unique(mat.events$family)
    mat.mean.by.parents=sapply(mat.parents,function(p){mean(mat.events[mat.events$family==p,"Freq"])})
    mat.mean.age.by.parents=sapply(mat.parents,function(p){mean(mat.events[mat.events$family==p,"age.at.birth"])})
    mat.mean.by.parent[[i]]=cbind(mat.mean.by.parents,mat.mean.age.by.parents)

    pat.parents=unique(pat.events$family)
    pat.mean.by.parents=sapply(pat.parents,function(p){mean(pat.events[pat.events$family==p,"Freq"])})
    pat.mean.age.by.parents=sapply(pat.parents,function(p){mean(pat.events[pat.events$family==p,"age.at.birth"])})
    pat.mean.by.parent[[i]]=cbind(pat.mean.by.parents,pat.mean.age.by.parents)

    mat.events$adjusted.nrec=mat.events$Freq-mat.mean.by.parents[mat.events$family]
    mat.events$adjusted.age=mat.events$age.at.birth-mat.mean.age.by.parents[mat.events$family]

    pat.events$adjusted.nrec=pat.events$Freq-pat.mean.by.parents[pat.events$family]
    pat.events$adjusted.age=pat.events$age.at.birth-pat.mean.age.by.parents[pat.events$family]
    
    data1.mat=rbind(data1.mat,mat.events)
    data1.pat=rbind(data1.pat,pat.events)
    
}
mat.cohort.count = table(data1.mat$cohort)
pat.cohort.count = table(data1.pat$cohort)

#only keep cohorts with more than 20
data1.mat=data1.mat[data1.mat$cohort %in% names(mat.cohort.count[mat.cohort.count>20]),]
data1.pat=data1.pat[data1.pat$cohort %in% names(pat.cohort.count[pat.cohort.count>20]),]

mat.key.model1=unique(data.frame(as.integer(as.factor(data1.mat$cohort)),as.factor(data1.mat$cohort)))
pat.key.model1=unique(data.frame(as.integer(as.factor(data1.pat$cohort)),as.factor(data1.pat$cohort)))
colnames(mat.key.model1)=c("Code","Cohort")
colnames(pat.key.model1)=c("Code","Cohort")
mat.key.model1=mat.key.model1[order(mat.key.model1$Code),]
pat.key.model1=pat.key.model1[order(pat.key.model1$Code),]
write.table(mat.key.model1,"RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.txt",quote=F,sep="\t",row.names=F)
write.table(pat.key.model1,"RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.txt",quote=F,sep="\t",row.names=F)


data1.mat2=list(y=data1.mat$Freq,cohort=as.integer(as.factor(data1.mat$cohort)),family=as.integer(as.factor(data1.mat$family)),Age=as.numeric(data1.mat$age.at.birth),J=nrow(data1.mat),I=length(unique(data1.mat$family)),
    C=length(unique(data1.mat$cohort)))
data1.pat2=list(y=data1.pat$Freq,cohort=as.integer(as.factor(data1.pat$cohort)),family=as.integer(as.factor(data1.pat$family)),Age=as.numeric(data1.pat$age.at.birth),J=nrow(data1.pat),I=length(unique(data1.pat$family)),
    C=length(unique(data1.pat$cohort)))

data1.adjusted.mat2=list(y=data1.mat$adjusted.nrec,cohort=as.integer(as.factor(data1.mat$cohort)),family=as.integer(as.factor(data1.mat$family)),Age=as.numeric(data1.mat$adjusted.age),
        J=nrow(data1.mat),I=length(unique(data1.mat$family)), C=length(unique(data1.mat$cohort)),nkids=data1.mat$nkids)
data1.adjusted.pat2=list(y=data1.pat$adjusted.nrec,cohort=as.integer(as.factor(data1.pat$cohort)),family=as.integer(as.factor(data1.pat$family)),Age=as.numeric(data1.pat$adjusted.age),
        J=nrow(data1.pat),I=length(unique(data1.pat$family)), C=length(unique(data1.pat$cohort)),nkids=data1.pat$nkids)

save.image("NFTOOLS_data_for_RSTAN.more_stringent.RData")
