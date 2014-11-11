load("NFTOOLS_results/NFTOOLS_counts.all_families_except_outliers.RData")

age.files=c("NTR/NTR_minus_MZs.with_parental_age.txt","NTR/NTR_minus_MZs.with_parental_age.txt","QTR/QTR_minus_MZs.370K.with_parental_age.txt","QTR/QTR_minus_MZs.610K.with_parental_age.txt","FC/FC.with_parental_age.txt","GPC/GPC.with_parental_age.txt",
    "VB/valborbera_b37-related.with_parental_age.txt")
names(age.files)=names(all.events)

mat.inf=NULL
pat.inf=NULL
mat.uninf=NULL
pat.uninf=NULL
#for(i in 1:length(all.events)){
for(i in c(1,3,4,5,6,7)){
    
    events=all.events[[i]]
    ages=read.delim(age.files[i],header=T,stringsAsFactors=F)
    ages$individual=gsub("-","_",ages$individual)
    ages$mother=gsub("-","_",ages$mother)
    ages$father=gsub("-","_",ages$father)
    mat.events.inf=events[events$sex=="female" & events$nkids>2,]
    pat.events.inf=events[events$sex=="male" & events$nkids>2,]

    mat.events.inf$age=ages[match(mat.events.inf$child,ages$individual),"maternal_age"]
    pat.events.inf$age=ages[match(pat.events.inf$child,ages$individual),"paternal_age"]
    
    mat.events.uninf=events[events$sex=="female" & events$nkids==2,]
    pat.events.uninf=events[events$sex=="male" & events$nkids==2,]

    mat.events.uninf$sib1.age=ages[match(mat.events.uninf$child,ages$individual),"maternal_age"]
    pat.events.uninf$sib1.age=ages[match(pat.events.uninf$child,ages$individual),"paternal_age"]

    mat.events.uninf$sib2=NA
    pat.events.uninf$sib2=NA

    mat.events.uninf$sib2.age=NA
    pat.events.uninf$sib2.age=NA

    
    for(j in 1:nrow(mat.events.uninf)){
        mat=ages[ages$individual==mat.events.uninf[j,"child"],"mother"]
        sibs=ages[ages$mother == mat,"individual"]
        sib=sibs[1]
        if(sibs[1]==mat.events.uninf[j,"child"]){
            sib=sibs[2]
        }
        mat.events.uninf$sib2[j]=sib
        mat.events.uninf$sib2.age[j]=ages[ages$individual==sib,"maternal_age"]        
    }
    
    for(j in 1:nrow(pat.events.uninf)){
        pat=ages[ages$individual==pat.events.uninf[j,"child"],"father"]
        sibs=ages[ages$father == pat,"individual"]
        sib=sibs[1]
        if(sibs[1]==pat.events.uninf[j,"child"]){
            sib=sibs[2]
        }
        pat.events.uninf$sib2[j]=sib
        pat.events.uninf$sib2.age[j]=ages[ages$individual==sib,"paternal_age"]        
    }
   
    mat.events.uninf$av.age=apply(mat.events.uninf[,c("sib1.age","sib2.age")],1,mean)
    pat.events.uninf$av.age=apply(pat.events.uninf[,c("sib1.age","sib2.age")],1,mean)
    mat.events.uninf$av.Freq=mat.events.uninf$Freq/2
    pat.events.uninf$av.Freq=pat.events.uninf$Freq/2

    mat.inf=rbind(mat.inf,mat.events.inf)
    pat.inf=rbind(pat.inf,pat.events.inf)

    mat.uninf=rbind(mat.uninf,mat.events.uninf)
    pat.uninf=rbind(pat.uninf,pat.events.uninf)

}
mat.inf2=mat.inf[,c("child","Freq","age")]
mat.uninf2=mat.uninf[,c("child","av.Freq","av.age")]
colnames(mat.uninf2)=c("child","nrec","age.at.birth")
colnames(mat.inf2)=c("child","nrec","age.at.birth")
pat.inf2=pat.inf[,c("child","Freq","age")]
pat.uninf2=pat.uninf[,c("child","av.Freq","av.age")]
colnames(pat.uninf2)=c("child","nrec","age.at.birth")
colnames(pat.inf2)=c("child","nrec","age.at.birth")
data1.2.mat=rbind(mat.inf2,mat.uninf2)
data1.2.pat=rbind(pat.inf2,pat.uninf2)


data1.2.mat$binned.age = 20
data1.2.mat$binned.age[data1.2.mat$age.at.birth >20 & data1.2.mat$age.at.birth<=25] = 25
data1.2.mat$binned.age[data1.2.mat$age.at.birth >25 & data1.2.mat$age.at.birth<=30] = 30
data1.2.mat$binned.age[data1.2.mat$age.at.birth >30 & data1.2.mat$age.at.birth<=35] = 35
data1.2.mat$binned.age[data1.2.mat$age.at.birth >35 & data1.2.mat$age.at.birth<=40] = 40
data1.2.mat$binned.age[data1.2.mat$age.at.birth >40] = 45

data1.2.pat$binned.age = 20
data1.2.pat$binned.age[data1.2.pat$age.at.birth >20 & data1.2.pat$age.at.birth<=25] = 25
data1.2.pat$binned.age[data1.2.pat$age.at.birth >25 & data1.2.pat$age.at.birth<=30] = 30
data1.2.pat$binned.age[data1.2.pat$age.at.birth >30 & data1.2.pat$age.at.birth<=35] = 35
data1.2.pat$binned.age[data1.2.pat$age.at.birth >35 & data1.2.pat$age.at.birth<=40] = 40
data1.2.pat$binned.age[data1.2.pat$age.at.birth >40] = 45

data1.2.mat$nrec.normalised=data1.2.mat$nrec-mean(data1.2.mat$nrec[data1.2.mat$binned.age==25])
data1.2.pat$nrec.normalised=data1.2.pat$nrec-mean(data1.2.pat$nrec[data1.2.pat$binned.age==25])

data1.2.mat.by.bin=split(data1.2.mat,data1.2.mat$binned.age)
data1.2.pat.by.bin=split(data1.2.pat,data1.2.pat$binned.age)

mat.binned.rec = as.data.frame(do.call("rbind",lapply(data1.2.mat.by.bin,function(rec){
    mean.rec=    mean(rec$nrec.normalised)
    lower.rec=mean.rec-qnorm(0.975)*sd(rec$nrec.normalised)/sqrt(length(rec$nrec))
    upper.rec=mean.rec+qnorm(0.975)*sd(rec$nrec.normalised)/sqrt(length(rec$nrec))
    return(c(mean.rec,lower.rec,upper.rec))
})))
colnames(mat.binned.rec)=c("mean","lower.95CI","upper.95CI")
pat.binned.rec = as.data.frame(do.call("rbind",lapply(data1.2.pat.by.bin,function(rec){
    mean.rec=    mean(rec$nrec.normalised)
    lower.rec=mean.rec-qnorm(0.975)*sd(rec$nrec.normalised)/sqrt(length(rec$nrec))
    upper.rec=mean.rec+qnorm(0.975)*sd(rec$nrec.normalised)/sqrt(length(rec$nrec))
    return(c(mean.rec,lower.rec,upper.rec))
})))
colnames(pat.binned.rec)=c("mean","lower.95CI","upper.95CI")
mat.binned.rec$x=c(17.5,22.5,27.5,32.5,37.5,42.5)
pat.binned.rec$x=c(17.5,22.5,27.5,32.5,37.5,42.5)+0.5
pdf("/home/hilary/maternal_age_recombination/crossover_count_by_binned_age.NFTOOLS.informative_only.pdf",height=5,width=5)
plot(mat.binned.rec$x,mat.binned.rec$mean,xlab="Parental age",ylab="Difference from 20-25 age group",col="red",ylim=c(-4,4))
points(pat.binned.rec$x,pat.binned.rec$mean,col="blue")
apply(mat.binned.rec,1,function(bounds){segments(x0=bounds[4],y0=bounds[2],x1=bounds[4],y1=bounds[3],col="red")})
apply(pat.binned.rec,1,function(bounds){segments(x0=bounds[4],y0=bounds[2],x1=bounds[4],y1=bounds[3],col="blue")})
abline(h=0)
dev.off()
