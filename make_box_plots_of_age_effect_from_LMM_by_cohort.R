#"duoHMM_data_for_RSTAN.NTR_v2.RData"
library(nlme)
load("NFTOOLS_data_for_RSTAN.NTR_v2.RData")


### test effect of grandparental age in NTR 
ntr.age=read.delim("NTR/NTR_v2_minus_MZs.with_parental_and_grandparental_age.txt",header=T)

ntr.mat.nftools=data1.mat[data1.mat$cohort=="NTR_v2",]
ntr.mat.nftools$maternal_grandmum_age=ntr.age[match(as.character(ntr.mat.nftools$child),as.character(ntr.age$individual)),"maternal_grandmum_age"]
ntr.mat.nftools$maternal_granddad_age=ntr.age[match(as.character(ntr.mat.nftools$child),as.character(ntr.age$individual)),"maternal_granddad_age"]


 summary(lme(Freq~maternal_grandmum_age,~1|as.factor(family),data=ntr.mat.nftools[!is.na(ntr.mat.nftools$maternal_grandmum_age),]))
summary(lme(Freq~maternal_granddad_age,~1|as.factor(family),data=ntr.mat.nftools[!is.na(ntr.mat.nftools$maternal_granddad_age),]))





##### test parental age in NFTOOLS data, all cohorts


data1.mat.nftools=data1.mat
data1.pat.nftools=data1.pat
data1.mat.nftools$cohort=gsub("NTR_v2","NTR",data1.mat.nftools$cohort)
data1.mat.nftools$cohort=gsub("QTR_370K","QTR370",data1.mat.nftools$cohort)
data1.mat.nftools$cohort=gsub("QTR_610K","QTR610",data1.mat.nftools$cohort)

data1.pat.nftools$cohort=gsub("NTR_v2","NTR",data1.pat.nftools$cohort)
data1.pat.nftools$cohort=gsub("QTR_370K","QTR370",data1.pat.nftools$cohort)
data1.pat.nftools$cohort=gsub("QTR_610K","QTR610",data1.pat.nftools$cohort)


#all.nftools.mat.lme=lme(Freq~age.at.birth+nkids,~1|as.factor(family),data=data1.mat.nftools)
all.nftools.mat.lme=lme(Freq~age.at.birth,~1|as.factor(family),data=data1.mat.nftools)
summary(all.nftools.mat.lme)
#all.nftools.pat.lme=lme(Freq~age.at.birth+nkids,~1|as.factor(family),data=data1.pat.nftools)
all.nftools.pat.lme=lme(Freq~age.at.birth,~1|as.factor(family),data=data1.pat.nftools)
summary(all.nftools.pat.lme)

#simple linear regression on adjusted counts

summary(lm(adjusted.nrec~adjusted.age,data=data1.mat.nftools))
summary(lm(adjusted.nrec~adjusted.age,data=data1.pat.nftools))

###Linear Mixed Model with family as random effect
mat.cohorts=unique(data1.mat.nftools$cohort)
nftools.mat.results.by.cohort=list()

nftools.mat.adjusted.results.by.cohort=list()

nftools.mat.adjusted.spearman.by.cohort=NULL

for(c in 1:length(mat.cohorts)){
#    cohort.nftools.mat.lme=lme(Freq~age.at.birth+nkids,~1|as.factor(family),data=data1.mat.nftools[data1.mat.nftools$cohort==mat.cohorts[c],])
        cohort.nftools.mat.lme=lme(Freq~age.at.birth,~1|as.factor(family),data=data1.mat.nftools[data1.mat.nftools$cohort==mat.cohorts[c],])
    cohort.nftools.mat.lm=summary(lm(adjusted.nrec~adjusted.age,data=data1.mat.nftools[data1.mat.nftools$cohort==mat.cohorts[c],]))
    cor(data1.mat.nftools[data1.mat.nftools$cohort==mat.cohorts[c],"adjusted.nrec"],data1.mat.nftools[data1.mat.nftools$cohort==mat.cohorts[c],"adjusted.age"],method="spearman")
    nftools.mat.results.by.cohort[[c]]=summary(cohort.nftools.mat.lme)
    nftools.mat.adjusted.results.by.cohort[[c]]=cohort.nftools.mat.lm
    spearman.cor=cor.test(data1.mat.nftools[data1.mat.nftools$cohort==mat.cohorts[c],"adjusted.nrec"],data1.mat.nftools[data1.mat.nftools$cohort==mat.cohorts[c],"adjusted.age"],method="spearman")
    nftools.mat.adjusted.spearman.by.cohort=rbind(nftools.mat.adjusted.spearman.by.cohort,c(spearman.cor$estimate,spearman.cor$p.value))
    }
all.mat.spearman.cor = cor.test(data1.mat.nftools[,"adjusted.nrec"],data1.mat.nftools[,"adjusted.age"],method="spearman")
nftools.mat.adjusted.spearman.by.cohort=rbind(nftools.mat.adjusted.spearman.by.cohort,c(all.mat.spearman.cor$estimate,all.mat.spearman.cor$p.value))

colnames(nftools.mat.adjusted.spearman.by.cohort)=c("Spearman.corr","p.value")
rownames(nftools.mat.adjusted.spearman.by.cohort)=c(mat.cohorts,"all")


summary.nftools.mat.results.by.cohort= sapply(nftools.mat.results.by.cohort,function(x){return(x$tTable[2,])})
summary.nftools.mat.results.by.cohort=cbind(summary.nftools.mat.results.by.cohort,summary(all.nftools.mat.lme)$tTable[2,])
colnames(summary.nftools.mat.results.by.cohort)=c(mat.cohorts,"all")


summary.nftools.mat.adjusted.results.by.cohort= sapply(nftools.mat.adjusted.results.by.cohort,function(x){return(x$coefficients[2,])})
summary.nftools.mat.adjusted.results.by.cohort=cbind(summary.nftools.mat.adjusted.results.by.cohort,summary(lm(adjusted.nrec~adjusted.age,data=data1.mat.nftools))$coefficients[2,])
colnames(summary.nftools.mat.adjusted.results.by.cohort)=c(mat.cohorts,"all")

pat.cohorts=unique(data1.pat.nftools$cohort)
nftools.pat.results.by.cohort=list()
nftools.pat.adjusted.results.by.cohort=list()
nftools.pat.adjusted.spearman.by.cohort=NULL

for(c in 1:length(pat.cohorts)){
#    cohort.nftools.pat.lme=lme(Freq~age.at.birth+nkids,~1|as.factor(family),data=data1.pat.nftools[data1.pat.nftools$cohort==pat.cohorts[c],])
        cohort.nftools.pat.lme=lme(Freq~age.at.birth,~1|as.factor(family),data=data1.pat.nftools[data1.pat.nftools$cohort==pat.cohorts[c],])
    nftools.pat.results.by.cohort[[c]]=summary(cohort.nftools.pat.lme)
    cohort.nftools.pat.lm=summary(lm(adjusted.nrec~adjusted.age,data=data1.pat.nftools[data1.pat.nftools$cohort==pat.cohorts[c],]))
    nftools.pat.adjusted.results.by.cohort[[c]]=cohort.nftools.pat.lm
    
    spearman.cor=  cor.test(data1.pat.nftools[data1.pat.nftools$cohort==pat.cohorts[c],"adjusted.nrec"],data1.pat.nftools[data1.pat.nftools$cohort==pat.cohorts[c],"adjusted.age"],method="spearman")
    nftools.pat.adjusted.spearman.by.cohort=rbind(nftools.pat.adjusted.spearman.by.cohort,c(spearman.cor$estimate,spearman.cor$p.value))
}
all.pat.spearman.cor = cor.test(data1.pat.nftools[,"adjusted.nrec"],data1.pat.nftools[,"adjusted.age"],method="spearman")
nftools.pat.adjusted.spearman.by.cohort=rbind(nftools.pat.adjusted.spearman.by.cohort,c(all.pat.spearman.cor$estimate,all.pat.spearman.cor$p.value))
colnames(nftools.pat.adjusted.spearman.by.cohort)=c("Spearman.corr","p.value")
rownames(nftools.pat.adjusted.spearman.by.cohort)=c(pat.cohorts,"all")


summary.nftools.pat.results.by.cohort= sapply(nftools.pat.results.by.cohort,function(x){return(x$tTable[2,])})
summary.nftools.pat.results.by.cohort=cbind(summary.nftools.pat.results.by.cohort,summary(all.nftools.pat.lme)$tTable[2,])
colnames(summary.nftools.pat.results.by.cohort)=c(pat.cohorts,"all")

summary.nftools.pat.adjusted.results.by.cohort= sapply(nftools.pat.adjusted.results.by.cohort,function(x){return(x$coefficients[2,])})
summary.nftools.pat.adjusted.results.by.cohort=cbind(summary.nftools.pat.adjusted.results.by.cohort,summary(lm(adjusted.nrec~adjusted.age,data=data1.pat.nftools))$coefficients[2,])
colnames(summary.nftools.pat.adjusted.results.by.cohort)=c(pat.cohorts,"all")

load("duoHMM_data_for_RSTAN.NTR_v2.RData")
data1.mat.duohmm=data1.mat
data1.pat.duohmm=data1.pat

data2.mat.duohmm=data2.mat
data2.pat.duohmm=data2.pat

data1.mat.duohmm2=data1.mat.duohmm[as.character(data1.mat.duohmm$CHILD) %in% as.character(data1.mat.nftools$child),]
data1.pat.duohmm2=data1.pat.duohmm[as.character(data1.pat.duohmm$CHILD) %in% as.character(data1.pat.nftools$child),]

data1.adjusted.mat.duohmm=data1.adjusted.mat
data1.adjusted.pat.duohmm=data1.adjusted.pat

data1.adjusted.mat.duohmm2=data1.adjusted.mat.duohmm[as.character(data1.adjusted.mat.duohmm$CHILD) %in% as.character(data1.mat.nftools$child),]
data1.adjusted.pat.duohmm2=data1.adjusted.pat.duohmm[as.character(data1.adjusted.pat.duohmm$CHILD) %in% as.character(data1.pat.nftools$child),]

#all.duohmm.mat.lme=lme(nrec~age.at.birth+nkid,~1|as.factor(PARENT),data=data1.mat.duohmm)
all.duohmm.mat.lme=lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data1.mat.duohmm)
summary(all.duohmm.mat.lme)
#all.duohmm.pat.lme=lme(nrec~age.at.birth+nkid,~1|as.factor(PARENT),data=data1.pat.duohmm)
all.duohmm.pat.lme=lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data1.pat.duohmm)
summary(all.duohmm.pat.lme)

mat.cohorts=unique(data1.mat.duohmm$cohort)
duohmm.mat.results.by.cohort=list()
duohmm.mat.adjusted.results.by.cohort=list()
duohmm.mat.adjusted.spearman.by.cohort=NULL
for(c in 1:length(mat.cohorts)){
#    cohort.duohmm.mat.lme=lme(nrec~age.at.birth+nkid,~1|as.factor(PARENT),data=data1.mat.duohmm[data1.mat.duohmm$cohort==mat.cohorts[c],])
        cohort.duohmm.mat.lme=lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data1.mat.duohmm[data1.mat.duohmm$cohort==mat.cohorts[c],])
    duohmm.mat.results.by.cohort[[c]]=summary(cohort.duohmm.mat.lme)

    cohort.duohmm.mat.lm=summary(lm(adjusted.nrec~adjusted.age,data=data1.adjusted.mat.duohmm[data1.adjusted.mat.duohmm$cohort==mat.cohorts[c],]))
    duohmm.mat.adjusted.results.by.cohort[[c]]=cohort.duohmm.mat.lm

spearman.cor=cor.test(data1.adjusted.mat.duohmm[data1.adjusted.mat.duohmm$cohort==mat.cohorts[c],"adjusted.nrec"],data1.adjusted.mat.duohmm[data1.adjusted.mat.duohmm$cohort==mat.cohorts[c],"adjusted.age"],method="spearman")
        duohmm.mat.adjusted.spearman.by.cohort=rbind(duohmm.mat.adjusted.spearman.by.cohort,c(spearman.cor$estimate,spearman.cor$p.value))
    
}
all.mat.spearman.cor = cor.test(data1.adjusted.mat.duohmm[,"adjusted.nrec"],data1.adjusted.mat.duohmm[,"adjusted.age"],method="spearman")
duohmm.mat.adjusted.spearman.by.cohort=rbind(duohmm.mat.adjusted.spearman.by.cohort,c(all.mat.spearman.cor$estimate,all.mat.spearman.cor$p.value))
colnames(duohmm.mat.adjusted.spearman.by.cohort)=c("Spearman.corr","p.value")
rownames(duohmm.mat.adjusted.spearman.by.cohort)=c(mat.cohorts,"all")

summary.duohmm.mat.results.by.cohort= sapply(duohmm.mat.results.by.cohort,function(x){return(x$tTable[2,])})
summary.duohmm.mat.results.by.cohort=cbind(summary.duohmm.mat.results.by.cohort,summary(all.duohmm.mat.lme)$tTable[2,])
colnames(summary.duohmm.mat.results.by.cohort)=c(mat.cohorts,"all")

summary.duohmm.mat.adjusted.results.by.cohort= sapply(duohmm.mat.adjusted.results.by.cohort,function(x){return(x$coefficients[2,])})
summary.duohmm.mat.adjusted.results.by.cohort=cbind(summary.duohmm.mat.adjusted.results.by.cohort,summary(lm(adjusted.nrec~adjusted.age,data=data1.adjusted.mat.duohmm))$coefficients[2,])
colnames(summary.duohmm.mat.adjusted.results.by.cohort)=c(mat.cohorts,"all")





pat.cohorts=unique(data1.pat.duohmm$cohort)
duohmm.pat.results.by.cohort=list()
duohmm.pat.adjusted.results.by.cohort=list()
duohmm.pat.adjusted.spearman.by.cohort=NULL
for(c in 1:length(pat.cohorts)){
#    cohort.duohmm.pat.lme=lme(nrec~age.at.birth+nkid,~1|as.factor(PARENT),data=data1.pat.duohmm[data1.pat.duohmm$cohort==pat.cohorts[c],])
        cohort.duohmm.pat.lme=lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data1.pat.duohmm[data1.pat.duohmm$cohort==pat.cohorts[c],])
    duohmm.pat.results.by.cohort[[c]]=summary(cohort.duohmm.pat.lme)

    cohort.duohmm.pat.lm=summary(lm(adjusted.nrec~adjusted.age,data=data1.adjusted.pat.duohmm[data1.adjusted.pat.duohmm$cohort==pat.cohorts[c],]))
    duohmm.pat.adjusted.results.by.cohort[[c]]=cohort.duohmm.pat.lm
    
spearman.cor=cor.test(data1.adjusted.pat.duohmm[data1.adjusted.pat.duohmm$cohort==pat.cohorts[c],"adjusted.nrec"],data1.adjusted.pat.duohmm[data1.adjusted.pat.duohmm$cohort==pat.cohorts[c],"adjusted.age"],method="spearman")
    duohmm.pat.adjusted.spearman.by.cohort=rbind(duohmm.pat.adjusted.spearman.by.cohort,c(spearman.cor$estimate,spearman.cor$p.value))

}

all.pat.spearman.cor = cor.test(data1.adjusted.pat.duohmm[,"adjusted.nrec"],data1.adjusted.pat.duohmm[,"adjusted.age"],method="spearman")
duohmm.pat.adjusted.spearman.by.cohort=rbind(duohmm.pat.adjusted.spearman.by.cohort,c(all.pat.spearman.cor$estimate,all.pat.spearman.cor$p.value))
colnames(duohmm.pat.adjusted.spearman.by.cohort)=c("Spearman.corr","p.value")
rownames(duohmm.pat.adjusted.spearman.by.cohort)=c(pat.cohorts,"all")

summary.duohmm.pat.results.by.cohort= sapply(duohmm.pat.results.by.cohort,function(x){return(x$tTable[2,])})
summary.duohmm.pat.results.by.cohort=cbind(summary.duohmm.pat.results.by.cohort,summary(all.duohmm.pat.lme)$tTable[2,])
colnames(summary.duohmm.pat.results.by.cohort)=c(pat.cohorts,"all")

summary.duohmm.pat.adjusted.results.by.cohort= sapply(duohmm.pat.adjusted.results.by.cohort,function(x){return(x$coefficients[2,])})
summary.duohmm.pat.adjusted.results.by.cohort=cbind(summary.duohmm.pat.adjusted.results.by.cohort,summary(lm(adjusted.nrec~adjusted.age,data=data1.adjusted.pat.duohmm))$coefficients[2,])
colnames(summary.duohmm.pat.adjusted.results.by.cohort)=c(pat.cohorts,"all")
####try including uninromative duos too

mat.cohorts=unique(data2.mat.duohmm$cohort.family.type)
uninform.duohmm.mat.results.by.cohort=list()

for(c in 1:length(mat.cohorts)){
#    cohort.duohmm.mat.lme=lme(nrec~age.at.birth+nkid,~1|as.factor(PARENT),data=data2.mat.duohmm[data2.mat.duohmm$cohort.family.type==mat.cohorts[c],])
        cohort.duohmm.mat.lme=lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data2.mat.duohmm[data2.mat.duohmm$cohort.family.type==mat.cohorts[c],])
    uninform.duohmm.mat.results.by.cohort[[c]]=summary(cohort.duohmm.mat.lme)
    
}
summary.uninform.duohmm.mat.results.by.cohort= sapply(uninform.duohmm.mat.results.by.cohort,function(x){return(x$tTable[2,])})
colnames(summary.uninform.duohmm.mat.results.by.cohort)=c(mat.cohorts)

###just duoHMM data for families in NFTOOLS results (phase informed)
#all.duohmm2.mat.lme=lme(nrec~age.at.birth+nkid,~1|as.factor(PARENT),data=data1.mat.duohmm2)
all.duohmm2.mat.lme=lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data1.mat.duohmm2)
summary(all.duohmm2.mat.lme)
#all.duohmm2.pat.lme=lme(nrec~age.at.birth+nkid,~1|as.factor(PARENT),data=data1.pat.duohmm2)
all.duohmm2.pat.lme=lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data1.pat.duohmm2)
summary(all.duohmm2.pat.lme)


mat.cohorts=unique(data1.mat.duohmm2$cohort)
duohmm2.mat.results.by.cohort=list()
duohmm2.mat.adjusted.results.by.cohort=list()
duohmm2.mat.adjusted.spearman.by.cohort=NULL

for(c in 1:length(mat.cohorts)){
#    cohort.duohmm2.mat.lme=lme(nrec~age.at.birth+nkid,~1|as.factor(PARENT),data=data1.mat.duohmm2[data1.mat.duohmm2$cohort==mat.cohorts[c],])
        cohort.duohmm2.mat.lme=lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data1.mat.duohmm2[data1.mat.duohmm2$cohort==mat.cohorts[c],])
    duohmm2.mat.results.by.cohort[[c]]=summary(cohort.duohmm2.mat.lme)

    cohort.duohmm.mat.lm=summary(lm(adjusted.nrec~adjusted.age,data=data1.adjusted.mat.duohmm2[data1.adjusted.mat.duohmm2$cohort==mat.cohorts[c],]))
    duohmm2.mat.adjusted.results.by.cohort[[c]]=cohort.duohmm.mat.lm
    
    spearman.cor=cor.test(data1.adjusted.mat.duohmm2[data1.adjusted.mat.duohmm2$cohort==mat.cohorts[c],"adjusted.nrec"],data1.adjusted.mat.duohmm2[data1.adjusted.mat.duohmm2$cohort==mat.cohorts[c],"adjusted.age"],
            method="spearman")
        duohmm2.mat.adjusted.spearman.by.cohort=rbind(duohmm2.mat.adjusted.spearman.by.cohort,c(spearman.cor$estimate,spearman.cor$p.value))
    
}
all.mat.spearman.cor = cor.test(data1.adjusted.mat.duohmm2[,"adjusted.nrec"],data1.adjusted.mat.duohmm2[,"adjusted.age"],method="spearman")
duohmm2.mat.adjusted.spearman.by.cohort=rbind(duohmm2.mat.adjusted.spearman.by.cohort,c(all.mat.spearman.cor$estimate,all.mat.spearman.cor$p.value))
colnames(duohmm2.mat.adjusted.spearman.by.cohort)=c("Spearman.corr","p.value")
rownames(duohmm2.mat.adjusted.spearman.by.cohort)=c(mat.cohorts,"all")

summary.duohmm2.mat.adjusted.results.by.cohort= sapply(duohmm2.mat.adjusted.results.by.cohort,function(x){return(x$coefficients[2,])})
summary.duohmm2.mat.adjusted.results.by.cohort=cbind(summary.duohmm2.mat.adjusted.results.by.cohort,summary(lm(adjusted.nrec~adjusted.age,data=data1.adjusted.mat.duohmm2))$coefficients[2,])
colnames(summary.duohmm2.mat.adjusted.results.by.cohort)=c(mat.cohorts,"all")


summary.duohmm2.mat.results.by.cohort= sapply(duohmm2.mat.results.by.cohort,function(x){return(x$tTable[2,])})
summary.duohmm2.mat.results.by.cohort=cbind(summary.duohmm2.mat.results.by.cohort,summary(all.duohmm2.mat.lme)$tTable[2,])
colnames(summary.duohmm2.mat.results.by.cohort)=c(mat.cohorts,"all")

pat.cohorts=unique(data1.pat.duohmm2$cohort)
duohmm2.pat.results.by.cohort=list()
duohmm2.pat.adjusted.results.by.cohort=list()
duohmm2.pat.adjusted.spearman.by.cohort=NULL

for(c in 1:length(pat.cohorts)){
#    cohort.duohmm2.pat.lme=lme(nrec~age.at.birth+nkid,~1|as.factor(PARENT),data=data1.pat.duohmm2[data1.pat.duohmm2$cohort==pat.cohorts[c],])
        cohort.duohmm2.pat.lme=lme(nrec~age.at.birth,~1|as.factor(PARENT),data=data1.pat.duohmm2[data1.pat.duohmm2$cohort==pat.cohorts[c],])
    duohmm2.pat.results.by.cohort[[c]]=summary(cohort.duohmm2.pat.lme)

    cohort.duohmm.pat.lm=summary(lm(adjusted.nrec~adjusted.age,data=data1.adjusted.pat.duohmm2[data1.adjusted.pat.duohmm2$cohort==pat.cohorts[c],]))
    duohmm2.pat.adjusted.results.by.cohort[[c]]=cohort.duohmm.pat.lm
    
spearman.cor=cor.test(data1.adjusted.pat.duohmm2[data1.adjusted.pat.duohmm2$cohort==pat.cohorts[c],"adjusted.nrec"],data1.adjusted.pat.duohmm2[data1.adjusted.pat.duohmm2$cohort==pat.cohorts[c],"adjusted.age"],
    method="spearman")
    duohmm2.pat.adjusted.spearman.by.cohort=rbind(duohmm2.pat.adjusted.spearman.by.cohort,c(spearman.cor$estimate,spearman.cor$p.value))


}

all.pat.spearman.cor = cor.test(data1.adjusted.pat.duohmm2[,"adjusted.nrec"],data1.adjusted.pat.duohmm2[,"adjusted.age"],method="spearman")
duohmm2.pat.adjusted.spearman.by.cohort=rbind(duohmm2.pat.adjusted.spearman.by.cohort,c(all.pat.spearman.cor$estimate,all.pat.spearman.cor$p.value))
colnames(duohmm2.pat.adjusted.spearman.by.cohort)=c("Spearman.corr","p.value")
rownames(duohmm2.pat.adjusted.spearman.by.cohort)=c(pat.cohorts,"all")

summary.duohmm2.pat.adjusted.results.by.cohort= sapply(duohmm2.pat.adjusted.results.by.cohort,function(x){return(x$coefficients[2,])})
summary.duohmm2.pat.adjusted.results.by.cohort=cbind(summary.duohmm2.pat.adjusted.results.by.cohort,summary(lm(adjusted.nrec~adjusted.age,data=data1.adjusted.pat.duohmm2))$coefficients[2,])
colnames(summary.duohmm2.pat.adjusted.results.by.cohort)=c(pat.cohorts,"all")

summary.duohmm2.pat.results.by.cohort= sapply(duohmm2.pat.results.by.cohort,function(x){return(x$tTable[2,])})
summary.duohmm2.pat.results.by.cohort=cbind(summary.duohmm2.pat.results.by.cohort,summary(all.duohmm2.pat.lme)$tTable[2,])
colnames(summary.duohmm2.pat.results.by.cohort)=c(pat.cohorts,"all")


    colnames(nftools.mat.adjusted.spearman.by.cohort) = paste0("nftools.",colnames(nftools.mat.adjusted.spearman.by.cohort))
    colnames(duohmm.mat.adjusted.spearman.by.cohort)=paste0("duohmm.",colnames(duohmm.mat.adjusted.spearman.by.cohort))
    colnames(duohmm2.mat.adjusted.spearman.by.cohort)=paste0("duohmm.in.nftools.",colnames(duohmm2.mat.adjusted.spearman.by.cohort))
    
    colnames(summary.nftools.mat.results.by.cohort) = paste0("nftools.",colnames(summary.nftools.mat.results.by.cohort))
    colnames(summary.duohmm.mat.results.by.cohort)=paste0("duohmm.",colnames(summary.duohmm.mat.results.by.cohort))
    colnames(summary.duohmm2.mat.results.by.cohort)=paste0("duohmm.in.nftools.",colnames(summary.duohmm2.mat.results.by.cohort))
    
    colnames(summary.nftools.mat.adjusted.results.by.cohort)= paste0("nftools.",colnames(summary.nftools.mat.adjusted.results.by.cohort))
    colnames(summary.duohmm.mat.adjusted.results.by.cohort)=paste0("duohmm.",colnames(summary.duohmm.mat.adjusted.results.by.cohort))
    colnames(summary.duohmm2.mat.adjusted.results.by.cohort)=paste0("duohmm.in.nftools.",colnames(summary.duohmm2.mat.adjusted.results.by.cohort))


write.table(summary.nftools.mat.results.by.cohort,"frequentist_analysis/nftools.mat.lme.results.by.cohort.no_nkids_term.txt",quote=F,sep="\t")
write.table(summary.duohmm.mat.results.by.cohort,"frequentist_analysis/duohmm.mat.lme.results.by.cohort.no_nkids_term.txt",quote=F,sep="\t")
write.table(summary.duohmm2.mat.results.by.cohort,"frequentist_analysis/duohmm_in_nftools.mat.lme.results.by.cohort.no_nkids_term.txt",quote=F,sep="\t")

write.table(summary.nftools.pat.results.by.cohort,"frequentist_analysis/nftools.pat.lme.results.by.cohort.no_nkids_term.txt",quote=F,sep="\t")
write.table(summary.duohmm.pat.results.by.cohort,"frequentist_analysis/duohmm.pat.lme.results.by.cohort.no_nkids_term.txt",quote=F,sep="\t")
write.table(summary.duohmm2.pat.results.by.cohort,"frequentist_analysis/duohmm_in_nftools.pat.lme.results.by.cohort.no_nkids_term.txt",quote=F,sep="\t")


mycols=c("black","blue","red","green","orange","red4","purple","darkgreen")
names(mycols)=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB")

library(metafor)
#fixed effects meta-analysis of betas
meta.fe=rma.uni(yi=summary.duohmm2.mat.results.by.cohort[1,1:6],sei=summary.duohmm2.mat.results.by.cohort[2,1:6],method="FE",weighted=T)

pdf("/home/hilary/maternal_age_recombination/frequentist_analysis/box_plot_of_beta_Age_from_LMM_on_informative_nuclear_families_duoHMM.no_nkids_term.pdf",height=4,width=5)
plot(1:7,c(summary.duohmm2.mat.results.by.cohort[1,1:6],meta.fe$b),main="beta_age from LMM",ylim=c(-0.7,1),col=c(rep("white",6),"black"),pch=21,xaxt="n",xlab="Cohort",ylab="beta_age",cex=0.5)
axis(1,at=1:7,labels=c(gsub("duohmm.in.nftools.","",colnames(summary.duohmm2.mat.results.by.cohort)[1:6],fixed=T),"meta-analysis"),tick=F,cex.axis=0.5)
abline(h=0)
sample.sizes= table(data1.mat.duohmm2$cohort)
for(i in 1:6){
my.x=i
my.y=summary.duohmm2.mat.results.by.cohort[1,i]
#scaling=0.05*sample.sizes[i]/sample.sizes[2]
scaling=0.05
#polygon(c(my.x-scaling,my.x-scaling,my.x+scaling,my.x+scaling),c(my.y-scaling,my.y+scaling,my.y+scaling,my.y-scaling))
offset=0.02*sample.sizes[i]/sample.sizes[2]
segments(x0=i,y0=summary.duohmm2.mat.results.by.cohort[1,i]-1.96*summary.duohmm2.mat.results.by.cohort[2,i],x1=i,y1=summary.duohmm2.mat.results.by.cohort[1,i]+1.96*summary.duohmm2.mat.results.by.cohort[2,i])
polygon(c(my.x-scaling,my.x-scaling,my.x+scaling,my.x+scaling),c(my.y-offset,my.y+offset,my.y+offset,my.y-offset),col=mycols[names(sample.sizes)[i]],border=mycols[names(sample.sizes)[i]])
}
segments(x0=7,y0=meta.fe$ci.lb,x1=7,y1=meta.fe$ci.ub)

dev.off()



####for all informative duos
meta.fe.all=rma.uni(yi=summary.duohmm.mat.results.by.cohort[1,1:8],sei=summary.duohmm.mat.results.by.cohort[2,1:8],method="FE",weighted=T)

pdf("/home/hilary/maternal_age_recombination/frequentist_analysis/box_plot_of_beta_Age_from_LMM_on_all_informative_duos_duoHMM.no_nkids_term.pdf",height=4,width=5)
plot(1:9,c(summary.duohmm.mat.results.by.cohort[1,1:8],meta.fe.all$b),main="beta_age from LMM",ylim=c(-0.7,1),col=c(rep("white",8),"black"),pch=21,xaxt="n",xlab="Cohort",ylab="beta_age",cex=0.5)
axis(1,at=1:9,labels=c(gsub("duohmm.","",colnames(summary.duohmm.mat.results.by.cohort)[1:8],fixed=T),"meta-analysis"),tick=F,cex.axis=0.5)
abline(h=0)
sample.sizes= table(data1.mat.duohmm$cohort)
for(i in 1:8){
my.x=i
my.y=summary.duohmm.mat.results.by.cohort[1,i]
#scaling=0.05*sample.sizes[i]/sample.sizes[2]
scaling=0.05
#polygon(c(my.x-scaling,my.x-scaling,my.x+scaling,my.x+scaling),c(my.y-scaling,my.y+scaling,my.y+scaling,my.y-scaling))
offset=0.02*sample.sizes[i]/sample.sizes[2]
segments(x0=i,y0=summary.duohmm.mat.results.by.cohort[1,i]-1.96*summary.duohmm.mat.results.by.cohort[2,i],x1=i,y1=summary.duohmm.mat.results.by.cohort[1,i]+1.96*summary.duohmm.mat.results.by.cohort[2,i])
polygon(c(my.x-scaling,my.x-scaling,my.x+scaling,my.x+scaling),c(my.y-offset,my.y+offset,my.y+offset,my.y-offset),col=mycols[names(sample.sizes)[i]],border=mycols[names(sample.sizes)[i]])
}
segments(x0=9,y0=meta.fe.all$ci.lb,x1=9,y1=meta.fe.all$ci.ub)

dev.off()

####for all duos from families with >2 kids

meta.fe.uninf=rma.uni(yi=summary.uninform.duohmm.mat.results.by.cohort[1,1:24],sei=summary.uninform.duohmm.mat.results.by.cohort[2,1:24],method="FE",weighted=T)

sample.sizes= table(data2.mat.duohmm$cohort.family.type)
key=unique(data2.mat[,c("cohort","cohort.family.type")])
colnames(key)=c("cohort","cohort.family.type")
key$family.type=NA
for(j in 1:nrow(key)){
key$flamily.type[j]=gsub(paste0(key[j,1],"."),"",key[j,2])
}
rownames(key)=key$cohort.family.type
family.types=unique(key$family.type)

my.lty=c(1,2,3,1,2,3)
names(my.lty)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent" )
my.lwd=c(3,1,3,1,3,1)
names(my.lwd)=names(my.lty)

pdf("/home/hilary/maternal_age_recombination/frequentist_analysis/box_plot_of_beta_Age_from_LMM_on_all_duos_from_fams_with_multiple_kids.duoHMM.no_nkids_term.pdf",height=4,width=6)
plot(1:25,c(summary.uninform.duohmm.mat.results.by.cohort[1,1:24],meta.fe.uninf$b),main="beta_age from LMM",ylim=c(-1.1,1.4),col=c(rep("white",24),"black"),pch=19,xaxt="n",
     xlab="",ylab="beta_age",cex=0.5)
abline(h=0)
axis(1,at=1:25,labels=c(colnames(summary.uninform.duohmm.mat.results.by.cohort)[1:24],"meta-analysis"),tick=F,cex.axis=0.3,las=3,padj=0)
for(i in 1:24){
my.x=i
my.y=summary.uninform.duohmm.mat.results.by.cohort[1,i]

this.lwd=my.lwd[key[names(sample.sizes)[i],"family.type"]]
this.lty=my.lty[key[names(sample.sizes)[i],"family.type"]]
offset=0.01*sample.sizes[i]/sample.sizes[2]
scaling=0.05
se = summary.uninform.duohmm.mat.results.by.cohort[2,i]
segments(x0=i,y0=summary.uninform.duohmm.mat.results.by.cohort[1,i]-1.96*se,x1=i,y1=summary.uninform.duohmm.mat.results.by.cohort[1,i]+1.96*se)
polygon(c(my.x-scaling,my.x-scaling,my.x+scaling,my.x+scaling),c(my.y-offset,my.y+offset,my.y+offset,my.y-offset),col=mycols[key[names(sample.sizes)[i],"cohort"]],border=mycols[key[names(sample.sizes)[i],"cohort"]])
                                        #mycols[key[names(sample.sizes),"cohort"]]
}
segments(x0=25,y0=meta.fe.uninf$ci.lb,x1=25,y1=meta.fe.uninf$ci.ub)
dev.off()


