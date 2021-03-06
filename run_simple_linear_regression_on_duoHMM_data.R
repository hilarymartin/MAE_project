
#load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.RData")
load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.NTR_v2.RData")
#for informative meioses only with nkid=>2, all cohorts, fit cohort effect and family effect
#cohort-specific age effect
#model1.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.stan",data=data1.mat2,iter=10000,chains=4)
#model1.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.stan",data=data1.pat2,iter=10000,chains=4)

summary(lm(data1.mat2$y~data1.mat2$Age))

summary(lm(data1.mat$nrec~data1.mat$age.at.birth))
summary(lm(data1.mat$nrec~data1.mat$age.at.birth+data1.mat$cohort+data1.mat$age.at.birth*data1.mat$cohort))
summary(aov(data1.mat$nrec~data1.mat$age.at.birth+as.factor(data1.mat$cohort)+data1.mat$age.at.birth*as.factor(data1.mat$cohort)))

library(nlme)

data1.mat$MOTHER=as.factor(data1.mat$PARENT)
lme(nrec~age.at.birth+cohort,random=~1|MOTHER,data=data1.mat)

summary(lm(data1.pat$nrec~data1.pat$age.at.birth))
summary(lm(data1.pat$nrec~data1.pat$age.at.birth+data1.pat$cohort+data1.pat$age.at.birth*data1.pat$cohort))
summary(aov(data1.pat$nrec~data1.pat$age.at.birth+as.factor(data1.pat$cohort)+data1.pat$age.at.birth*as.factor(data1.pat$cohort)))


#for all informative meioses
summary(lm(data1.2.mat$nrec~data1.2.mat$age.at.birth))
summary(lm(data1.2.mat$nrec~data1.2.mat$age.at.birth+data1.2.mat$cohort))
summary(lm(data1.2.mat$nrec~data1.2.mat$age.at.birth+data1.2.mat$cohort+data1.2.mat$age.at.birth*data1.2.mat$cohort))
summary(aov(data1.2.mat$nrec~data1.2.mat$age.at.birth+as.factor(data1.2.mat$cohort)+data1.2.mat$age.at.birth*as.factor(data1.2.mat$cohort)))
summary(aov(data1.2.mat$nrec~data1.2.mat$age.at.birth+as.factor(data1.2.mat$cohort)+data1.2.mat$age.at.birth*as.factor(data1.2.mat$cohort)+Error(data1.2.mat$PARENT)))

summary(lm(data1.2.pat$nrec~data1.2.pat$age.at.birth))
summary(lm(data1.2.pat$nrec~data1.2.pat$age.at.birth+data1.2.pat$cohort))
summary(lm(data1.2.pat$nrec~data1.2.pat$age.at.birth+data1.2.pat$cohort+data1.2.pat$age.at.birth*data1.2.pat$cohort))
summary(aov(data1.2.pat$nrec~data1.2.pat$age.at.birth+as.factor(data1.2.pat$cohort)+data1.2.pat$age.at.birth*as.factor(data1.2.pat$cohort)))


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
#pdf("/home/hilary/maternal_age_recombination/crossover_count_by_binned_age.informative_only.pdf",height=5,width=5)

pdf("/home/hilary/maternal_age_recombination/crossover_count_by_binned_age.informative_only.with_NTR_v2.pdf",height=5,width=5)
plot(mat.binned.rec$x,mat.binned.rec$mean,xlab="Parental age",ylab="Difference from 20-25 age group",col="red",ylim=c(-4,4))
points(pat.binned.rec$x,pat.binned.rec$mean,col="blue")
apply(mat.binned.rec,1,function(bounds){segments(x0=bounds[4],y0=bounds[2],x1=bounds[4],y1=bounds[3],col="red")})
apply(pat.binned.rec,1,function(bounds){segments(x0=bounds[4],y0=bounds[2],x1=bounds[4],y1=bounds[3],col="blue")})
abline(h=0)
dev.off()

summary(aov(data1.2.mat$nrec~data1.2.mat$binned.age)) #significant
summary(aov(data1.2.mat$nrec~data1.2.mat$binned.age+as.factor(data1.2.mat$cohort)))  #significant       
summary(aov(data1.2.mat$nrec~as.factor(data1.2.mat$binned.age))) #not significant
summary(aov(data1.2.mat$nrec~as.factor(data1.2.mat$binned.age)+as.factor(data1.2.mat$cohort)))  #not significant    

summary(aov(data1.2.mat$nrec~data1.2.mat$age.at.birth)) #significant
summary(aov(data1.2.mat$nrec~data1.2.mat$age.at.birth+as.factor(data1.2.mat$cohort)))  #significant

summary(lm(data1.2.mat$nrec~data1.2.mat$age.at.birth)) #significant
summary(lm(data1.2.mat$nrec~data1.2.mat$age.at.birth+as.factor(data1.2.mat$cohort))) # not significant


####is the age distribution different for different cohorts

data1.2.mat.by.cohort=split(data1.2.mat,data1.2.mat$cohort)
data1.2.pat.by.cohort=split(data1.2.pat,data1.2.pat$cohort)

cols = c("black","blue","darkgreen","red","pink2","purple","orange","orange3","red4","pink3")
mylty=c(1,2,4,1,2,4,1,2,4,1)
pdf("/home/hilary/maternal_age_recombination/parental_age_distribution_by_cohort.informative_only.pdf",height=5,width=10)
par(mfrow=c(1,2))
for(i in 1:length(data1.2.mat.by.cohort)){
    if(i==1){
        plot(density(data1.2.mat.by.cohort[[i]]$age.at.birth),main="Maternal age",xlab="Maternal age",col=cols[i],lty=mylty[i],xlim=c(12,50),ylim=c(0,0.1),lwd=2)
    } else {
        lines(density(data1.2.mat.by.cohort[[i]]$age.at.birth),col=cols[i],lty=mylty[i],lwd=2)
    }
}
legend("topright",paste(names(data1.2.mat.by.cohort),", n = ",unlist(lapply(data1.2.mat.by.cohort,nrow)),sep=" "),col=cols,lty=mylty,cex=0.6,lwd=2)
for(i in 1:length(data1.2.pat.by.cohort)){
    if(i==1){
        plot(density(data1.2.pat.by.cohort[[i]]$age.at.birth),main="Paternal age",xlab="Paternal age",col=cols[i],lty=mylty[i],xlim=c(15,60),ylim=c(0,0.1),lwd=2)
    } else {
        lines(density(data1.2.pat.by.cohort[[i]]$age.at.birth),col=cols[i],lty=mylty[i],lwd=2)
    }
}
legend("topright",paste(names(data1.2.pat.by.cohort),", n = ",unlist(lapply(data1.2.pat.by.cohort,nrow)),sep=" "),col=cols,lty=mylty,cex=0.6,lwd=2)
dev.off()
    


### Model 2.
#for informative meioses only, all cohorts, fit cohort effect, BUT NO FAMILY EFFECT
#cohort-specific age effect
#common age effect



#### Model 3.
#for nkid=>2, family effect, and fit N_called_crossovers~Binom(Y,p=f(cohort*family_type)); only for those with 
#data2.mat2=list(y=data2.mat$nrec,cohort=as.factor(data2.mat$cohort),family=as.factor(data2.mat$PARENT),age=as.numeric(data2.mat$age.at.birth),J=nrow(data2.mat),I=length(unique(data2.mat$PARENT)))
#data2.pat2=list(y=data2.pat$nrec,cohort=as.factor(data2.pat$cohort),family=as.factor(data2.pat$PARENT),age=as.numeric(data2.pat$age.at.birth),J=nrow(data2.pat),I=length(unique(data2.pat$PARENT)))

#model3.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.stan",data=data2.mat2,iter=100,chains=4)
#model3.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.stan",data=data2.pat2,iter=10000,chains=4)


### Model 4. N.B. NOT SURE SHOULD INCLUDE FAMILIES WITH ONLY 1 KID - DOESN'T REALLY LOOK BINOMIAL
#for all meioses, NO FAMILY EFFECT; fit N_called_crossovers~Binom(Y,p=f(cohort*family_type)); only for those with
#data3.mat2=list(y=data3.mat$nrec,cohort=as.factor(data3.mat$cohort),age=as.numeric(data3.mat$age.at.birth),J=nrow(data3.mat),I=length(unique(data3.mat$PARENT)))
#data3.pat2=list(y=data3.pat$nrec,cohort=as.factor(data3.pat$cohort),age=as.numeric(data3.pat$age.at.birth),J=nrow(data3.pat),I=length(unique(data3.pat$PARENT)))

