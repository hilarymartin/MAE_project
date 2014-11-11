library(rstan)
library(scales)
library(parallel)
load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.RData")
argv <- commandArgs(trailingOnly = TRUE)
print(argv)    
### Model 1.
#for informative meioses only with nkid=>2, all cohorts, fit cohort effect and family effect
#cohort-specific age effect
mycols=c("black","blue","red","green","orange","red4","purple","darkgreen")
names(mycols)=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB")
if(argv[1]==1){
#model1.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.stan",data=data1.mat2,iter=10000,chains=4)
#print(model1.mat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.mat.RData")
load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.mat.RData")
mysim<-extract(model1.mat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T)
my.ylim=7
pdf("RSTAN_output/model1.maternal.posteriors.beta_Age.pdf",height=5,width=5)

for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
abline(v=0,lwd=2)
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")
legend("topleft",as.character(cohort.codes[order(cohort.codes[,1]),2]),col=mycols[as.character(cohort.codes[,2])],lty=1,lwd=2,cex=0.5)
legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
}
dev.off()
}else if(argv[1] ==2){
#model1.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.stan",data=data1.pat2,iter=10000,chains=4)
#print(model1.pat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.pat.RData")
load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.pat.RData")
mysim<-extract(model1.pat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
my.ylim=15
pdf("RSTAN_output/model1.paternal.posteriors.beta_Age.pdf",height=5,width=5)

for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
legend("topleft",as.character(cohort.codes[order(cohort.codes[,1]),2]),col=mycols[as.character(cohort.codes[,2])],lty=1,lwd=2,cex=0.5)
abline(v=0,lwd=2)
}
dev.off()

}else if(argv[1] ==3){
#common age effect
#model1.2.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.2.stan",data=data1.mat2,iter=10000,chains=4)
#print(model1.2.mat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.2.mat.RData")
load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.2.mat.RData")
mysim<-extract(model1.2.mat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T)
my.ylim=12
pdf("RSTAN_output/model1.2.maternal.posteriors.beta_Age.pdf",height=5,width=5)
    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
abline(v=0,lwd=2)
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")
legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
dev.off()

    
}else if(argv[1] ==4){
#model1.2.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.2.stan",data=data1.pat2,iter=10000,chains=4)
#print(model1.2.pat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.2.pat.RData")
load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.2.pat.RData")
mysim<-extract(model1.2.pat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
my.ylim=12
pdf("RSTAN_output/model1.2.paternal.posteriors.beta_Age.pdf",height=5,width=5)
    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col="black",lwd=2)
abline(v=0,lwd=2)

dev.off()

}else if(argv[1] == 9){
#cohort-specific age effects, but drawn from N(beta_global,sigmasq_global); with cohort effect
#model1.3.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.stan",data=data1.mat2,iter=10000,chains=4)
#print(model1.3.mat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.3.mat.RData")
load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.3.mat.RData")
#print(model1.3.mat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
mysim<-extract(model1.3.mat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T)
pdf("RSTAN_output/model1.3.maternal.posteriors.beta_Age.pdf",height=5,width=5)
my.ylim=7
for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
}
dev.off()
pdf("RSTAN_output/model1.3.maternal.posteriors.sigmasq_global.pdf",height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()

            
}else if(argv[1] ==10){
#model1.3.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.stan",data=data1.pat2,iter=10000,chains=4)
#print(model1.3.pat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.3.pat.RData")

load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.3.pat.RData")
#print(model1.3.mat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
mysim<-extract(model1.3.pat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
my.ylim=14
pdf("RSTAN_output/model1.3.paternal.posteriors.beta_Age.pdf",height=5,width=5)


for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
}
dev.off()
pdf("RSTAN_output/model1.3.paternal.posteriors.sigmasq_global.pdf",height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()

    
}else if(argv[1] == 11){
#cohort-specific age effects, but drawn from N(beta_global,sigmasq_global); no cohort effect
#model1.4.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.4.stan",data=data1.mat2,iter=10000,chains=4)
#print(model1.4.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.4.mat.RData")
load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.4.mat.RData")

mysim<-extract(model1.4.mat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T)
my.ylim=12
pdf("RSTAN_output/model1.4.maternal.posteriors.beta_Age.pdf",height=5,width=5)

for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
}
dev.off()
pdf("RSTAN_output/model1.4.maternal.posteriors.sigmasq_global.pdf",height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()



}else if(argv[1] ==12){
#model1.4.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.4.stan",data=data1.pat2,iter=10000,chains=4)
#print(model1.4.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.4.pat.RData")

load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.4.pat.RData")

mysim<-extract(model1.4.pat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
my.ylim=15
pdf("RSTAN_output/model1.4.paternal.posteriors.beta_Age.pdf",height=5,width=5)

for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
}
dev.off()
pdf("RSTAN_output/model1.4.paternal.posteriors.sigmasq_global.pdf",height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()



    
}else if(argv[1] ==5){
### Model 2.
#for informative meioses only, all cohorts, fit cohort effect, BUT NO FAMILY EFFECT
#cohort-specific age effect
#model2.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.stan",data=data1.2.mat2,iter=10000,chains=4)
#print(model2.mat,pars=c("beta_Age","beta_Cohort","tausq","mu_m"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.mat.RData")

load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.mat.RData")
mysim<-extract(model2.mat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_2.txt",header=T)
my.ylim=8
pdf("RSTAN_output/model2.maternal.posteriors.beta_Age.pdf",height=5,width=5)

for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
abline(v=0,lwd=2)
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")
legend("topleft",as.character(cohort.codes[order(cohort.codes[,1]),2]),col=mycols[as.character(cohort.codes[,2])],lty=1,lwd=2,cex=0.5)
legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
}
dev.off()

    
}else if(argv[1] ==6){
#model2.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.stan",data=data1.2.pat2,iter=10000,chains=4)
#print(model2.pat,pars=c("beta_Age","beta_Cohort","tausq","mu_m"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.pat.RData")

load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.pat.RData")
mysim<-extract(model2.pat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_2.txt",header=T)
my.ylim=15
pdf("RSTAN_output/model2.paternal.posteriors.beta_Age.pdf",height=5,width=5)

for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
legend("topleft",as.character(cohort.codes[order(cohort.codes[,1]),2]),col=mycols[as.character(cohort.codes[,2])],lty=1,lwd=2,cex=0.5)
abline(v=0,lwd=2)
}
dev.off()


    
}else if(argv[1] ==7){
#common age effect
#model2.2.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.2.stan",data=data1.2.mat2,iter=10000,chains=4)
#print(model2.2.mat,pars=c("beta_Age","beta_Cohort","tausq","mu_m"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.2.mat.RData")

load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.2.mat.RData")
mysim<-extract(model2.2.mat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_2.txt",header=T)
pdf("RSTAN_output/model2.2.maternal.posteriors.beta_Age.pdf",height=5,width=5)
my.ylim=12

    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
abline(v=0,lwd=2)
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")

legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
dev.off()

   
}else if(argv[1] ==8){
#model2.2.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.2.stan",data=data1.2.pat2,iter=10000,chains=4)
#print(model2.2.pat,pars=c("beta_Age","beta_Cohort","tausq","mu_m"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.2.pat.RData")

load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.2.pat.RData")
mysim<-extract(model2.2.pat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_2.txt",header=T)
pdf("RSTAN_output/model2.2.paternal.posteriors.beta_Age.pdf",height=5,width=5)
my.ylim=12

    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col="black",lwd=2)

abline(v=0,lwd=2)

dev.off()

}

#### Model 3.
#for nkid=>2, family effect, and fit N_called_crossovers~Binom(Y,p=f(cohort*family_type)); only for those with 
#data2.mat2=list(y=data2.mat$nrec,cohort=as.factor(data2.mat$cohort),family=as.factor(data2.mat$PARENT),age=as.numeric(data2.mat$age.at.birth),J=nrow(data2.mat),I=length(unique(data2.mat$PARENT)))
#data2.pat2=list(y=data2.pat$nrec,cohort=as.factor(data2.pat$cohort),family=as.factor(data2.pat$PARENT),age=as.numeric(data2.pat$age.at.birth),J=nrow(data2.pat),I=length(unique(data2.pat$PARENT)))
if(argv[1] >=13 & argv[1] <17){
###to run these on cluster1, need to run each chain separately, save, then load in a new R session as a list, and combine in single stanfit object using sflist2stanfit
model3.mat.chain = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.stan",data=data2.mat2,iter=10000,chains=1,
    pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m"))
####N.B. chain should be 14-as.numeric(argv[1])
save(model3.mat.chain ,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.v2.chain",as.numeric(argv[1])-13+1,".mat.RData"))
}else if(argv[1] >=17 & argv[1] <=20){
model3.pat.chain = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.stan",data=data2.pat2,iter=10000,chains=1,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m"))
save(model3.pat.chain ,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.v2.chain",as.numeric(argv[1])-17+1,".pat.RData"))
#model3.pat <- sflist2stanfit(model3.pat.list)
#print(model3.pat,pars=c("beta_Age","p_by_cohort","mu_m","omega"),digits=4)
#save.image("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model3.pat.RData")
}

#if(argv[1]==22|argv[1]==24|argv[1]==25|argv[1]==26){
if(argv[1]==25|argv[1]==26){
#combine maternal model 3
if(argv[1]==22){    
    model3.mat.list=list()
    for(i in 1:4){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.v2.chain",i,".mat.RData"))
        model3.mat.list[[i]]=model3.mat.chain
        rm(model3.mat.chain)
    }
    model3.mat <- sflist2stanfit(model3.mat.list)
    print(model3.mat,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)

    mysim<-extract(model3.mat,permuted=T)
    parent="maternal"
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.txt",header=T)[,1:5]
}
 
if(argv[1]==25){
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.adjusted_priors.all_chains.mat.including_alpha.RData")
    print(model3.mat,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
    mysim<-extract(model3.mat,permuted=T)
    parent="adjusted_priors.maternal"
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.txt",header=T)[,1:5]
}

if(argv[1]==26){
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.adjusted_priors.all_chains.pat.including_alpha.RData")
    print(model3.pat,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
    mysim<-extract(model3.pat,permuted=T)
    parent="adjusted_priors.paternal"
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.txt",header=T)[,1:5]
}

if(argv[1]==24){    
    model3.pat.list=list()
    for(i in 1:4){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.v2.chain",i,".pat.RData"))
        model3.pat.list[[i]]=model3.pat.chain
        rm(model3.pat.chain)
    }
    model3.pat <- sflist2stanfit(model3.pat.list)
    print(model3.pat,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
    mysim<-extract(model3.pat,permuted=T)
    parent="paternal"
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.txt",header=T)[,1:5]
}

    cohort.codes=cohort.codes[order(cohort.codes[,1]),]
    cohort.codes$Pop=unlist(lapply(strsplit(cohort.codes$Cohort,".",fixed=T),function(x){return(x[[1]])}))
    cohort.codes$Fam.type=sapply(1:nrow(cohort.codes),function(x){gsub(paste0(cohort.codes$Pop[x],"."),"",cohort.codes$Cohort[x])})
    beta.parameters=beta.parameters[order(beta.parameters$Cohort.type),]
    mylty=c(1,2,4,1,2,4)
    names(mylty)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
    mytype=c("l","l","l","b","b","b")
    names(mytype)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
    mylwd=c(2,2,2,1,1,1)
    names(mylwd)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
        mypch=c(NA,NA,NA,19,19,19)
    names(mypch)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")

    pdf(paste0("RSTAN_output/model3.",parent,".posteriors.beta_Age.pdf"),height=5,width=5)
    my.ylim=27
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
            plot(density(mysim$beta_Age[,i]),xlim=c(-0.2,0.2),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],ylim=c(0,my.ylim),
                 lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])])
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,"Pop"])],lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])])
        }
        lines(density(mysim$beta_global),col="black",lwd=3,lty=3)
          legend("topleft",c(names(mycols),"global","informative, 2 generations","informative, 3 generations","uninformative, 2 kids","both parents","one parent"),col=c(mycols,rep("black",6)),lty=c(rep(1,length(mycols)),3,1,2,4,1,1),
                 lwd=c(rep(2,length(mycols)),3,2,2,2,2,1),cex=0.5)
        abline(v=0,lwd=2)
        if(FALSE & parent=="maternal"){
            polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
            abline(v=0.067,lty=2,lwd=2)
            polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
            abline(v=0.082,lty=3,lwd=2)
            polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
            abline(v=0.19,lty=4,lwd=2)
            abline(v=-0.42,lty=4,lwd=2,col="grey")
            legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
        }
    }
    dev.off()

    pdf(paste0("RSTAN_output/model3.",parent,".posteriors.p.pdf"),height=5,width=5)
    for(i in 1:ncol(mysim$p_by_cohort)){
        mylabels=c(paste(c("informative, 2 generations","informative, 3 generations","uninformative, 2 kids"),"both parents",sep=", "),paste(c("informative, 2 generations","informative, 3 generations","uninformative, 2 kids"),"1 parent",sep=", "))
        names(mylabels)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
        sample.size= sum(data2.mat$cohort.family.type==cohort.codes[cohort.codes[,1]==i,"Cohort"])
        plot(density(mysim$p_by_cohort[,i]),xlab="p",col=mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],xlim=c(0,1),
             lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],main=paste0(cohort.codes[cohort.codes[,1]==i,"Pop"],", ",mylabels[cohort.codes[cohort.codes[,1]==i,"Fam.type"]]))
        curve(dbeta(x,beta.parameters[beta.parameters$Cohort.type==i,"alpha"],beta.parameters[beta.parameters$Cohort.type==i,"beta"]),add=T)
        abline(v=median(mysim$p_by_cohort[,i]),col=mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])])
        abline(v=beta.parameters[beta.parameters$Cohort.type==i,"mean.relative.to.all.informative"],col="black")
        legend("topleft",c("posterior","prior",paste0("n = ",sample.size)),col=c(mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],"black",NA),lty=c(mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],1,NA),
               lwd=c(mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],1,NA))
    }
dev.off()
                                                                                                       
    pdf(paste0("RSTAN_output/model3.",parent,".posteriors.sigmasq_global.pdf"),height=5,width=5)
      plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
    dev.off()
    
}

####model 3.2
if(argv[1] ==21){
cat("Model 3.2 maternal\n")
#library(parallel)
#model3.2.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m"))
#model3.2.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.mat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m")))
#model3.2.mat <- sflist2stanfit(model3.2.mat.list)
#save(model3.2.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.2.all_chains.mat.RData"))
}

if(argv[1] ==23){
cat("Model 3.2 paternal\n")
#library(parallel)
#model3.2.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m"))
#model3.2.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m")))
#model3.2.pat <- sflist2stanfit(model3.2.pat.list)
#save(model3.2.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.2.all_chains.pat.RData"))
}

if(FALSE){
if(argv[1] ==25){
cat("Model 3 maternal, with adjusted priors\n")
library(parallel)
#model3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m"))
model3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
model3.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.mat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
model3.mat <- sflist2stanfit(model3.mat.list)
#save(model3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.adjusted_priors.all_chains.mat.RData"))
save(model3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.adjusted_priors.all_chains.mat.including_alpha.RData"))
}

if(argv[1] ==26){
#cat("Model 3 paternal, with adjusted priors\n")
library(parallel)
#model3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m"))
model3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
model3.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
model3.pat <- sflist2stanfit(model3.pat.list)
#save(model3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.adjusted_priors.all_chains.pat.RData"))
save(model3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.adjusted_priors.all_chains.pat.including_alpha.RData"))
}

if(argv[1] ==27){
cat("Model 3.2 maternal, with adjusted priors\n")
library(parallel)
#model3.2.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m"))
model3.2.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
#model3.2.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.mat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m")))
model3.2.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.mat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
model3.2.mat <- sflist2stanfit(model3.2.mat.list)
#save(model3.2.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.2.adjusted_priors.all_chains.mat.RData"))
save(model3.2.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.2.adjusted_priors.all_chains.mat.including_alpha.RData"))
}

if(argv[1] ==28){
cat("Model 3.2 paternal, with adjusted priors\n")
library(parallel)
#model3.2.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m"))
model3.2.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
#model3.2.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m")))
model3.2.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
model3.2.pat <- sflist2stanfit(model3.2.pat.list)
#save(model3.2.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.2.adjusted_priors.all_chains.pat.RData"))
save(model3.2.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.2.adjusted_priors.all_chains.pat.including_alpha.RData"))
}
}

#if(argv[1] %in% c(21,23,27,28)){
if(argv[1] %in% c(27,28)){
    if(argv[1]==21){
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.2.all_chains.mat.RData")
        print(model3.2.mat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
        mysim<-extract(model3.2.mat,permuted=T)
        parent="maternal"
        cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.txt",header=T,stringsAsFactors=F)
        beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.txt",header=T)[,1:5]
    }
    if(argv[1]==23){
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.2.all_chains.pat.RData")
        print(model3.2.pat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
        mysim<-extract(model3.2.pat,permuted=T)
        parent="paternal"
        cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.txt",header=T,stringsAsFactors=F)
        beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.txt",header=T)[,1:5]
    }
    if(argv[1]==27){
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.2.adjusted_priors.all_chains.mat.including_alpha.RData")
        print(model3.2.mat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
        mysim<-extract(model3.2.mat,permuted=T)
        parent="adjusted_priors.maternal"
        cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.txt",header=T,stringsAsFactors=F)
        beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.txt",header=T)[,1:5]
    }
    if(argv[1]==28){
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.2.adjusted_priors.all_chains.pat.including_alpha.RData")
        print(model3.2.pat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
        mysim<-extract(model3.2.pat,permuted=T)
        parent="adjusted_priors.paternal"
        cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.txt",header=T,stringsAsFactors=F)
        beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.txt",header=T)[,1:5]
    }

    cohort.codes=cohort.codes[order(cohort.codes[,1]),]
    cohort.codes$Pop=unlist(lapply(strsplit(cohort.codes$Cohort,".",fixed=T),function(x){return(x[[1]])}))
    cohort.codes$Fam.type=sapply(1:nrow(cohort.codes),function(x){gsub(paste0(cohort.codes$Pop[x],"."),"",cohort.codes$Cohort[x])})
    beta.parameters=beta.parameters[order(beta.parameters$Cohort.type),]
    mylty=c(1,2,4,1,2,4)
    names(mylty)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
    mytype=c("l","l","l","b","b","b")
    names(mytype)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
    mylwd=c(2,2,2,1,1,1)
    names(mylwd)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
        mypch=c(NA,NA,NA,19,19,19)
    names(mypch)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")

    pdf(paste0("RSTAN_output/model3.2.",parent,".posteriors.beta_Age.pdf"),height=5,width=5)
    plot(density(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age")
#    legend("topleft",c(names(mycols),"global","informative, 2 generations","informative, 3 generations","uninformative, 2 kids","both parents","one parent"),col=c(mycols,rep("black",6)),lty=c(rep(1,length(mycols)),3,1,2,4,1,1),
#                 lwd=c(rep(2,length(mycols)),3,2,2,2,2,1),cex=0.5)
    abline(v=0,lwd=2)
    dev.off()

    pdf(paste0("RSTAN_output/model3.2.",parent,".posteriors.p.pdf"),height=5,width=5)
    for(i in 1:ncol(mysim$p_by_cohort)){
        mylabels=c(paste(c("informative, 2 generations","informative, 3 generations","uninformative, 2 kids"),"both parents",sep=", "),paste(c("informative, 2 generations","informative, 3 generations","uninformative, 2 kids"),"1 parent",sep=", "))
        names(mylabels)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
        sample.size= sum(data2.mat$cohort.family.type==cohort.codes[cohort.codes[,1]==i,"Cohort"])
        plot(density(mysim$p_by_cohort[,i]),xlab="p",col=mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],xlim=c(0,1),
             lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],main=paste0(cohort.codes[cohort.codes[,1]==i,"Pop"],", ",mylabels[cohort.codes[cohort.codes[,1]==i,"Fam.type"]]))
        curve(dbeta(x,beta.parameters[beta.parameters$Cohort.type==i,"alpha"],beta.parameters[beta.parameters$Cohort.type==i,"beta"]),add=T)
        abline(v=median(mysim$p_by_cohort[,i]),col=mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])])
        abline(v=beta.parameters[beta.parameters$Cohort.type==i,"mean.relative.to.all.informative"],col="black")
        legend("topleft",c("posterior","prior",paste0("n = ",sample.size)),col=c(mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],"black",NA),lty=c(mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],1,NA),
               lwd=c(mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],1,NA))
    }
    dev.off()

    
    pdf(paste0("RSTAN_output/model3.2.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
#if(myname=="maternal"){my.ylim=0.3}else {my.ylim=0.5}
    plot(density(mysim$mu_m),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",lwd=2)#,ylim=c(0,my.ylim))
    curve(dnorm(x,36,sqrt(6)),lty=2,lwd=2,add=T)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()
library(pscl)
    pdf(paste0("RSTAN_output/model3.2.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
#if(myname=="maternal"){my.ylim=0.3}else {my.ylim=0.5}
    plot(density(mysim$sigmasq_m),xlim=range(mysim$sigmasq_m),xlab="sigmasq_m",main="Posterior for sigmasq_m",lwd=2)#,ylim=c(0,my.ylim))
    lines(density(rigamma(10000,5,5)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

        pdf(paste0("RSTAN_output/model3.2.",myname,".posteriors.omega.pdf"),height=5,width=5)
    plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
    lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()
 
dev.off()


 
}

if(argv[1] ==29 | argv[1]==30){
#### Model 1.6 -- beta_Age drawn from distribution; alphas drawn from N(mu_cohort,sigmasq)
if(argv[1]==29){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
#    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.mat2,iter=10000,chains=0)
#    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
#    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
   
#    print(model1.6.mat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
#    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.6.all_chains.mat.including_alpha.RData"))
     load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.6.all_chains.mat.including_alpha.RData")


    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
#    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.pat2,iter=10000,chains=0)
#    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
#    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
#    print(model1.6.pat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
#    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.6.all_chains.pat.including_alpha.RData"))
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.6.all_chains.pat.including_alpha.RData")

    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output/model1.6.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
if(myname=="maternal"){my.ylim=7}else {my.ylim=12}
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
   lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
   abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
   if(myname=="maternal"){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
   }
    dev.off()
    
    pdf(paste0("RSTAN_output/model1.6.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
if(myname=="maternal"){my.ylim=0.3}else {my.ylim=0.5}
for(i in 1:ncol(mysim$mu_m)){
    if(i==1){
        plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()


pdf(paste0("RSTAN_output/model1.6.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}


if(argv[1] ==31 | argv[1]==32){
#### Model 1.7 -- beta_Age drawn from distribution; with family effect and cohort effect; using t distribution for Y and prior on beta_Age
if(argv[1]==31){
    data1.mat2$df_y <- 5
    data1.mat2$df_beta_age <- 5

#    model1.7.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.7.stan",data=data1.mat2,iter=10000,chains=0)
#    model1.7.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.7.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
#    model1.7.mat <- sflist2stanfit(model1.7.mat.list)
#    print(model1.7.mat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
#    save(model1.7.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.7.all_chains.mat.including_alpha.RData"))
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.7.all_chains.mat.including_alpha.RData")

    mysim<-extract(model1.7.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
    data1.pat2$df_y <- 5
    data1.pat2$df_beta_age <- 5

#    model1.7.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.7.stan",data=data1.pat2,iter=10000,chains=0)
#    model1.7.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.7.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
#    model1.7.pat <- sflist2stanfit(model1.7.pat.list)
#    print(model1.7.pat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
#    save(model1.7.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.7.all_chains.pat.including_alpha.RData"))
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.7.all_chains.pat.including_alpha.RData")

    mysim<-extract(model1.7.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output/model1.7.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
for(i in 1:ncol(mysim$beta_Age)){
    if(myname=="maternal"){my.ylim=7}else{my.ylim=12}
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
   if(myname=="maternal"){
       my.ylim=20
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
   }
    dev.off()
    
pdf(paste0("RSTAN_output/model1.7.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}

if(argv[1] ==33 | argv[1]==34){
#### Model 1.5 -- negative binomial model for informative families only
if(argv[1]==33){
#    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.stan",data=data1.mat2,iter=10000,chains=0)
#    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
#    print(model1.5.mat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
#    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.all_chains.mat.including_alpha.RData"))
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.all_chains.mat.including_alpha.RData")

    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
#    model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.stan",data=data1.pat2,iter=10000,chains=0)
#    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.pat <- sflist2stanfit(model1.5.pat.list)
#    print(model1.5.pat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
#    save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.all_chains.pat.including_alpha.RData"))
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.all_chains.pat.including_alpha.RData")
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output/model1.5.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
my.ylim=55
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
dev.off()
    

pdf(paste0("RSTAN_output/model1.5.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}

if(argv[1] ==35 | argv[1]==36){
#### Model 1.5, with different link function    -- negative binomial model for informative families only, with different link function
    cat("#### Model 1.5, with different link function\n")

if(argv[1]==35){
#    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.stan",data=data1.mat2,iter=10000,chains=0)
#    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
##    print(model1.5.mat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
 #   save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.different_link.all_chains.mat.including_alpha.RData"))
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.different_link.all_chains.mat.including_alpha.RData")

    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
#    model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.stan",data=data1.pat2,iter=10000,chains=0)
#    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.pat <- sflist2stanfit(model1.5.pat.list)
 #    print(model1.5.pat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
 #   save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.different_link.all_chains.pat.including_alpha.RData"))
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.different_link.all_chains.pat.including_alpha.RData")
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output/model1.5.different_link.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=200}else{my.ylim=300}
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
            plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
    abline(v=0,lwd=2)
    dev.off()
    

pdf(paste0("RSTAN_output/model1.5.different_link.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}

if(argv[1] ==37 | argv[1]==38){
#if(argv[1] ==37 ){
#### Model 1.5, with different link function  and weaker priors  -- negative binomial model for informative families only, with different link function
    cat("#### Model 1.5, with different link function\n")
if(argv[1]==37){
#    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.weaker_priors.stan",data=data1.mat2,iter=10000,chains=0)
#    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
 #    print(model1.5.mat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
 #   save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.different_link.weaker_priors.all_chains.mat.including_alpha.RData"))
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.different_link.weaker_priors.all_chains.mat.including_alpha.RData")

    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
####Need to run this for longer - didn't converge
#    model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.weaker_priors.stan",data=data1.pat2,iter=10000,chains=0)
#    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=20000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.pat <- sflist2stanfit(model1.5.pat.list)
#    
#    print(model1.5.pat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
 #   save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.different_link.weaker_priors.all_chains.pat.including_alpha.RData"))
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.5.different_link.weaker_priors.all_chains.pat.including_alpha.RData")
    
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output/model1.5.different_link.weaker_priors.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=200}else{my.ylim=200}
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
   lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
   abline(v=0,lwd=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()
    
pdf(paste0("RSTAN_output/model1.5.different_link.weaker_priors.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}

####################################################

if(argv[1] ==39 | argv[1]==40){
#### Model 1.3, weaker priors on age effect -- beta_Age drawn from distribution
cat("    #### Model 1.3, weaker priors on age effect -- beta_Age drawn from distribution\n")
if(argv[1]==39){
    model1.3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.weaker_priors_on_age_effect.stan",data=data1.mat2,iter=10000,chains=0)
    model1.3.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.3.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
    model1.3.mat <- sflist2stanfit(model1.3.mat.list)
    print(model1.3.mat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.3.weaker_priors_on_age_effect.all_chains.mat.including_alpha.RData"))
#    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.3.weaker_priors_on_age_effect.all_chains.mat.including_alpha.RData")

    mysim<-extract(model1.3.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
    model1.3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.weaker_priors_on_age_effect.stan",data=data1.pat2,iter=10000,chains=0)
    model1.3.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.3.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
    model1.3.pat <- sflist2stanfit(model1.3.pat.list)
    print(model1.3.pat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.3.weaker_priors_on_age_effect.all_chains.pat.including_alpha.RData"))
#    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.3.weaker_priors_on_age_effect.all_chains.pat.including_alpha.RData")
    mysim<-extract(model1.3.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output/model1.3.weaker_priors_on_age_effect.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
for(i in 1:ncol(mysim$beta_Age)){
if(myname=="maternal"){    my.ylim=7}else {my.ylim=12}
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
if(myname=="maternal"){
    polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")

legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
}
dev.off()

pdf(paste0("RSTAN_output/model1.3.weaker_priors_on_age_effect.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()

}


if(argv[1] ==41 | argv[1]==42){
#### Model 1.3, uniform_prior_on_sigmasq_m. -- beta_Age drawn from distribution
cat(" #### Model 1.3, uniform_prior_on_sigmasq_m. -- beta_Age drawn from distribution\n")
if(argv[1]==41){
    model1.3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.uniform_prior_on_sigmasq_m.stan",data=data1.mat2,iter=10000,chains=0)
    model1.3.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.3.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
    model1.3.mat <- sflist2stanfit(model1.3.mat.list)
    print(model1.3.mat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.3.uniform_prior_on_sigmasq_m.all_chains.mat.including_alpha.RData"))
#    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.3.uniform_prior_on_sigmasq_m.all_chains.mat.including_alpha.RData")

    mysim<-extract(model1.3.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
    model1.3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.uniform_prior_on_sigmasq_m.stan",data=data1.pat2,iter=10000,chains=0)
    model1.3.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.3.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
    model1.3.pat <- sflist2stanfit(model1.3.pat.list)
    print(model1.3.pat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.3.uniform_prior_on_sigmasq_m.all_chains.pat.including_alpha.RData"))
#    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.modell1.3.uniform_prior_on_sigmasq_m.all_chains.pat.including_alpha.RData")
    mysim<-extract(model1.3.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output/model1.3.uniform_prior_on_sigmasq_m.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
for(i in 1:ncol(mysim$beta_Age)){
if(myname=="maternal"){    my.ylim=7} else {my.ylim=15}
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
if(myname=="maternal"){
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")

legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
}
dev.off()

pdf(paste0("RSTAN_output/model1.3.uniform_prior_on_sigmasq_m.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()

pdf(paste0("RSTAN_output/model1.3.uniform_prior_on_sigmasq_m.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_m",main="Posterior for sigmasq_m",lwd=2)
dev.off()

}

if(FALSE){
###################### HAVEN'T FINISHED THIS SECTION OF CODE - MAYBE VARYING P_BY_COHORT IS NOT SO IMPORTANT
    
if(argv[1] ==43 | argv[1]==44){
    #Model 3, uniform_p_priors, with adjusted priors
cat("Model 3, uniform_p_priors, with adjusted priors\n")
    if(argv[1] ==43){
        model3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.uniform_p_priors.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
        model3.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.mat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
        model3.mat <- sflist2stanfit(model3.mat.list)
        save(model3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.uniform_p_priors.all_chains.mat.including_alpha.RData"))
        cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.txt",header=T)
    }else {
    model3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.uniform_p_priors.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
    model3.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
    model3.pat <- sflist2stanfit(model3.pat.list)
    save(model3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output/RSTAN.model3.uniform_p_priors.all_chains.pat.including_alpha.RData"))
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.txt",header=T)
}

pdf("RSTAN_output/model1.3.maternal.posteriors.beta_Age.pdf",height=5,width=5)
my.ylim=7
for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
}
dev.off()
pdf("RSTAN_output/model1.3.maternal.posteriors.sigmasq_global.pdf",height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()

}
}
### Model 4. N.B. NOT SURE SHOULD INCLUDE FAMILIES WITH ONLY 1 KID - DOESN'T REALLY LOOK BINOMIAL
#for all meioses, NO FAMILY EFFECT; fit N_called_crossovers~Binom(Y,p=f(cohort*family_type)); only for those with
#data3.mat2=list(y=data3.mat$nrec,cohort=as.factor(data3.mat$cohort),age=as.numeric(data3.mat$age.at.birth),J=nrow(data3.mat),I=length(unique(data3.mat$PARENT)))
#data3.pat2=list(y=data3.pat$nrec,cohort=as.factor(data3.pat$cohort),age=as.numeric(data3.pat$age.at.birth),J=nrow(data3.pat),I=length(unique(data3.pat$PARENT)))

