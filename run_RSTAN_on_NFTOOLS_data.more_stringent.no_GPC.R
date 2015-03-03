library(rstan)
library(scales)
library(parallel)
library(pscl)
load("/well/donnelly/hilary/maternal_age_and_recombination/NFTOOLS_data_for_RSTAN.more_stringent.no_GPC.RData")
argv <- commandArgs(trailingOnly = TRUE)
print(argv)    
mycols=c("black","blue","red","green","orange","red4","purple","darkgreen","violetred2")
names(mycols)=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB","ORCADES")
print.only=as.logical(as.numeric(argv[2]))


if(argv[1] ==29 | argv[1]==30){
#### Model 1.6 -- beta_Age drawn from distribution; alphas drawn from N(mu_cohort,sigmasq)
if(argv[1]==29){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family

    mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=100
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    
if(!print.only){
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.mat2,iter=10000,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
 print(model1.6.mat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
   mean_alpha_prior=27
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=64
    data1.pat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
 
if(!print.only){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.pat2,iter=10000,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
 load( paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))

}
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
if(myname=="maternal"){my.ylim=6}else {my.ylim=12}
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=range(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
   lines(density(mysim$beta_global),col="black",lwd=3,lty=4)
curve(dnorm(x,0,1),add=T,lty=2,lwd=2)
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
              abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
              legend("topright",c("Campbell","Kong","Coop","Hussin","Bleazard"),lty=c(2,3,4,4,4),lwd=2,col=c("black","black","black","grey","darkgreen"),cex=0.5)


   }
    dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.mu_m.pdf"),height=5,width=5)
if(myname=="maternal"){my.ylim=0.2}else {my.ylim=0.5}
for(i in 1:ncol(mysim$mu_m)){
    if(i==1){
        plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.15
        my.xlim=100
    }else {
        my.ylim=0.6
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i],n=256),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i],n=256),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    emp.variances=read.delim("variance_in_mean_duoHMM_count_by_parent.informative_meioses_only.min_2_kids.txt",header=T,stringsAsFactors=F)
    for(z in 1:nrow(cohort.codes)){
        if(argv[1]==29){
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"maternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        } else {
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"paternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        }
    }
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior","empirical variance in parent means"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,2),lwd=c(rep(2,nrow(cohort.codes)),3,1),cex=0.5)
    dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
    lines(density(rigamma(10000,3,0.5)),lwd=2,lty=2)
legend("topright",c("posterior","prior"),lty=c(1,2),lwd=c(2,2))
dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.tausq.pdf"),height=5,width=5)
plot(density(mysim$tausq),xlab="tausq",main="Posterior for tausq",lwd=2)
lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()

}

if(argv[1] ==33 | argv[1]==34){
#### Model 1.5 -- negative binomial model for informative families only
if(argv[1]==33){
  mu_alpha_prior = 37
  sigmasq_alpha_prior = 6
  sigmasq_m_alpha = 5
  sigmasq_m_beta = 5
  data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.mat2$sigmasq_m_beta = sigmasq_m_beta
  data1.mat2$mean_alpha_prior = mu_alpha_prior
  data1.mat2$sigmasq_alpha_prior = sigmasq_alpha_prior

  cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
      cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
  data1.mat2$cohort_by_family = cohort_by_family
  
    
if(!print.only ){
    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.stan",data=data1.mat2,iter=10000,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
print(model1.5.mat,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.all_chains.mat.including_alpha.RData"))
} else {
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.all_chains.mat.including_alpha.RData")
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
  mu_alpha_prior = 32
  sigmasq_alpha_prior = 6
  sigmasq_m_alpha = 5
  sigmasq_m_beta = 5
  data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.pat2$sigmasq_m_beta = sigmasq_m_beta
  data1.pat2$mean_alpha_prior = mu_alpha_prior
  data1.pat2$sigmasq_alpha_prior = sigmasq_alpha_prior
  cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
      cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
  data1.pat2$cohort_by_family = cohort_by_family

  if(!print.only){
    model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.stan",data=data1.pat2,iter=10000,chains=0)
    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.pat <- sflist2stanfit(model1.5.pat.list)
print(model1.5.pat,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.all_chains.pat.including_alpha.RData"))
} else {
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.all_chains.pat.including_alpha.RData")
}
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
my.ylim=18
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
lines(density(mysim$beta_global),col="black",lwd=3,lty=4)
curve(dnorm(x,0,sqrt(0.05)),add=T,lty=2,lwd=2)
abline(v=0,lwd=2)
dev.off()
    

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}

if(argv[1] ==35 | argv[1]==36){
#### Model 1.5, with different link function    -- negative binomial model for informative families only, with different link function
cat("#### Model 1.5, with different link function\n")
if(argv[1]==35){
if(!print.only){
    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.stan",data=data1.mat2,iter=10000,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
print(model1.5.mat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
   save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.different_link.all_chains.mat.including_alpha.RData"))
} else {
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.different_link.all_chains.mat.including_alpha.RData")
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
 if(!print.only){
     model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.stan",data=data1.pat2,iter=10000,chains=0)
    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.pat <- sflist2stanfit(model1.5.pat.list)
 print(model1.5.pat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.different_link.all_chains.pat.including_alpha.RData"))
 } else {
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.different_link.all_chains.pat.including_alpha.RData")
}
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.different_link.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
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
    

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.different_link.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}

if(argv[1] ==37 | argv[1]==38){
#### Model 1.5, with different link function  and weaker priors  -- negative binomial model for informative families only, with different link function
    cat("#### Model 1.5, with different link function\n")
if(argv[1]==37){
if(!print.only){
    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.weaker_priors.stan",data=data1.mat2,iter=10000,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
 print(model1.5.mat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.different_link.weaker_priors.all_chains.mat.including_alpha.RData"))
} else {
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.different_link.weaker_priors.all_chains.mat.including_alpha.RData")
}

    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
if(!print.only){
####Need to run this for longer - didn't converge
    model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.weaker_priors.stan",data=data1.pat2,iter=10000,chains=0)
    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.pat <- sflist2stanfit(model1.5.pat.list)
print(model1.5.pat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.different_link.weaker_priors.all_chains.pat.including_alpha.RData"))
} else {
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.different_link.weaker_priors.all_chains.pat.including_alpha.RData")
}
    
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.different_link.weaker_priors.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
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
    
pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.different_link.weaker_priors.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}


if(argv[1] ==43 | argv[1]==44){
#### Model 1.62 -- common beta_Age for all cohorts; N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort
    df.alpha=5
    if(argv[1]==43){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    data1.mat2$df_alpha = df.alpha

    mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=100
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
if(!print.only){
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.stan",data=data1.mat2,iter=10000,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,
                ".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha

     mean_alpha_prior=27
     data1.pat2$mean_alpha_prior = mean_alpha_prior
     sigmasq_m_alpha = 2
     sigmasq_m_beta = 15
     data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
     data1.pat2$sigmasq_m_beta = sigmasq_m_beta
     sigmasq_mu_m_prior=8^2
     data1.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior
     if(!print.only){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.stan",data=data1.pat2,iter=10000,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
     save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",
                 sigmasq_m_beta,".RData"))
 }
     mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}
if(FALSE){
pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".traceplots.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    traceplot(model1.6.mat,pars=c("beta_Age"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("mu_m"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("sigmasq_m"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("tausq"),window=c(10^2,10^4))
    dev.off()
}
    
pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    my.ylim=100
#    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
        plot(density(mysim$beta_Age),xlim=c(-0.6,0.6),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
    curve(dnorm(x,0,1),lty=2,lwd=2,add=T)
    legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2)
    if(myname=="maternal"){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
       legend("topright",c("Campbell","Kong","Coop","Hussin","Bleazard"),lty=c(2,3,4,4,4),lwd=2,col=c("black","black","black","grey","darkgreen"),cex=0.5)
        
   }
    dev.off()
    
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=0.3}else {my.ylim=0.7}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
#//        curve(dnorm(x,mean_alpha_prior,sqrt(5)),add=T,lwd=3,lty=2)
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.13
        my.xlim=100
    }else {
        my.ylim=0.7
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    emp.variances=read.delim("variance_in_mean_duoHMM_count_by_parent.informative_meioses_only.min_2_kids.txt",header=T,stringsAsFactors=F)
    for(z in 1:nrow(cohort.codes)){
        if(argv[1]==43){
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"maternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        } else {
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"paternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        }
    }
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior","empirical variance in parent means"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,2),lwd=c(rep(2,nrow(cohort.codes)),3,1),cex=0.5)
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),lwd=2,main="Posterior for tausq",xlab="tausq")
    lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
    dev.off()
}

if(argv[1] ==45 | argv[1]==46){
#### Model 1.5 -- negative binomial model for informative families only; common beta_Age for all cohorts; alpha drawn from N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort     
if(argv[1]==45){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    mean_alpha_prior=36
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 5
    sigmasq_m_beta = 5
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
if(!print.only){
    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.stan",data=data1.mat2,iter=10000,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
print(model1.5.mat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.4.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.4.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    mean_alpha_prior=33
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 5
    sigmasq_m_beta = 5
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
    if(!print.only){
        model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.stan",data=data1.pat2,iter=10000,chains=0)
        model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
        model1.5.pat <- sflist2stanfit(model1.5.pat.list)
        print(model1.5.pat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
        save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.4.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.4.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.4.",myname,".posteriors.beta_Age.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
curve(dnorm(x,0,sqrt(0.05)),add=T,lwd=2,lty=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
abline(v=0,lwd=2)
dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.4.",myname,".posteriors.mu_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
      if(myname=="maternal"){my.ylim=1.6}else {my.ylim=1.6}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sd=sqrt(6)),lty=2,lwd=3,add=T)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()


    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.4.",myname,".posteriors.sigmasq_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)

if(myname=="maternal"){
    my.ylim=3
    my.xlim=4
}else{
    my.ylim=5
    my.xlim=4
}
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lty=2,lwd=3)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.4.",myname,".posteriors.omega.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()


}


if(argv[1] ==47 | argv[1]==48|argv[1]==55|argv[1]==56|argv[1]==57|argv[1]==58|argv[1]==59|argv[1]==60|argv[1]==61|argv[1]==62){
#### Model 1.62 -- common beta_Age for all cohorts; N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently; uniform priors on all parameters
    df.alpha=5
    if(argv[1]==47|argv[1]==55|argv[1]==57|argv[1]==59|argv[1]==61){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    data1.mat2$df_alpha = df.alpha
        mean_alpha_prior=38    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
if(!print.only){
    if(argv[1]==47){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors.stan",data=data1.mat2,iter=10000,chains=0)
    }
    if(argv[1]==55){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_beta_Age.stan",data=data1.mat2,iter=10000,chains=0)
    }
    if(argv[1]==57){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_mu_m_cohort.stan",data=data1.mat2,iter=10000,chains=0)
    }
    if(argv[1]==59){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_tau.stan",data=data1.mat2,iter=10000,chains=0)
    }
    if(argv[1]==61){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_sigma_m.stan",data=data1.mat2,iter=10000,chains=0)
    }
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=30000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)

    if(argv[1]==47){
                myname="uniform_priors.maternal"
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==55){
        myname="uniform_priors_on_beta_Age.maternal"
        save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==57){
        myname="uniform_priors_on_mu_m_cohort.maternal"
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==59){
        myname="uniform_priors_on_tau.maternal"
        save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==61){
        myname="uniform_priors_on_sigma_m.maternal"
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.mat.including_alpha.RData"))
    }
} else{
    if(argv[1]==47){
        myname="uniform_priors.maternal"
             load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors.all_chains.mat.including_alpha.RData")
        }
    if(argv[1]==55){
        myname="uniform_priors_on_beta_Age.maternal"
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==57){
        myname="uniform_priors_on_mu_m_cohort.maternal"
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.mat.including_alpha.RData"))
    }
      if(argv[1]==59){
          myname="uniform_priors_on_tau.maternal"
          load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.mat.including_alpha.RData"))
      }
        if(argv[1]==61){
            myname="uniform_priors_on_sigma_m.maternal"
            load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.mat.including_alpha.RData"))
        }
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)

} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha

        mean_alpha_prior=33
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
if(!print.only){
    if(argv[1]==48){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors.stan",data=data1.pat2,iter=10000,chains=0)
    }
    if(argv[1]==56){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_beta_Age.stan",data=data1.pat2,iter=10000,chains=0)
    }
    if(argv[1]==58){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_mu_m_cohort.stan",data=data1.pat2,iter=10000,chains=0)
    }
    if(argv[1]==60){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_tau.stan",data=data1.pat2,iter=10000,chains=0)
    }
    if(argv[1]==62){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_sigma_m.stan",data=data1.pat2,iter=10000,chains=0)
    }

    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=30000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)

    if(argv[1]==48){
        myname="uniform_priors.paternal"
    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors.v4.all_chains.pat.including_alpha.RData"))
                                        #    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors.all_chains.pat.including_alpha.RData")
#    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors.v4.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==56){
        myname="uniform_priors_on_beta_Age.paternal"
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==58){
        myname="uniform_priors_on_mu_m_cohort.paternal"
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==60){
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.pat.including_alpha.RData"))
        myname="uniform_priors_on_tau.paternal"
    }
    if(argv[1]==62){
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.pat.including_alpha.RData"))
        myname="uniform_priors_on_sigma_m.paternal"
    }

} else {

   if(argv[1]==48){
       myname="uniform_priors.paternal"
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors.all_chains.pat.including_alpha.RData")

    }
    if(argv[1]==56){
        myname="uniform_priors_on_beta_Age.paternal"
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==58){
        myname="uniform_priors_on_mu_m_cohort.paternal"
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==60){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.pat.including_alpha.RData"))
        myname="uniform_priors_on_tau.paternal"
    }
    if(argv[1]==62){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.pat.including_alpha.RData"))
        myname="uniform_priors_on_sigma_m.paternal"
    }

}
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
}
#pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.uniform_priors.v4.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
    my.ylim=100
    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
    #    if(argv[1] %in% c(47,48,55,56,57,58,59,60,61,62)){
    if(argv[1] %in% c(57,58,59,60,61,62)){
        curve(dnorm(x,0,1),add=T,lwd=2,lty=2)
        legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2)
    }
    if(argv[1] %in% c(47,55,57,59,61)){
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
    
        pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
        if(argv[1] %in% c(47,55,57,59,61)){
        my.ylim=0.6
        my.xlim=c(32,46)
    }else {
        my.ylim=0.7
        my.xlim=c(23,32)
    }

    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=my.xlim,xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    if(argv[1] %in% c(55,56,59,60,61,62)){
        curve(dnorm(x,mean_alpha_prior,sqrt(5)),add=T,lwd=2,lty=2)
    }
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),2),cex=0.5)
    dev.off()


        pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
        if(argv[1] %in% c(47,55,57,59,61)){
        my.ylim=0.2
        my.xlim=100
    }else {
        my.ylim=0.7
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    if(argv[1] %in% c(55,56,57,58,59,60)){
        lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    }

    emp.variances=read.delim("variance_in_mean_duoHMM_count_by_parent.informative_meioses_only.min_2_kids.txt",header=T,stringsAsFactors=F)
    for(z in 1:nrow(cohort.codes)){
        if(argv[1]==47|argv[1]==55|argv[1]==57|argv[1]==59|argv[1]==61){
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"maternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        } else {
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"paternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        }
    }

    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),2),cex=0.5)
    dev.off()


        pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.tausq.pdf"),height=5,width=5)
    plot(density(mysim$tausq),main="Posterior for tausq",lwd=2,xlab="tausq")
    if(argv[1] %in% c(55,56,57,58,61,62)){
        lines(density(rigamma(10000,2,70)),lwd=3,lty=2)
    }
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)

    dev.off()
    
}


if(argv[1] ==49 | argv[1]==50){
#### Model 1.5 -- negative binomial model for informative families only; common beta_Age for all cohorts; alpha drawn from N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently  for each cohort ; uniform priors
if(argv[1]==49){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
if(!print.only){
   model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.uniform_priors.stan",data=data1.mat2,iter=10000,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
print(model1.5.mat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.4.uniform_priors.all_chains.mat.including_alpha.RData"))
} else {
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.4.uniform_priors.all_chains.mat.including_alpha.RData")
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family

    data1.pat2$mean_alpha_prior = 33
    if(!print.only){
        model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.uniform_priors.stan",data=data1.pat2,iter=20000,chains=0)
        model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
        model1.5.pat <- sflist2stanfit(model1.5.pat.list)
        print(model1.5.pat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
        save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.4.uniform_priors.all_chains.pat.including_alpha.RData"))
    } else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.5.4.uniform_priors.v4.all_chains.pat.including_alpha.RData"))
}
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.4.uniform_priors.v4.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
abline(v=0,lwd=2)
dev.off()


pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.4.uniform_priors.v4.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
if(myname=="maternal"){
    my.ylim=1.7
} else {
    my.ylim=1.7
}
                                        #    if(myname=="maternal"){my.ylim=0.3}else {my.ylim=0.5}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.4.uniform_priors.v4.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=1
        my.xlim=10
    }else {
        my.ylim=10
        my.xlim=7
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2])),col=c(mycols[as.character(cohort.codes[,2])]),lty=c(rep(1,nrow(cohort.codes))),lwd=c(rep(2,nrow(cohort.codes))),cex=0.5)
dev.off()
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.5.4.uniform_priors.v4.",myname,".posteriors.omega.pdf"),height=5,width=5)
plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()


}


if(argv[1] ==51 | argv[1]==52){
#### Model 1.62 -- common beta_Age for all cohorts; N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort
    df.alpha=5
    if(argv[1]==51){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    data1.mat2$df_alpha = df.alpha
    mean_alpha_prior=38    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    
if(!print.only){
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas.stan",data=data1.mat2,iter=10000,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".RData"))
} else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha
    mean_alpha_prior=30
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
 if(!print.only){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas.stan",data=data1.pat2,iter=10000,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas.df_",df.alpha,".RData"))
}else{
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas.df_",df.alpha,".RData"))
 }
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.beta_Age.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
    my.ylim=100
    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
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
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.mu_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=0.6}else {my.ylim=0.7}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sqrt(5)),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.sigmasq_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.2
        my.xlim=100
    }else {
        my.ylim=0.7
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.tausq.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),xlab="tausq",main="Posterior for tausq",lwd=2)
    lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()
}

if(argv[1] ==53 | argv[1]==54){
#### Model 1.62 -- common beta_Age for all cohorts; N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort
    df.alpha=5
    if(argv[1]==53){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    data1.mat2$df_alpha = df.alpha

    ##mean_alpha_prior=38
    mean_alpha_prior=41
    sigmasq_mu_m_prior=100
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
if(!print.only){
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas_and_Y.stan",data=data1.mat2,iter=10000,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".RData"))
}else {
load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".RData"))
}
    
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha

    mean_alpha_prior=27
    sigmasq_mu_m_prior=64
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    data1.pat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior

     data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
 if(!print.only){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas_and_Y.stan",data=data1.pat2,iter=10000,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
         save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas_and_Y.df_",df.alpha,".RData"))
} else {
load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas_and_Y.df_",df.alpha,".RData"))
}
     mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
    my.ylim=100
    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
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
    
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=0.6}else {my.ylim=0.7}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.2
        my.xlim=100
    }else {
        my.ylim=0.7
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),xlab="tausq",main="Posterior for tausq",lwd=2)
    lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

}
if(argv[1] ==67 | argv[1]==68){
#### Model 1.62 -- common beta_Age for all cohorts; N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort
    df.alpha=5
    if(argv[1]==67){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    data1.mat2$df_alpha = df.alpha
##    mean_alpha_prior=38
    mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=100
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    if(!print.only){   
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.with_agesq_term.stan",data=data1.mat2,iter=10000,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=20000,pars=c("beta_Age","beta_Agesq","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    print(model1.6.mat,pars=c("beta_Age","beta_Agesq","tausq","sigmasq_m","mu_m"),digits=4)
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.with_agesq_term.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.with_agesq_term.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,                ".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="with_agesq_term.maternal"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha
     mean_alpha_prior=27
     data1.pat2$mean_alpha_prior = mean_alpha_prior
     sigmasq_m_alpha = 2
     sigmasq_m_beta = 15
     data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
     data1.pat2$sigmasq_m_beta = sigmasq_m_beta
     sigmasq_mu_m_prior=8^2
     data1.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior
     if(!print.only){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.with_agesq_term.stan",data=data1.pat2,iter=20000,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=20000,pars=c("beta_Age","beta_Agesq","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","beta_Agesq","tausq","sigmasq_m","mu_m"),digits=4)
     save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
 }
     mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="with_agesq_term.paternal"
}

pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    my.ylim=100
    plot(density(mysim$beta_Age),xlim=c(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
    curve(dnorm(x,0,1),lty=2,lwd=2,add=T)
    legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2)
    if(argv[1]==67){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
       legend("topright",c("Campbell","Kong","Coop","Hussin","Bleazard"),lty=c(2,3,4,4,4),lwd=2,col=c("black","black","black","grey","darkgreen"),cex=0.5)
   }
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.beta_Agesq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$beta_Agesq),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
    curve(dnorm(x,0,5),lty=2,lwd=2,add=T)
    legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2)
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=0.4}else {my.ylim=0.7}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }

    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.17
        my.xlim=100
    }else {
        my.ylim=0.7
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    emp.variances=read.delim("variance_in_mean_duoHMM_count_by_parent.informative_meioses_only.min_2_kids.txt",header=T,stringsAsFactors=F)
    for(z in 1:nrow(cohort.codes)){
        if(argv[1]==43){
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"maternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        } else {
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"paternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        }
    }
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior","empirical variance in parent means"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,2),lwd=c(rep(2,nrow(cohort.codes)),3,1),cex=0.5)
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/model1.6.2.",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),lwd=2,main="Posterior for tausq",xlab="tausq")
    lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
    dev.off()
}


if(argv[1] ==69 | argv[1]==70|argv[1] ==71 | argv[1]==72|argv[1] ==73 | argv[1]==74|argv[1] ==75 | argv[1]==76){
 ######## add in 73-76 - using data0 (only informative nuclear families)
#### Model 1.82 -- common beta_Age for all cohorts; N(mu_m,sigmasq_m) i.e. like 1.6.2 but no cohort effect - so we can compare to model 1.8 which has cohort-specific age effects but needs to avoid identifiability problem
    if(as.numeric(argv[1]) %%2 !=0){
#    cohort_by_family=rep(NA,data1.mat2$I)
#    for(i in 1:data1.mat2$I){
#        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
#    }
#    data1.mat2$cohort_by_family = cohort_by_family
    mean_alpha_prior=41    
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    sigmasq_mu_m_prior=100
    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    if(argv[1]==69){
       myname="model1.8.2.maternal"
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
       if(!print.only){       
           model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data1.mat2,iter=10000,chains=0)
           model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
       }
   }
    if(argv[1]==71){
        myname="model1.8.maternal"
        cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
        if(!print.only){
            model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data1.mat2,iter=10000,chains=0)
            model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0","beta_global","sigmasq_global")))
        }
        
    }

    if(!print.only){ 
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
        save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.",myname,".all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.6.mat,permuted=T)
    if("beta_global" %in% names(mysim)){
        print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
    } else {
        print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    }

} else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.",myname,".all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",                sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.6.mat,permuted=T)
}

} else {
 #    cohort_by_family=rep(NA,data1.pat2$I)
 #   for(i in 1:data1.pat2$I){
 #       cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
 #   }
  #  data1.pat2$cohort_by_family = cohort_by_family
     mean_alpha_prior=27
     sigmasq_m_alpha = 2
     sigmasq_m_beta = 15
     sigmasq_mu_m_prior=8^2
     
     data1.pat2$mean_alpha_prior = mean_alpha_prior
     data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
     data1.pat2$sigmasq_m_beta = sigmasq_m_beta
     data1.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior
     if(argv[1] ==70){
         myname="model1.8.2.paternal"
         cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
         if(!print.only){
             model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data1.pat2,iter=10000,chains=0)
             model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
         }
     }
     if(argv[1]==72){
         myname="model1.8.paternal"
         cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
         if(!print.only){
             model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data1.pat2,iter=10000,chains=0)
             model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
         }
     }
     if(!print.only){
         model1.6.pat <- sflist2stanfit(model1.6.pat.list)
         save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.",myname,".all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
         mysim<-extract(model1.6.pat,permuted=T)
    if("beta_global" %in% names(mysim)){
        print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
    } else {
        print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    }

     } else {
         load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/RSTAN.",myname,".all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
         mysim<-extract(model1.6.pat,permuted=T)
     }

 }
if(FALSE){
pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/",myname,".traceplots.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    traceplot(model1.6.mat,pars=c("beta_Age"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("mu_m"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("sigmasq_m"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("tausq"),window=c(10^2,10^4))
    if("beta_global" %in% names(mysim)){
        traceplot(model1.6.mat,pars=c("beta_global"),window=c(10^2,10^4))
        traceplot(model1.6.mat,pars=c("sigmasq_global"),window=c(10^2,10^4))
    }
    dev.off()
}
    
pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    my.ylim=100
    if(!"beta_global" %in% names(mysim)){
        plot(density(mysim$beta_Age),xlim=c(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
        abline(v=0,lwd=2)
        curve(dnorm(x,0,1),lty=2,lwd=2,add=T)
        legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2)
    } else {
        for(i in 1:ncol(mysim$beta_Age)){
            if(i==1){
                plot(density(mysim$beta_Age[,i]),xlim=c(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,10))
            }else {
                lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
            }
        }
        curve(dnorm(x,0,1),lty=2,lwd=2,add=T,col="black")
        abline(v=0,lwd=2)
        if("beta_global" %in% names(mysim)){
            lines(density(mysim$beta_global),lty=3,lwd=2,col="red")
            legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global","prior"),col=c(mycols[as.character(cohort.codes[,2])],"red","black"),lty=c(rep(1,nrow(cohort.codes)),3,2),lwd=2,cex=0.5)
        } else {
            legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=2,cex=0.5)
        }
    }

    if(as.numeric(argv[1]) %%2 !=0){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
          legend("topright",c("Campbell","Kong","Coop","Hussin","Bleazard"),lty=c(2,3,4,4,4),lwd=2,col=c("black","black","black","grey","darkgreen"),cex=0.5)
    }
    dev.off()
    
    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$mu_m),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",lwd=2)
    curve(dnorm(x,mean_alpha_prior,sigmasq_mu_m_prior),lwd=2,lty=2,add=T)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()


        pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$sigmasq_m),xlab="sigmasq_m",main="Posterior for sigmasq_m",lwd=2)
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS_more_stringent/no_GPC/",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),lwd=2,main="Posterior for tausq",xlab="tausq")
    lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
    dev.off()
}

