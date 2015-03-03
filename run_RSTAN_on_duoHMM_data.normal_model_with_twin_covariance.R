library(rstan)
library(scales)
library(parallel)
library(pscl)
load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.more_stringent.twin_cohorts_inc_family_matrix.RData")
argv <- commandArgs(trailingOnly = TRUE)
print(argv)    
mycols=c("black","blue","red","green","orange","red4","purple","darkgreen","violetred2")
names(mycols)=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB","ORCADES")

print.only=as.logical(as.numeric(argv[2]))

#### Model 3.

if(argv[1] ==1 | argv[1]==2){
#### Model 1.62 -- common beta_Age for all cohorts; N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort

    if(argv[1]==1){
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
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.normal_model_with_twin_covariance.stan",data=data1.mat2,iter=10000,chains=0)
#model1.6.mat =stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0","lambda"))
#model1.6.mat =stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.normal_model_with_twin_covariance.stan", data = data1.mat2,chains = 1, iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0","lambda"))
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0","lambda")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m","lambda"),digits=4)
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.normal_model_with_twin_covariance.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.normal_model_with_twin_covariance.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
#    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.txt",header=T)
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
     sigmasq_mu_m_prior=8^2
     data1.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior
     if(!print.only){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.normal_model_with_twin_covariance.stan",data=data1.pat2,iter=10000,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0","lambda")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m","lambda"),digits=4)
     save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.normal_model_with_twin_covariance.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.normal_model_with_twin_covariance.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
 }
     mysim<-extract(model1.6.pat,permuted=T)
 #   cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.txt",header=T)
    myname="paternal"
}

if(FALSE){
pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/model1.6.2.",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    my.ylim=100
#    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
        plot(density(mysim$beta_Age),xlim=c(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
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
    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/model1.6.2.",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
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

    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/model1.6.2.",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.15
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
        if(argv[1]==1){
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"maternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        } else {
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"paternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        }
    }
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior","empirical variance in parent means"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,2),lwd=c(rep(2,nrow(cohort.codes)),3,1),cex=0.5)
    dev.off()

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/model1.6.2.",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),lwd=2,main="Posterior for tausq",xlab="tausq")
    lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
    dev.off()
}
}

