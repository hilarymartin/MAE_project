library(rstan)
library(parallel)

#description="model1.6.2.maternal.mu_m_N_41_100.sigmasq_m_IG_2_40"
description="model1.6.2.paternal.mu_m_N_27_64.sigmasq_m_IG_2_15"
mydir="RSTAN_output_on_duoHMM_more_stringent/model1.62"
               #first prepare data for input
if(FALSE){
#    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_41_100.sigmasq_m_IG_2_40.RData")
      load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_27_64.sigmasq_m_IG_2_15.RData")
#    mysim<-extract(model1.6.mat,permuted=T)
    mysim<-extract(model1.6.pat,permuted=T)
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp_)])
        } else {
            return(x[which.max(mysim$lp_),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
    save(best.draw,file=paste0(mydir,"/",description,".best_draw.RData"))


argv <- commandArgs(trailingOnly = TRUE)
print(argv)

setwd("/well/donnelly/hilary/maternal_age_and_recombination/")
#load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.NTR_v2.RData")
load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.more_stringent.RData")
#load("RSTAN_output_with_NTR_v2/model1.62/model1.6.2.maternal.mu_m_N_41_100.sigmasq_m_IG_2_40.best_draw.RData")
load(paste0(mydir,"/",description,".best_draw.RData"))

beta_Age = as.numeric(argv[1])
#beta_Age = 0.06
sim.number =as.numeric(argv[2])
#sim.number =1
set.seed(sim.number)

family.effects=c()
####simulate data based on parameters from best draw
for(c in 1:length(unique(data1.mat2$cohort_by_family))){
    family.effects=c(family.effects,rnorm(n=sum(data1.mat2$cohort_by_family ==c),mean=best.draw$mu_m[c],sd=sqrt(best.draw$sigmasq_m[c])))
}

a0=family.effects[data1.mat2$family]
fake.y=c()

for(i in 1:length(a0)){
fake.y=c(fake.y,rnorm(1,a0[i]+beta_Age * data1.mat2$Age[i],sd=sqrt(best.draw$tausq)))
}
summary(fake.y)
###prepare data
mean_alpha_prior=41
sigmasq_m_alpha = 2
sigmasq_m_beta = 40
sigmasq_mu_m_prior=100

data1.mat2$mean_alpha_prior = mean_alpha_prior
data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
data1.mat2$sigmasq_m_beta = sigmasq_m_beta
data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior

sim.data=data1.mat2
sim.data$y=fake.y

### run model
model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.stan",data=sim.data,iter=10000,chains=0)
model1.6.mat = stan(fit = model1.6.mat.compiled, data = sim.data,chains = 4, iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))

print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
mysim<-extract(model1.6.mat,permuted=T)
myquantiles=quantile(mysim$beta_Age,c(1/3,2/3,0.025,0.975,0.005,0.995))
beta.posterior=mysim$beta_Age
save( beta.posterior,file=paste0(mydir,"/simulation_with_true_beta_Age_",beta_Age,".beta_Age_posteriors.",description,"simulation_",sim.number,".RData"))
write.table(myquantiles,paste0(mydir,"/simulation_with_true_beta_Age_",beta_Age,".beta_Age_credible_intervals.",description,".simulation_",sim.number,".txt"))
}
if(FALSE){
coverage.by.beta.age=list()
exclude=list()
for(i in 1:3){
    beta_Age = c(0,0.067,0.082)[i]
    all.intervals=as.data.frame(matrix(NA,ncol=1000,nrow=6))
    for(sim.number in 1:1000){
        intervals=read.delim(paste0(mydir,"/simulation_with_true_beta_Age_",beta_Age,".beta_Age_credible_intervals.",description,".simulation_",sim.number,".txt"),header=T,sep="")
        all.intervals[,sim.number]=intervals
    }
    exclude[[i]]=which(colSums(is.na(all.intervals))!=0)
    all.intervals=all.intervals[,colSums(is.na(all.intervals))==0]
    coverage.by.beta.age[[i]]=all.intervals
    coverage=c(sum(all.intervals[1,]<beta_Age & all.intervals[2,] > beta_Age), sum(all.intervals[1,]<0 & all.intervals[2,] > 0),sum(all.intervals[3,]<beta_Age & all.intervals[4,] > beta_Age),
        sum(all.intervals[3,]<0 & all.intervals[4,] > 0),sum(all.intervals[5,]<beta_Age & all.intervals[6,] > beta_Age), sum(all.intervals[5,]<0 & all.intervals[6,] > 0))/ncol(all.intervals)
    names(coverage)=c("coverage_beta_Age_67%","coverage_0_67%","coverage_beta_Age_95%","coverage_0_95%","coverage_beta_Age_99%","coverage_0_99%")
    write.table(coverage,paste0("/home/hilary/maternal_age_recombination/frequentist_properties/coverage_probabilities_of_beta_Age.simulated_with_beta_Age_",beta_Age,".",description,".txt"),quote=F,sep="\t")
}

for(i in 1:3){
    beta_Age = c(0,0.067,0.082)[i]
        beta.posteriors=NULL
    for(sim.number in 1:1000){
                if(!sim.number %in% exclude[[i]]){
                    load(paste0(mydir,"/simulation_with_true_beta_Age_",beta_Age,".beta_Age_posteriors.",description,"simulation_",sim.number,".RData"))
                    beta.posteriors=cbind(beta.posteriors,beta.posterior)
                    rm(beta.posterior)
                }
            }
    pdf(paste0("/home/hilary/maternal_age_recombination/frequentist_properties/distribution_of_beta_Age_posteriors.simulated_with_beta_Age_",beta_Age,".",description,".pdf"),height=5,width=5)
    for(sim.number in 1:ncol(beta.posteriors)){
        if(sim.number==1){
            plot(density(beta.posteriors[,sim.number]),xlab="beta_Age",main="posteriors of beta_Age",ylim=c(0,12))
        } else {
            lines(density(beta.posteriors[,sim.number]))
        }
        abline(v=0,col="blue",lwd=3)
        abline(v=beta_Age,col="red",lwd=3)
    }
    dev.off()
}
}
