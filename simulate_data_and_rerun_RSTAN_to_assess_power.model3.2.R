library(rstan)
library(parallel)

description="model3.2.maternal.mu_m_N_36_6.sigmasq_m_IG_5_5"
mydir="RSTAN_output_on_duoHMM_more_stringent/model3.2"

#first prepare data for input
if(FALSE){
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.model3.2.adjusted_priors.all_chains.mat.including_alpha.RData")
    mysim<-extract(model3.2.mat,permuted=T)
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp_)])
        } else {
            return(x[which.max(mysim$lp_),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
    best.draw[["inv_omega"]]=1/best.draw[["omega"]]
    save(best.draw,file=paste0(mydir,"/",description,".best_draw.RData"))
}

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


####simulate data based on parameters from best draw

family.effects=rnorm(n=length(unique(data2.mat2$family)),mean=best.draw$mu_m,sd=sqrt(best.draw$sigmasq_m))

a0=family.effects[data2.mat2$family]

#negative binomial
fake.y=c()

#        return(rnbinom(1,size=exp(0.1*(a0[i]+beta_Age*age[i]))/(omega-1),prob=(1/(p_by_ cohort[cohort[i]] * (omega-1) + 1))))

for(i in 1:length(a0)){
    fake.y=c(fake.y,rnbinom(1,size=exp(0.1*(a0[i] + beta_Age*data2.mat2$Age[i]))/(best.draw$omega-1),prob=(1/(best.draw$p_by_cohort[data2.mat2$cohort[i]] * (best.draw$omega-1) + 1))))
}
summary(fake.y)
###prepare data
sim.data=data2.mat2
sim.data$y=fake.y

### run model
model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan",data=sim.data,iter=10000,chains=0)
model1.6.mat = stan(fit = model1.6.mat.compiled, data = sim.data,chains = 4, iter=10000,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0"))

print(model1.6.mat,pars=c("beta_Age","omega","sigmasq_m","mu_m"),digits=4)
mysim<-extract(model1.6.mat,permuted=T)
myquantiles=quantile(mysim$beta_Age,c(1/3,2/3,0.025,0.975,0.005,0.995))
beta.posterior=mysim$beta_Age
save( beta.posterior,file=paste0(mydir,"/simulation_with_true_beta_Age_",beta_Age,".beta_Age_posteriors.",description,"simulation_",sim.number,".RData"))
write.table(myquantiles,paste0(mydir,"/simulation_with_true_beta_Age_",beta_Age,".beta_Age_credible_intervals.",description,".simulation_",sim.number,".txt"))

if(FALSE){
coverage.by.beta.age=list()
for(i in 1:2){
#    beta_Age = c(0.067,0.19)[i]
    beta_Age = c(0,0.082)[i]
    all.intervals=as.data.frame(matrix(NA,ncol=1000,nrow=6))
#    for(sim.number in 1:1000){
x=0
    for(sim.number in 1001:2000){
        if(sim.number!=1450){
        x=x+1
        intervals=read.delim(paste0(mydir,"/simulation_with_true_beta_Age_",beta_Age,".beta_Age_credible_intervals.",description,".simulation_",sim.number,".txt"),header=T,sep="")
        all.intervals[,x]=intervals
    }
    }
    all.intervals=all.intervals[,colSums(is.na(all.intervals))==0]
    coverage.by.beta.age[[i]]=all.intervals
    coverage=c(sum(all.intervals[1,]<beta_Age & all.intervals[2,] > beta_Age), sum(all.intervals[1,]<0 & all.intervals[2,] > 0),sum(all.intervals[3,]<beta_Age & all.intervals[4,] > beta_Age),
        sum(all.intervals[3,]<0 & all.intervals[4,] > 0),sum(all.intervals[5,]<beta_Age & all.intervals[6,] > beta_Age), sum(all.intervals[5,]<0 & all.intervals[6,] > 0))/ncol(all.intervals)
    names(coverage)=c("coverage_beta_Age_67%","coverage_0_67%","coverage_beta_Age_95%","coverage_0_95%","coverage_beta_Age_99%","coverage_0_99%")
    write.table(coverage,paste0("/home/hilary/maternal_age_recombination/coverage_probabilities_of_beta_Age.simulated_with_beta_Age_",beta_Age,
".",description,".txt"),quote=F,sep="\t")
}

for(i in 1:2){
    beta_Age = c(0,0.082)[i]
    beta.posteriors=NULL
#    for(sim.number in 1:1000){
    for(sim.number in 1001:2000){
                if(sim.number!=1450){
                    load(paste0(mydir,"/simulation_with_true_beta_Age_",beta_Age,
                                ".beta_Age_posteriors.",description,".simulation_",sim.number,".RData"))
        beta.posteriors=cbind(beta.posteriors,beta.posterior)
                    rm(beta.posterior)
                }
            }
    pdf(paste0("/home/hilary/maternal_age_recombination/distribution_of_beta_Age_posteriors.simulated_with_beta_Age_",beta_Age,
               ".",description,".pdf"),height=5,width=5)
#    for(sim.number in 1:1000){
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
