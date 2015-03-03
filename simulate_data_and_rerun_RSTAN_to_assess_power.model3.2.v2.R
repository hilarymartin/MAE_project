library(rstan)
library(parallel)

#description="model3.2.maternal.mu_m_N_36_6.sigmasq_m_IG_5_5"
description="model3.2.adjusted_priors.different_link"
mydir="RSTAN_output_on_duoHMM_more_stringent/with_new_QTR.no_GPC_different_link_for_NB_changing_sigmasq_m_prior/model3.2"

#first prepare data for input

#    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.model3.2.adjusted_priors.all_chains.mat.including_alpha.RData")
#        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/longer_chains/RSTAN.model3.2.adjusted_priors.all_chains.mat.including_alpha.RData")
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/with_new_QTR.no_GPC_different_link_for_NB_changing_sigmasq_m_prior/RSTAN.model3.2.adjusted_priors.different_link.all_chains.mat.including_alpha.RData")
    mysim<-extract(model3.2.mat,permuted=T)

        # -0.01698482  0.04226615
#    mybins=c(seq(from=0,to=0.025,by=0.005),max(mysim$beta_Age))
#            mybins=c(-0.01698482,seq(from=-0.002,to=0.03,by=0.004),0.04226615)
            mybins=c(min(mysim$beta_Age),seq(from=0.0011,to=0.0037,by=0.0003),max(mysim$beta_Age))
    n.in.bins=c()
    for(i in 1:(length(mybins)-1)){
        n.in.bins=c(n.in.bins,sum(mysim$beta_Age > mybins[i] & mysim$beta_Age <=mybins[i+1]))
    }
if(FALSE){
    mysim$bin.beta.age = 0
    mysim$bin.beta.age[mysim$beta_Age <=0] = -1
    samples.by.bin = list()
    for(i in 1:(length(mybins)-1)){
        mysim$bin.beta.age[mysim$beta_Age > mybins[i] & mysim$beta_Age <=mybins[i+1]] = mybins[i]
        indexes = which(mysim$beta_Age > mybins[i] & mysim$beta_Age <=mybins[i+1])
#        samples = sample(indexes,1000,replace=T)
        samples = sample(indexes,100,replace=T)
        samples.by.bin[[i]] = samples
    }
    
    save(mysim,file=paste0(mydir,"/",description,".draws_with_beta_age_binned.RData"))
#            save(mysim,file=paste0(mydir,"/",description,".longer_chains.draws_with_beta_age_binned.RData"))

    save(samples.by.bin,file=paste0(mydir,"/",description,".draws_to_sample_by_beta_age_bin.RData"))

    
}
if(FALSE){
argv <- commandArgs(trailingOnly = TRUE)
print(argv)
load(file=paste0(mydir,"/",description,".draws_with_beta_age_binned.RData"))
#load(file=paste0(mydir,"/",description,".longer_chains.draws_with_beta_age_binned.RData"))
load(paste0(mydir,"/",description,".draws_to_sample_by_beta_age_bin.RData"))
#load(file=paste0(mydir,"/",description,".longer_chains.draws_to_sample_by_beta_age_bin.RData"))

setwd("/well/donnelly/hilary/maternal_age_and_recombination/")
#load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.more_stringent.RData")
load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.more_stringent.with_new_QTR.no_GPC.RData")
#load(paste0(mydir,"/",description,".best_draw.RData"))

bin.number= as.numeric(argv[1])
sim.number =as.numeric(argv[2])
set.seed(sim.number)

#need to use a coherent set of parameters
index=samples.by.bin[[bin.number]][sim.number]

####simulate data based on parameters 
family.effects=rnorm(n=length(unique(data2.mat2$family)),mean=mysim$mu_m[index],sd=sqrt(mysim$sigmasq_m[index]))

a0=family.effects[data2.mat2$family]
#negative binomial
fake.y=c()
for(i in 1:length(a0)){
#    fake.y=c(fake.y,rnbinom(1,size=exp(0.1*(a0[i] + mysim$beta_Age[index]*data2.mat2$Age[i]))/(mysim$omega[index]-1),prob=(1/(mysim$p_by_cohort[index,data2.mat2$cohort[i]] * (mysim$omega[index]-1) + 1))))
    fake.y=c(fake.y,rnbinom(1,size=exp((a0[i] + mysim$beta_Age[index]*data2.mat2$Age[i]))/(mysim$omega[index]-1),prob=(1/(mysim$p_by_cohort[index,data2.mat2$cohort[i]] * (mysim$omega[index]-1) + 1))))
}
summary(fake.y)
###prepare data
sim.data=data2.mat2
sim.data$y=fake.y

#sim.data$mu_alpha_prior = 37
sim.data$mu_alpha_prior = 3.7
#sim.data$sigmasq_alpha_prior = 6
sim.data$sigmasq_alpha_prior = 0.2


### run model
#model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan",data=sim.data,iter=10000,chains=0)
model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.different_link.stan",data=sim.data,iter=10000,chains=0)
model1.6.mat = stan(fit = model1.6.mat.compiled, data = sim.data,chains = 4, iter=10000,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0"))
print(model1.6.mat,pars=c("beta_Age","omega","sigmasq_m","mu_m"),digits=4)
mysim<-extract(model1.6.mat,permuted=T)
myquantiles=quantile(mysim$beta_Age,c(1/3,2/3,0.25,0.75,1/6,5/6,0.025,0.975,0.005,0.995))
beta.posterior=mysim$beta_Age
save( beta.posterior,file=paste0(mydir,"/simulation_with_true_beta_Age_bin_",bin.number,".beta_Age_posteriors.",description,".simulation_",sim.number,".RData"))
#save(beta.posterior,file=paste0(mydir,"/longer_chains/simulation_with_true_beta_Age_bin_",bin.number,".beta_Age_posteriors.",description,".simulation_",sim.number,".RData"))
write.table(myquantiles,paste0(mydir,"/simulation_with_true_beta_Age_bin_",bin.number,".beta_Age_credible_intervals.",description,".simulation_",sim.number,".txt"))

}

if(FALSE){
coverage.by.beta.age=list()
exclude=list()

#load(file=paste0(mydir,"/",description,".longer_chains.draws_to_sample_by_beta_age_bin.RData"))
load(file=paste0(mydir,"/",description,".draws_to_sample_by_beta_age_bin.RData"))
#load(file=paste0(mydir,"/",description,".longer_chains.draws_with_beta_age_binned.RData"))
load(file=paste0(mydir,"/",description,".draws_with_beta_age_binned.RData"))

coverage.by.bin=NULL
for(b in 1:10){
    all.coverage=as.data.frame(matrix(NA,ncol=100,nrow=10))
    for(sim.number in 1:100){
        file.exists=FALSE
#        if(file.exists(paste0(mydir,"/longer_chains/simulation_with_true_beta_Age_bin_",b,".beta_Age_posteriors.",description,".simulation_",sim.number,".RData"))){
                if(file.exists(paste0(mydir,"/simulation_with_true_beta_Age_bin_",b,".beta_Age_posteriors.",description,".simulation_",sim.number,".RData"))){
          file.exists=TRUE
#          intervals=read.delim(paste0(mydir,"/longer_chains/simulation_with_true_beta_Age_bin_",b,".beta_Age_credible_intervals.",description,".simulation_",sim.number,".txt"),header=T,sep="")
                    intervals=read.delim(paste0(mydir,"/simulation_with_true_beta_Age_bin_",b,".beta_Age_credible_intervals.",description,".simulation_",sim.number,".txt"),header=T,sep="")
            }
        if(file.exists){
            beta_Age = mysim$beta_Age[samples.by.bin[[b]][sim.number]]
            coverage=c(            intervals[1,1]<beta_Age & intervals[2,1] > beta_Age, intervals[1,1]<0 & intervals[2,1] > 0,
              intervals[3,1]<beta_Age & intervals[4,1] > beta_Age,intervals[3,1]<0 & intervals[4,1] > 0,
              intervals[5,1]<beta_Age & intervals[6,1] > beta_Age, intervals[5,1]<0 & intervals[6,1] > 0,
              intervals[7,1]<beta_Age & intervals[8,1] > beta_Age, intervals[7,1]<0 & intervals[8,1] > 0,
              intervals[9,1]<beta_Age & intervals[10,1] > beta_Age, intervals[9,1]<0 & intervals[10,1] > 0)
             
            all.coverage[,sim.number]=coverage
        }
    }
    all.coverage=all.coverage[,colSums(is.na(all.coverage))==0]

    total.coverage = rowSums(all.coverage)/ncol(all.coverage)
    myquantiles=c("33%","50%","66%","95%","99%")
    names(total.coverage)= c(sapply(myquantiles,function(x){paste(c("coverage_beta_age_","coverage_0_"),x,sep="")}))
      coverage.by.bin = rbind(coverage.by.bin,total.coverage)
}
write.table(coverage.by.bin,paste0(mydir,"/coverage_and_power_by_bin_for_beta_Age.",description,".txt"),quote=F,sep="\t")

mycolors=c("orange","blue","darkgreen","purple","red")

#pdf(paste0(mydir,"/coverage_plot_model3.2.maternal.pdf"),height=4,width=5,useDingbats=FALSE)
pdf(paste0(mydir,"/coverage_plot_model3.2.maternal.exponentiated.pdf"),height=4,width=5,useDingbats=FALSE)
#myx=mybins[1:(length(mybins)-1)]+(mybins[2:length(mybins)]-mybins[1:(length(mybins)-1)])/2
myx=exp(mybins[1:(length(mybins)-1)]+(mybins[2:length(mybins)]-mybins[1:(length(mybins)-1)])/2)
for(i in 1:length(myquantiles)){
  quant=myquantiles[i]
  if(i==1){
#    plot(myx,coverage.by.bin[,paste0("coverage_beta_age_",quant)]*100,col=mycolors[i],xlab="center of beta_Age bin",ylab="Coverage (%)",ylim=c(0,100),xlim=range(mybins),type="b",pch=20,
          plot(myx,coverage.by.bin[,paste0("coverage_beta_age_",quant)]*100,col=mycolors[i],xlab="center of beta_Age bin",ylab="Coverage (%)",ylim=c(0,100),xlim=range(exp(mybins)),type="b",pch=20,
         main="Model 3")
} else {
  lines(myx,coverage.by.bin[,paste0("coverage_beta_age_",quant)]*100,col=mycolors[i],type="b",pch=20)
}
}
abline(v=0)
abline(h=80,lty=2)
  legend("left",paste(myquantiles,"CI"),col=mycolors,lty=1,cex=0.6,pch=20)
dev.off()


#pdf(paste0(mydir,"/power_plot_model3.2.maternal.pdf"),height=4,width=5,useDingbats=FALSE)
pdf(paste0(mydir,"/power_plot_model3.2.maternal.exponentiated.pdf"),height=4,width=5,useDingbats=FALSE)
myx=exp(mybins[1:(length(mybins)-1)]+(mybins[2:length(mybins)]-mybins[1:(length(mybins)-1)])/2)
for(i in 1:length(myquantiles)){
  quant=myquantiles[i]
  if(i==1){
#    plot(myx,100-coverage.by.bin[,paste0("coverage_0_",quant)]*100,col=mycolors[i],xlab="center of beta_Age bin",ylab="Power (%)",ylim=c(0,100),xlim=range(mybins),type="b",pch=20,main="Model 3")
          plot(myx,100-coverage.by.bin[,paste0("coverage_0_",quant)]*100,col=mycolors[i],xlab="center of beta_Age bin",ylab="Power (%)",ylim=c(0,100),xlim=range(exp(mybins)),type="b",pch=20,main="Model 3")
} else {
  lines(myx,100-coverage.by.bin[,paste0("coverage_0_",quant)]*100,col=mycolors[i],type="b",pch=20)
}
}
abline(v=0)
abline(h=80,lty=2)

legend("left",paste(myquantiles,"CI"),col=mycolors,lty=1,cex=0.6,pch=20)

dev.off()

}
