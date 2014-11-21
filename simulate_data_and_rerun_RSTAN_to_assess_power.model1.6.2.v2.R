library(rstan)
library(parallel)

description="model1.6.2.maternal.mu_m_N_41_100.sigmasq_m_IG_2_40"
#description="model1.6.2.paternal.mu_m_N_27_64.sigmasq_m_IG_2_15"
mydir="RSTAN_output_on_duoHMM_more_stringent/model1.62"
               #first prepare data for input
if(FALSE){
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_41_100.sigmasq_m_IG_2_40.RData")
#      load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_27_64.sigmasq_m_IG_2_15.RData")
    mysim<-extract(model1.6.mat,permuted=T)
 
    mybins=c(seq(from=0,to=0.10,by=0.02),max(mysim$beta_Age))
    n.in.bins=c()
    for(i in 1:(length(mybins)-1)){
        n.in.bins=c(n.in.bins,sum(mysim$beta_Age > mybins[i] & mysim$beta_Age <=mybins[i+1]))
    }

    mysim$bin.beta.age = 0
    mysim$bin.beta.age[mysim$beta_Age <=0] = -1
    samples.by.bin = list()
    for(i in 1:(length(mybins)-1)){
        mysim$bin.beta.age[mysim$beta_Age > mybins[i] & mysim$beta_Age <=mybins[i+1]] = mybins[i]

        indexes = which(mysim$beta_Age > mybins[i] & mysim$beta_Age <=mybins[i+1])
        samples = sample(indexes,1000,replace=T)
        samples.by.bin[[i]] = samples
    }

    save(mysim,file=paste0(mydir,"/",description,".draws_with_beta_age_binned.RData"))
    save(samples.by.bin,file=paste0(mydir,"/",description,".draws_to_sample_by_beta_age_bin.RData"))

}

load(file=paste0(mydir,"/",description,".draws_with_beta_age_binned.RData"))
load(file=paste0(mydir,"/",description,".draws_to_sample_by_beta_age_bin.RData"))
load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.more_stringent.RData")

setwd("/well/donnelly/hilary/maternal_age_and_recombination/")


argv <- commandArgs(trailingOnly = TRUE)
print(argv)

bin.number= as.numeric(argv[1])
sim.number =as.numeric(argv[2])
set.seed(sim.number+rnorm(1,100,2))

#need to use a coherent set of parameters
index=samples.by.bin[[bin.number]][sim.number]

family.effects=c()
####simulate data based on parameters from best draw
for(c in 1:length(unique(data1.mat2$cohort_by_family))){
    family.effects=c(family.effects,rnorm(n=sum(data1.mat2$cohort_by_family ==c),mean=mysim$mu_m[index,c],sd=sqrt(mysim$sigmasq_m[index,c])))
}

a0=family.effects[data1.mat2$family]
fake.y=c()

for(i in 1:length(a0)){
    fake.y=c(fake.y,rnorm(1,a0[i]+mysim$beta_Age[index] * data1.mat2$Age[i],sd=sqrt(mysim$tausq[index])))
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
save(beta.posterior,file=paste0(mydir,"/simulation_with_true_beta_Age_bin_",bin.number,".beta_Age_posteriors.",description,".simulation_",sim.number,".RData"))
write.table(myquantiles,paste0(mydir,"/simulation_with_true_beta_Age_bin_",bin.number,".beta_Age_credible_intervals.",description,".simulation_",sim.number,".txt"))
}

if(FALSE){
coverage.by.beta.age=list()
exclude=list()
load(file=paste0(mydir,"/",description,".draws_with_beta_age_binned.RData"))
load(file=paste0(mydir,"/",description,".draws_to_sample_by_beta_age_bin.RData"))

coverage.by.bin=NULL
for(b in 1:6){
    all.coverage=as.data.frame(matrix(NA,ncol=100,nrow=8))
    for(sim.number in 1:100){
        load(file=paste0(mydir,"/simulation_with_true_beta_Age_bin_",b,".beta_Age_posteriors.",description,".simulation_",sim.number,".RData"))
        myquantiles=quantile(beta.posterior,c(1/6,5/6))
        intervals=read.delim(paste0(mydir,"/simulation_with_true_beta_Age_bin_",b,".beta_Age_credible_intervals.",description,".simulation_",sim.number,".txt"),header=T,sep="")
        beta_Age = mysim$beta_Age[samples.by.bin[[b]][sim.number]]
        coverage=c(            intervals[1,1]<beta_Age & intervals[2,1] > beta_Age, intervals[1,1]<0 & intervals[2,1] > 0,
            myquantiles[1] <beta_Age & myquantiles[2] > beta_Age,myquantiles[1] <0 & myquantiles[2] > 0,
            intervals[3,1]<beta_Age & intervals[4,1] > beta_Age,intervals[3,1]<0 & intervals[4,1] > 0,
            intervals[5,1]<beta_Age & intervals[6,1] > beta_Age, intervals[5,1]<0 & intervals[6,1] > 0)
        all.coverage[,sim.number]=coverage
    }
                                        #    exclude[[i]]=which(colSums(is.na(all.intervals))!=0)
                                        #    all.intervals=all.intervals[,colSums(is.na(all.intervals))==0]
 #   coverage.by.beta.age[[i]]=all.intervals
#    coverage=c(sum(all.intervals[1,]<beta_Age & all.intervals[2,] > beta_Age), sum(all.intervals[1,]<0 & all.intervals[2,] > 0),sum(all.intervals[3,]<beta_Age & all.intervals[4,] > beta_Age),
#        sum(all.intervals[3,]<0 & all.intervals[4,] > 0),sum(all.intervals[5,]<beta_Age & all.intervals[6,] > beta_Age), sum(all.intervals[5,]<0 & all.intervals[6,] > 0))/ncol(all.intervals)
    all.coverage=all.coverage[,colSums(is.na(all.coverage))==0]

    total.coverage = rowSums(all.coverage)/ncol(all.coverage)
    
    names(total.coverage)=c("coverage_beta_Age_33%","coverage_0_33%","coverage_beta_Age_67%","coverage_0_67%","coverage_beta_Age_95%","coverage_0_95%","coverage_beta_Age_99%","coverage_0_99%")
    coverage.by.bin = rbind(coverage.by.bin,total.coverage)
#    write.table(coverage,paste0("/home/hilary/maternal_age_recombination/frequentist_properties/coverage_probabilities_of_beta_Age.simulated_with_beta_Age_",beta_Age,".",description,".txt"),quote=F,sep="\t")
}
write.table(coverage.by.bin,paste0(mydir,"/coverage_and_power_by_bin_for_beta_Age.",description,".txt"),quote=F,sep="\t")
            
mybins=c(seq(from=0,to=0.10,by=0.02),max(mysim$beta_Age))
   
#pdf("/home/hilary/maternal_age_recombination/frequentist_properties/coverage_plot_model1.6.2.maternal.pdf",height=4,width=5,useDingbats=FALSE)
pdf("coverage_plot_model1.6.2.maternal.pdf",height=4,width=5,useDingbats=FALSE)
myx=mybins[1:6]+(mybins[2:7]-mybins[1:6])/2
plot(myx,coverage.by.bin[,"coverage_beta_Age_67%"]*100,col="blue",xlab="center of beta_Age bin",ylab="Coverage (%)",ylim=c(0,100),xlim=c(0,max(mysim$beta_Age)),type="b",pch=20)
for(i in 1:6){
#segments(x0=mybins[i],y0=coverage.by.bin[i,"coverage_beta_Age_67%"]*100,x1=mybins[i+1],y1=coverage.by.bin[i,"coverage_beta_Age_67%"]*100,col="blue",lwd=2)
#segments(x0=mybins[i],x1=mybins[i],y0=coverage.by.bin[i,"coverage_beta_Age_67%"]*100-2,y1=coverage.by.bin[i,"coverage_beta_Age_67%"]*100+2,col="blue")
#segments(x0=mybins[i+1],x1=mybins[i+1],y0=coverage.by.bin[i,"coverage_beta_Age_67%"]*100-2,y1=coverage.by.bin[i,"coverage_beta_Age_67%"]*100+2,col="blue")
}

  lines(myx,coverage.by.bin[,"coverage_beta_Age_33%"]*100,type="b",col="darkgreen",pch=20)

  lines(myx,coverage.by.bin[,"coverage_beta_Age_95%"]*100,type="b",col="purple",pch=20)

for(i in 1:6){
#segments(x0=mybins[i],y0=coverage.by.bin[i,"coverage_beta_Age_95%"]*100,x1=mybins[i+1],y1=coverage.by.bin[i,"coverage_beta_Age_95%"]*100,col="purple",lwd=2)
#segments(x0=mybins[i],x1=mybins[i],y0=coverage.by.bin[i,"coverage_beta_Age_95%"]*100-2,y1=coverage.by.bin[i,"coverage_beta_Age_95%"]*100+2,col="purple")
#segments(x0=mybins[i+1],x1=mybins[i+1],y0=coverage.by.bin[i,"coverage_beta_Age_95%"]*100-2,y1=coverage.by.bin[i,"coverage_beta_Age_95%"]*100+2,col="purple")
}
lines(myx,coverage.by.bin[,"coverage_beta_Age_99%"]*100,type="b",col="red",pch=20)
for(i in 1:6){
#segments(x0=mybins[i],y0=coverage.by.bin[i,"coverage_beta_Age_99%"]*100,x1=mybins[i+1],y1=coverage.by.bin[i,"coverage_beta_Age_99%"]*100,col="red",lwd=2)
#segments(x0=mybins[i],x1=mybins[i],y0=coverage.by.bin[i,"coverage_beta_Age_99%"]*100-2,y1=coverage.by.bin[i,"coverage_beta_Age_99%"]*100+2,col="red")
#segments(x0=mybins[i+1],x1=mybins[i+1],y0=coverage.by.bin[i,"coverage_beta_Age_99%"]*100-2,y1=coverage.by.bin[i,"coverage_beta_Age_99%"]*100+2,col="red")
}
legend("right",c("33% CI","67% CI","95% CI","99% CI"),col=c("darkgreen","blue","purple","red"),lty=1,cex=0.6,pch=20)
dev.off()

#pdf("/home/hilary/maternal_age_recombination/frequentist_properties/power_plot_model1.6.2.maternal.pdf",height=4,width=5,useDingbats=FALSE)
pdf("power_plot_model1.6.2.maternal.pdf",height=4,width=5,useDingbats=FALSE)
myx=mybins[1:6]+(mybins[2:7]-mybins[1:6])/2
plot(myx,100-coverage.by.bin[,"coverage_0_67%"]*100,col="blue",xlab="center for beta_Age bin",ylab="Power (%)",ylim=c(0,100),xlim=c(0,max(mysim$beta_Age)),type="b",pch=20)
for(i in 1:6){
#segments(x0=mybins[i],y0=100-coverage.by.bin[i,"coverage_0_67%"]*100,x1=mybins[i+1],y1=100-coverage.by.bin[i,"coverage_0_67%"]*100,col="blue",lwd=2)
}
lines(myx,100-coverage.by.bin[,"coverage_0_33%"]*100,type="b",col="darkgreen",pch=20)
lines(myx,100-coverage.by.bin[,"coverage_0_95%"]*100,col="purple",type="b",pch=20)
for(i in 1:6){
#segments(x0=mybins[i],y0=100-coverage.by.bin[i,"coverage_0_95%"]*100,x1=mybins[i+1],y1=100-coverage.by.bin[i,"coverage_0_95%"]*100,col="purple",lwd=2)
}
lines(myx,100-coverage.by.bin[,"coverage_0_99%"]*100,col="red",type="b",pch=20)
for(i in 1:6){
#segments(x0=mybins[i],y0=100-coverage.by.bin[i,"coverage_0_99%"]*100,x1=mybins[i+1],y1=100-coverage.by.bin[i,"coverage_0_99%"]*100,col="red",lwd=2)
}
legend("topleft",c("33% CI","67% CI","95% CI","99% CI"),col=c("darkgreen","blue","purple","red"),lty=1,cex=0.6,pch=20)
dev.off()



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
