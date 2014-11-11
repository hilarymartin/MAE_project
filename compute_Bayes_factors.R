if(FALSE){
library(rstan)
load("/well/donnelly/hilary/TwinRecombination/Bayesian_modeling/RData_files_all_cohorts/normal.fit.cohort.effect.RData")
load("/well/donnelly/hilary/TwinRecombination/Bayesian_modeling/RData_files_all_cohorts/normal.fit.cohort.specific.b.RData")

#load("/home/hilary/TwinRecombination/Bayesian_modeling/input_data/NFTOOLS_counts_for_all_cohorts.RData")
load("/well/donnelly/hilary/TwinRecombination/Bayesian_modeling/input_data/NFTOOLS_counts_for_all_cohorts.RData")
print(normal.model.cohort.effect,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
model1.sim<-extract(normal.model.cohort.effect,permuted=F)
model2.sim<-extract(normal.fit.cohort.specific.b,permuted=F)

model1.sim.2=rbind(model1.sim[,1,],model1.sim[,2,],model1.sim[,3,],model1.sim[,4,])
model2.sim.2=rbind(model2.sim[,1,],model2.sim[,2,],model2.sim[,3,],model2.sim[,4,])

model1.keep=model1.sim.2[ which((1:nrow(model1.sim.2)) %%200 ==0),]
model2.keep=model2.sim.2[ which((1:nrow(model2.sim.2)) %%200 ==0),]

#look autocorrelation
pdf("/home/hilary/TwinRecombination/Bayesian_modeling/results_from_4_cohorts/normal_model_cohort_effect/autocorrelation.beta_Age.pdf",height=6,width=6)
acf( model1.sim[,1,"beta_Age"], lag.max=200)
acf( model1.sim[,1,"tausq"], lag.max=200)
acf( model1.sim[,1,"beta_Cohort[1]"], lag.max=200)
dev.off()
#lots of autocorrelation - need to pick every 1/200 to avoid this, so will have only 100 samples --> could run chains for longer
}



model1.probs=(1/sqrt(2*pi*model1.keep[,"tausq"]))*exp(-0.5*(t(matrix(rep(twin_data$y,100),ncol=100,nrow=length(twin_data$y)))-(model1.keep[,grep("a0_bychild",colnames(model1.keep))]
        +c(model1.keep[,"beta_Age"])%*% t(twin_data$Age)+model1.keep[,grep("my_beta_Cohort",colnames(model1.keep))]))^2/   model1.keep[,"tausq"])
model1.log.likelihoods=rowSums(log(model1.probs))
            
model2.probs=(1/sqrt(2*pi*model2.keep[,"tausq"]))*exp(-0.5*(t(matrix(rep(twin_data$y,100),ncol=100,nrow=length(twin_data$y)))-(model2.keep[,grep("a0_bychild",colnames(model2.keep))]+
               model2.keep[,grep("beta_Age",colnames(model2.keep))][,twin_data$cohort+1] %*% diag(twin_data$Age)+model2.keep[,grep("my_beta_Cohort",colnames(model2.keep))]))^2/   model2.keep[,"tausq"])
model2.log.likelihoods=rowSums(log(model2.probs))

          
compare.log.likes=rbind(model1.log.likelihoods,model2.log.likelihoods)

### use RYan's matrixlogadd function to do the sums in log space, since probabilities are too small to compute
### this is calculating sum_over_draws(exp(sum_over_individuals(log(model.probs))))
log.likes =matrixlogadd(compare.log.likes)

BF.1.vs.2 = exp(log.likes [1]-log.likes [2])
BF.2.vs.1= exp(log.likes [2]-log.likes [1])

#> BF.1.vs.2
#            0.03700851
#> BF.2.vs.1
#              27.02081 
##feed matrixlogadd a matrix with entries in row 1 being sum_over_n_indivs(log(P(y_i|theta_Model1_sample_s))),row 2 sum_over_n_indivs(log(P(y_i|theta_Model2_sample_s)))
            
matrixlogadd<-function(v){
      #Does Addition in Log Space down each row to get a sum for each row
      N<-dim(v)[2]
        y<-v[max.col(v)]
        return(y+log(rowSums(exp(v-matrix(rep(y,N),nrow=length(y),N)))))
  }
