library(rstan)
library(scales)
library(parallel)
load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.NTR_v2.RData")
argv <- commandArgs(trailingOnly = TRUE)
print(argv)
setwd("/well/donnelly/hilary/maternal_age_and_recombination/")


predict_y<-function(a0=NA,mu_m=NA,beta_Age,age,cohort,beta_Cohort=NA,p_by_cohort=NA,model){
    if(model==1|model==1.3|model==1.32|model==1.33|model==1.7){
        my_beta_Cohort <- c(0,beta_Cohort)
        return(a0+beta_Age[cohort]*age+my_beta_Cohort[cohort])
    }else if(model==1.2){
        my_beta_Cohort <- c(0,beta_Cohort)
        return(a0+beta_Age*age+my_beta_Cohort[cohort])
    }else if(model==1.4){
        return(a0+beta_Age[cohort]*age)
    }   else if(model ==1.5){
        return(exp(0.1*(a0+beta_Age[cohort]*age)))
    }   else if(model ==1.52|model==1.53){
        return(exp(a0+beta_Age[cohort]*age))
    }else if(model==1.6){
        return(a0+beta_Age[cohort]*age)
    }else if(model==2){
        my_beta_Cohort <- c(0,beta_Cohort)
        return(mu_m + beta_Age[cohort]*age+my_beta_Cohort[cohort])
    }else if(model==2.2){
        my_beta_Cohort <- c(0,beta_Cohort)
        return(mu_m + beta_Age*age+my_beta_Cohort[cohort])
    }else if(model ==3){
        return(exp(0.1*(a0+beta_Age[cohort]*age))*p_by_cohort[cohort])
    }else if(model==3.2){
        return(exp(0.1*(a0+beta_Age*age))*p_by_cohort[cohort])
    }
}



simulate_y<-function(N,a0=NA,mu_m=NA,beta_Age,age,cohort,beta_Cohort=NA,p_by_cohort=NA,tausq=NA,omega=NA,model){
    if(model==1|model==1.3|model==1.32|model==1.33){#cohort-specific age effects, from different or from common distribution
        my_beta_Cohort <- c(0,beta_Cohort)
        return(rnorm(N,mean=(a0+beta_Age[cohort]*age+my_beta_Cohort[cohort]),sd=sqrt(tausq)))
    }else if(model==1.2){#commmon age effect
        my_beta_Cohort <- c(0,beta_Cohort)
        return(rnorm(N,mean=(a0+beta_Age*age+my_beta_Cohort[cohort]),sd=sqrt(tausq)))
    }else if(model==1.4){#no cohort effect
        return(rnorm(N,mean=(a0+beta_Age[cohort]*age),sd=sqrt(tausq)))
    }else if(model ==1.5){ #don't think we can parallelize this
            my_beta_Cohort <- c(0,beta_Cohort)
        draws=sapply(1:N,function(i){
            return(rnbinom(1,size=exp(0.1*(a0[i]+beta_Age[cohort[i]]*age[i]+my_beta_Cohort[cohort[i]]))/(omega-1),prob=1/omega))
        })
        return(draws)
    }else if(model ==1.52|model==1.53){ #don't think we can parallelize this
            my_beta_Cohort <- c(0,beta_Cohort)
            draws=sapply(1:N,function(i){
                return(rnbinom(1,size=exp(a0[i]+beta_Age[cohort[i]]*age[i]+my_beta_Cohort[cohort[i]])/(omega-1),prob=1/omega))
            })
        return(draws)
    } else if(model==1.6){#cohort-specific age effects, from different or from common distribution
        return(rnorm(N,mean=(a0+beta_Age[cohort]*age),sd=sqrt(tausq)))
    }else if(model==1.7){#like model 1.3, but Y~student_t
        #see Gelman p 581 on how to simulate a non-central t distribution
        nu=5
        z=rnorm(N,mean=0,sd=1)
        x=rchisq(N,df=nu)
        my_beta_Cohort <- c(0,beta_Cohort)
        means=a0+beta_Age[cohort]*age+my_beta_Cohort[cohort] 
        draws = means + sqrt(tausq)*z*sqrt(nu/x)
        return(draws)
    }else if(model==2){#cohort-specific age effects from different distirbutions, no family effect
        my_beta_Cohort <- c(0,beta_Cohort)
        return(rnorm(N,mean=(mu_m+beta_Age[cohort]*age+my_beta_Cohort[cohort]),sd=sqrt(tausq)))
    }else if(model==2.2){#cohort-specific age effects from same distribution, no family effect
        my_beta_Cohort <- c(0,beta_Cohort)
        return(rnorm(N,mean=(mu_m+beta_Age*age+my_beta_Cohort[cohort]),sd=sqrt(tausq)))
    }else if(model ==3){ #don't think we can parallelize this
        draws=sapply(1:N,function(i){
            return(rnbinom(1,size=exp(0.1*(a0[i]+beta_Age[cohort[i]]*age[i]))/(omega-1),prob=(1/(p_by_cohort[cohort[i]] * (omega-1) + 1))))
        })
        return(draws)
    }else if(model==3.2){
        draws=sapply(1:N,function(i){
            return(rnbinom(1,size=exp(0.1*(a0[i]+beta_Age*age[i]))/(omega-1),prob=(1/(p_by_cohort[cohort[i]] * (omega-1) + 1))))
        })
        return(draws)
    }
}

simulate_y_for_estimation<-function(N,family=NA,mu_m=NA,beta_Age,age,cohort,beta_Cohort=NA,p_by_cohort=NA,tausq=NA,omega=NA,model,sigmasq_m=NA){
    if(model==1.2){#commmon age effect
        my_beta_Cohort <- c(0,beta_Cohort)
        a0.by.family=rnorm(length(unique(family)),mean=mu_m,sd=sqrt(sigmasq_m))
        a0=a0.by.family[family]
        return(rnorm(N,mean=(a0+beta_Age*age+my_beta_Cohort[cohort]),sd=sqrt(tausq)))
    } else if(model==1.6){#cohort-specific age effects, from different or from common distribution
        a0.by.family=rnorm(length(unique(family)),mean=mu_m[cohort_by_family],sd=sqrt(sigmasq_m))
        a0=a0.by.family[family]
        return(rnorm(N,mean=(a0+beta_Age[cohort]*age),sd=sqrt(tausq)))
    }
}

y=data1.mat2$y
J=length(y)
cohort=data1.mat2$cohort
cohort.codes=read.delim("RSTAN_output_with_NTR_v2/key_for_maternal_cohorts_to_include_in_model_1.txt",header=T,stringsAsFactors=F)
Age=data1.mat2$Age

#true.beta.age=0.06
true.beta.age=as.numeric(argv[2])

    
if(argv[1]<4){
#Model 1.2 #common age effect
    description="model1.2.maternal"
    my.model=1.2
    mydir="RSTAN_output_with_NTR_v2/model1.2"
    all.estimates=read.delim(paste0(mydir,"/parameter_estimates.",description,".txt"),header=T)
        

    if(argv[1]==1){
        my.estim="MAP"
    } else if (argv[1]==2){
        my.estim="highest_lp"
    } else if(argv[1]==3){
        my.estim="mean_of_posterior"
    }
#simulate the data
    y.sim <- simulate_y_for_estimation(N=J,family=data1.mat2$family,mu_m=all.estimates["mu_m",my.estim],beta_Age=true.beta.age,age=Age,cohort=cohort,
                                       beta_Cohort=all.estimates[grep("beta_Cohort",rownames(all.estimates)),my.estim],tausq=all.estimates["tausq",my.estim],model=my.model,sigmasq_m=all.estimates["sigmasq_m",my.estim])
    
    simulated.data <- data1.mat2
    simulated.data$y <- y.sim

    #re-estimate parameters under that same model
    model1.2.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.2.stan",data=simulated.data,iter=10,chains=0,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0"))
    model1.2.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.2.mat.compiled, data = simulated.data,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0")))
    model1.2.mat <- sflist2stanfit(model1.2.mat.list)
    print(model1.2.mat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save.image(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NTR_v2/simulations/RSTAN.model1.2.mat.fit_on_data_simulated_using_",my.estim,"_and_beta_Age_",true.beta.age,".RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NTR_v2/simulations/RSTAN.model1.2.mat.fit_on_data_simulated_using_",my.estim,".RData"))

    stanmodel=stan_model("/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.2.stan")

    ###look at the estimates
    mysim<-extract(model1.2.mat,permuted=T)
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp_)])
        } else {
            return(x[which.max(mysim$lp_),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
    umap <- optimizing(stanmodel, data=simulated.data, init=best.draw)
    new.estimates=cbind(c(umap$par["beta_Age"],umap$par[grep("^beta_Cohort",names(umap$par),perl=T)],tausq=umap$par["tausq"],umap$par["sigmasq_m"],umap$par["mu_m"]),
        c(best.draw$beta_Age,best.draw$beta_Cohort, best.draw$tausq,best.draw$sigmasq_m,best.draw$mu_m),
        c(mean(mysim$beta_Age),apply(mysim$beta_Cohort,2,mean),mean(mysim$tausq),mean(mysim$sigmasq_m),mean(mysim$mu_m)),
        c(quantile(mysim$beta_Age,probs=0.025),apply(mysim$beta_Cohort,2,quantile,probs=0.025),quantile(mysim$tausq,probs=0.025),quantile(mysim$sigmasq_m,probs=0.025),quantile(mysim$mu_m,probs=0.025)),
        c(quantile(mysim$beta_Age,probs=0.975),apply(mysim$beta_Cohort,2,quantile,probs=0.975),quantile(mysim$tausq,probs=0.975),quantile(mysim$sigmasq_m,probs=0.975),quantile(mysim$mu_m,probs=0.975)))
    rownames(new.estimates)=c("beta_Age",paste0("beta_Cohort[",1:ncol(mysim$beta_Cohort),"]"),"tausq","sigmasq_m","mu_m")
    colnames(new.estimates)=c("MAP","highest_lp","mean_of_posterior","2.5th_percentile_posterior","97.5th_percentile_posterior")
write.table(new.estimates,paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NTR_v2/simulations/parameter_estimates.",description,".fit_on_data_simulated_using_",my.estim,"_and_beta_Age_",true.beta.age,".txt"),quote=F,sep="\t")
                       
} else if(argv[1]<7){    
#### Model 1.6 -- beta_Age drawn from distribution; alphas drawn from N(mu_cohort,sigmasq)
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    
    description="model1.6.maternal"
    my.model=1.6

    mydir="RSTAN_output_with_NTR_v2/model1.6"
    all.estimates=read.delim(paste0(mydir,"/parameter_estimates.",description,".txt"),header=T)
    if(argv[1]==4){
        my.estim="MAP"
    } else if (argv[1]==5){
        my.estim="highest_lp"
    } else if(argv[1]==6){
        my.estim="mean_of_posterior"
    }
#simulate the data
    y.sim <- simulate_y_for_estimation(N=J,family=data1.mat2$family,mu_m=all.estimates[grep("mu_m",rownames(all.estimates)),my.estim],beta_Age=rep(true.beta.age,length(grep("beta_Age",rownames(all.estimates)))),age=Age,
cohort=cohort,tausq=all.estimates["tausq",my.estim],model=my.model,sigmasq_m=all.estimates["sigmasq_m",my.estim])
    
    simulated.data <- data1.mat2
    simulated.data$y <- y.sim

    #re-estimate parameters under that same model
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=simulated.data,iter=10000,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = simulated.data,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    
    print(model1.6.mat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)

  save.image(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NTR_v2/simulations/RSTAN.model1.6.mat.fit_on_data_simulated_using_",my.estim,"_and_beta_Age_",true.beta.age,".RData"))
 #   load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NTR_v2/simulations/RSTAN.model1.6.mat.fit_on_data_simulated_using_",my.estim,".RData"))  
    stanmodel=stan_model("/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan")

    ###look at the estimates
    mysim<-extract(model1.6.mat,permuted=T)
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp_)])
        } else {
            return(x[which.max(mysim$lp_),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
    umap <- optimizing(stanmodel, data=simulated.data, init=best.draw)
    new.estimates=cbind(c(umap$par[grep("beta_Age[",names(umap$par),fixed=T)],tausq=umap$par["tausq"],umap$par["sigmasq_m"],umap$par[grep("^mu_m",names(umap$par),perl=T)],umap$par["beta_global"],
        umap$par["sigmasq_global"]), c(best.draw$beta_Age, best.draw$tausq,best.draw$sigmasq_m,best.draw$mu_m,best.draw$beta_global,best.draw$sigmasq_global),
        c(apply(mysim$beta_Age,2,mean),mean(mysim$tausq),mean(mysim$sigmasq_m),apply(mysim$mu_m,2,mean),mean(mysim$beta_global),mean(mysim$sigmasq_global)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.025),quantile(mysim$tausq,probs=0.025),quantile(mysim$sigmasq_m,probs=0.025),apply(mysim$mu_m,2,quantile,probs=0.025),quantile(mysim$beta_global,probs=0.025),
          quantile(mysim$sigmasq_global,probs=0.025)), c(apply(mysim$beta_Age,2,quantile,probs=0.975),quantile(mysim$tausq,probs=0.975),quantile(mysim$sigmasq_m,probs=0.975),apply(mysim$mu_m,2,quantile,probs=0.975),
                                            quantile(mysim$beta_global,probs=0.975),quantile(mysim$sigmasq_global,probs=0.975)))
    rownames(new.estimates)=c(paste0("beta_Age[",1:ncol(mysim$beta_Age),"]"),"tausq","sigmasq_m",paste0("mu_m[",1:ncol(mysim$mu_m),"]"),"beta_global","sigmasq_global")
    colnames(new.estimates)=c("MAP","highest_lp","mean_of_posterior","2.5th_percentile_posterior","97.5th_percentile_posterior")
    
write.table(new.estimates,paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NTR_v2/simulations/parameter_estimates.",description,".fit_on_data_simulated_using_",my.estim,"_and_beta_Age_",true.beta.age,".txt"),quote=F,sep="\t")

}


