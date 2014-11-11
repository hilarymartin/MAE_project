library(rstan)
library(scales)
load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.more_stringent.RData")
argv <- commandArgs(trailingOnly = TRUE)
print(argv)
setwd("/well/donnelly/hilary/maternal_age_and_recombination/")

if(argv[1]==1){
    load("RSTAN_output_on_duoHMM_more_stringent/RSTAN.model3.adjusted_priors.all_chains.mat.including_alpha.RData")
    mysim<-extract(model3.mat,permuted=T)
    y=data2.mat2$y
    cohort=data2.mat2$cohort
    mother=data2.mat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data2.mat2$Age
    description="model3.maternal"
    my.model=3
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model3"
    
    #stanmodel=stan_model("/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan")
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
    best.draw[["inv_omega"]]=1/best.draw[["omega"]]
    #umap <- optimizing(stanmodel, data=data2.mat2, init=best.draw)
    
}else if(argv[1]==2){
    load("RSTAN_output_on_duoHMM_more_stringent/RSTAN.model3.adjusted_priors.all_chains.pat.including_alpha.RData")
    mysim<-extract(model3.pat,permuted=T)
    y=data2.pat2$y
    cohort=data2.pat2$cohort
    mother=data2.pat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data2.pat2$Age
    description="model3.paternal"
    my.model=3
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model3"


    #stanmodel=stan_model("/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan")
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
    best.draw[["inv_omega"]]=1/best.draw[["omega"]]
    #umap <- optimizing(stanmodel, data=data2.pat2, init=best.draw)
    
    
}else if(argv[1]==3){
    load("RSTAN_output_on_duoHMM_more_stringent/RSTAN.model3.2.adjusted_priors.all_chains.mat.including_alpha.RData")
    mysim<-extract(model3.2.mat,permuted=T)
    y=data2.mat2$y
    cohort=data2.mat2$cohort
    mother=data2.mat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data2.mat2$Age
    description="model3.2.maternal"
    my.model=3.2
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model3.2"

    #stanmodel=stan_model("/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan")
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
    best.draw[["inv_omega"]]=1/best.draw[["omega"]]
    #umap <- optimizing(stanmodel, data=data2.mat2, init=best.draw)

}else if(argv[1]==4){
    load("RSTAN_output_on_duoHMM_more_stringent/RSTAN.model3.2.adjusted_priors.all_chains.pat.including_alpha.RData")
    mysim<-extract(model3.2.pat,permuted=T)
    y=data2.pat2$y
    cohort=data2.pat2$cohort
    mother=data2.pat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data2.pat2$Age
    description="model3.2.paternal"
    my.model=3.2
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model3.2"

    #stanmodel=stan_model("/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan")
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
    best.draw[["inv_omega"]]=1/best.draw[["omega"]]
    #umap <- optimizing(stanmodel, data=data2.pat2, init=best.draw)

} else if(argv[1]==5){        
#### Model 1.5 -- negative binomial model for informative families only
   load("RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.5.all_chains.mat.including_alpha.RData")
    mysim<-extract(model1.5.mat,permuted=T)
    y=data1.mat2$y
    cohort=data1.mat2$cohort
    mother=data1.mat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data1.mat2$Age
    description="model1.5.maternal"
    my.model=1.5
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model1.5"

    
    #stanmodel=stan_model("/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.stan")

    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
   best.draw[["inv_omega"]]=1/best.draw[["omega"]]
   #umap <- optimizing(stanmodel, data=data1.mat2, init=best.draw)

    
} else if(argv[1]==6){
        
#### Model 1.5 -- negative binomial model for informative families only
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.5.all_chains.pat.including_alpha.RData")
    mysim<-extract(model1.5.pat,permuted=T)
    y=data1.pat2$y
    cohort=data1.pat2$cohort
    mother=data1.pat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data1.pat2$Age
    description="model1.5.paternal"
    my.model=1.5
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model1.5"

    
    #stanmodel=stan_model("/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.stan")
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
      best.draw[["inv_omega"]]=1/best.draw[["omega"]]
    #umap <- optimizing(stanmodel, data=data1.pat2, init=best.draw)
   
} else if(argv[1]==7){    
#### Model 1.6 -- beta_Age drawn from distribution; alphas drawn from N(mu_cohort,sigmasq)
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.6.all_chains.mat.including_alpha.mu_m_N_41_100.sigmasq_m_IG_2_40.RData")
    mysim<-extract(model1.6.mat,permuted=T)
    y=data1.mat2$y
    cohort=data1.mat2$cohort
    mother=data1.mat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data1.mat2$Age
    description="model1.6.maternal"
    my.model=1.6
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model1.6"

     #stanmodel=stan_model("/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan")
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
    #umap <- optimizing(stanmodel, data=data1.mat2, init=best.draw)
        
} else if(argv[1]==8){
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.6.all_chains.pat.including_alpha.mu_m_N_27_64.sigmasq_m_IG_2_15.RData")
    mysim<-extract(model1.6.pat,permuted=T)
    y=data1.pat2$y
    cohort=data1.pat2$cohort
    mother=data1.pat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data1.pat2$Age
    description="model1.6.paternal"
    my.model=1.6
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model1.6"

    #stanmodel=stan_model("/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan")
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
    #umap <- optimizing(stanmodel, data=data1.pat2, init=best.draw)

} else if(argv[1]==9){    
#### Model 1.6.2 -- common beta_Age; alphas drawn from N(mu_cohort,sigmasq)
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family

    mean_alpha_prior=38
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta

    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_41_100.sigmasq_m_IG_2_40.RData"))
    
    mysim<-extract(model1.6.mat,permuted=T)
    y=data1.mat2$y
    cohort=data1.mat2$cohort
    mother=data1.mat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data1.mat2$Age
    description="model1.6.2.maternal.mu_m_N_41_100.sigmasq_m_IG_2_40"
    my.model=1.62
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model1.62"

    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
    save(best.draw,file=paste0(mydir,"/",description,".best_draw.RData"))
    
} else if(argv[1]==10){
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    mean_alpha_prior=30
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta

    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_27_64.sigmasq_m_IG_2_15.RData"))
    mysim<-extract(model1.6.pat,permuted=T)
    y=data1.pat2$y
    cohort=data1.pat2$cohort
    mother=data1.pat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data1.pat2$Age
    description="model1.6.2.paternal.mu_m_N_27_64.sigmasq_m_IG_2_15"
    my.model=1.62
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model1.62"

    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]

} else if(argv[1]==11){        
#### Model 1.5.4 -- negative binomial model for informative families only, no cohort effect but family effects drawn from distribution with cohort-specific mean; common beta_Age
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
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.5.4.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))

    mysim<-extract(model1.5.mat,permuted=T)
    y=data1.mat2$y
    cohort=data1.mat2$cohort
    mother=data1.mat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data1.mat2$Age
    description="model1.5.4.maternal"
    my.model=1.54
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model1.54"
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
   best.draw[["inv_omega"]]=1/best.draw[["omega"]]
} else if(argv[1]==12){
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
#### Model 1.5 -- negative binomial model for informative families only
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/RSTAN.modell1.5.4.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.5.pat,permuted=T)
    y=data1.pat2$y
    cohort=data1.pat2$cohort
    mother=data1.pat2$family
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.txt",header=T,stringsAsFactors=F)
    Age=data1.pat2$Age
    description="model1.5.4.paternal"
    my.model=1.54
    J=length(y)
    mydir="RSTAN_output_on_duoHMM_more_stringent/model1.54"

   
    best.draw=sapply(1:(length(mysim)-1),function(i){
        x=mysim[[i]]
        if(class(x)=="array"){
            return(x[which.max(mysim$lp__)])
        } else {
            return(x[which.max(mysim$lp__),])
        }
    })
    names(best.draw)=names(mysim)[1:(length(mysim)-1)]
     best.draw[["inv_omega"]]=1/best.draw[["omega"]]
}
    

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
    }   else if(model ==1.54){
        return(exp(0.1*(a0+beta_Age*age)))
    }   else if(model ==1.52|model==1.53){
        return(exp(a0+beta_Age[cohort]*age))
    }else if(model==1.6){
        return(a0+beta_Age[cohort]*age)
    }else if(model==1.62){
        return(a0+beta_Age*age)
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
    }else if(model ==1.54){ #don't think we can parallelize this
        draws=sapply(1:N,function(i){
            return(rnbinom(1,size=exp(0.1*(a0[i]+beta_Age*age[i]))/(omega-1),prob=1/omega))
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
    }else if(model==1.62){#common age effect, family effects from distribution with cohort-specific mean
        return(rnorm(N,mean=(a0+beta_Age*age),sd=sqrt(tausq)))
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

#test statistic sensitive to asymmetry in y
test_10pc_90pc <- function (y){
return (abs(quantile(y,0.9)-median(y))-abs(quantile(y,0.1)-median(y)))
}


cat("write out different estimates\n")
if(my.model==1){
    all.estimates=cbind(        c(best.draw$beta_Age,best.draw$beta_Cohort, best.draw$tausq,best.draw$sigmasq_m,best.draw$mu_m),
        c(apply(mysim$beta_Age,2,mean),apply(mysim$beta_Cohort,2,mean),mean(mysim$tausq),mean(mysim$sigmasq_m),mean(mysim$mu_m)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.025),apply(mysim$beta_Cohort,2,quantile,probs=0.025),quantile(mysim$tausq,probs=0.025),quantile(mysim$sigmasq_m,probs=0.025),quantile(mysim$mu_m,probs=0.025)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.975),apply(mysim$beta_Cohort,2,quantile,probs=0.975),quantile(mysim$tausq,probs=0.975),quantile(mysim$sigmasq_m,probs=0.975),quantile(mysim$mu_m,probs=0.975)))
    rownames(all.estimates)=c(paste0("beta_Age[",1:ncol(mysim$beta_Age),"]"),paste0("beta_Cohort[",1:ncol(mysim$beta_Cohort),"]"),"tausq","sigmasq_m","mu_m")
}else if(my.model==1.2){
    all.estimates=cbind(        c(best.draw$beta_Age,best.draw$beta_Cohort, best.draw$tausq,best.draw$sigmasq_m,best.draw$mu_m),
        c(mean(mysim$beta_Age),apply(mysim$beta_Cohort,2,mean),mean(mysim$tausq),mean(mysim$sigmasq_m),mean(mysim$mu_m)),
        c(quantile(mysim$beta_Age,probs=0.025),apply(mysim$beta_Cohort,2,quantile,probs=0.025),quantile(mysim$tausq,probs=0.025),quantile(mysim$sigmasq_m,probs=0.025),quantile(mysim$mu_m,probs=0.025)),
        c(quantile(mysim$beta_Age,probs=0.975),apply(mysim$beta_Cohort,2,quantile,probs=0.975),quantile(mysim$tausq,probs=0.975),quantile(mysim$sigmasq_m,probs=0.975),quantile(mysim$mu_m,probs=0.975)))
    rownames(all.estimates)=c("beta_Age",paste0("beta_Cohort[",1:ncol(mysim$beta_Cohort),"]"),"tausq","sigmasq_m","mu_m")
}else if(my.model==1.3|my.model==1.32|my.model==1.33|my.model==1.7){
    all.estimates=cbind(        c(best.draw$beta_Age,best.draw$beta_Cohort, best.draw$tausq,best.draw$sigmasq_m,best.draw$mu_m,best.draw$beta_global,best.draw$sigmasq_global),
        c(apply(mysim$beta_Age,2,mean),apply(mysim$beta_Cohort,2,mean),mean(mysim$tausq),mean(mysim$sigmasq_m),mean(mysim$mu_m),mean(mysim$beta_global),mean(mysim$sigmasq_global)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.025),apply(mysim$beta_Cohort,2,quantile,probs=0.025),quantile(mysim$tausq,probs=0.025),quantile(mysim$sigmasq_m,probs=0.025),
          quantile(mysim$mu_m,probs=0.025),quantile(mysim$beta_global,probs=0.025),quantile(mysim$sigmasq_global,probs=0.025)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.975),apply(mysim$beta_Cohort,2,quantile,probs=0.975),quantile(mysim$tausq,probs=0.975),quantile(mysim$sigmasq_m,probs=0.975),
          quantile(mysim$mu_m,probs=0.975),quantile(mysim$beta_global,probs=0.975),quantile(mysim$sigmasq_global,probs=0.975)))
    rownames(all.estimates)=c(paste0("beta_Age[",1:ncol(mysim$beta_Age),"]"),paste0("beta_Cohort[",1:ncol(mysim$beta_Cohort),"]"),"tausq","sigmasq_m","mu_m","beta_global","sigmasq_global")
}else if(my.model==1.4){
    all.estimates=cbind(        c(best.draw$beta_Age,best.draw$tausq,best.draw$sigmasq_m,best.draw$mu_m,best.draw$beta_global,best.draw$sigmasq_global),
        c(apply(mysim$beta_Age,2,mean),mean(mysim$tausq),mean(mysim$sigmasq_m),mean(mysim$mu_m),mean(mysim$beta_global),mean(mysim$sigmasq_global)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.025),quantile(mysim$tausq,probs=0.025),quantile(mysim$sigmasq_m,probs=0.025),quantile(mysim$mu_m,probs=0.025),quantile(mysim$beta_global,probs=0.025),quantile(mysim$sigmasq_global,probs=0.025)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.975),quantile(mysim$tausq,probs=0.975),quantile(mysim$sigmasq_m,probs=0.975),quantile(mysim$mu_m,probs=0.975),quantile(mysim$beta_global,probs=0.975),quantile(mysim$sigmasq_global,probs=0.975)))
    rownames(all.estimates)=c(paste0("beta_Age[",1:ncol(mysim$beta_Age),"]"),"tausq","sigmasq_m","mu_m","beta_global","sigmasq_global")
}else if (my.model==1.5|my.model==1.52|my.model==1.53){
    all.estimates=cbind(        c(best.draw$beta_Age,best.draw$beta_global,best.draw$sigmasq_global,best.draw$beta_Cohort,best.draw$mu_m,best.draw$sigmasq_m,best.draw$omega),
        c(apply(mysim$beta_Age,2,mean),mean(mysim$beta_global),mean(mysim$sigmasq_global),apply(mysim$beta_Cohort,2,mean),mean(mysim$mu_m),mean(mysim$sigmasq_m),mean(mysim$omega)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.025),quantile(mysim$beta_global,probs=0.025),quantile(mysim$sigmasq_global,probs=0.025),apply(mysim$beta_Cohort,2,quantile,probs=0.025),quantile(mysim$mu_m,probs=0.025),quantile(mysim$sigmasq_m,probs=0.025),
          quantile(mysim$omega,probs=0.025)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.975),quantile(mysim$beta_global,probs=0.975),quantile(mysim$sigmasq_global,probs=0.975),apply(mysim$beta_Cohort,2,quantile,probs=0.975),quantile(mysim$mu_m,probs=0.975),quantile(mysim$sigmasq_m,probs=0.975),
          quantile(mysim$omega,probs=0.975)))
    rownames(all.estimates)=c(paste0("beta_Age[",1:ncol(mysim$beta_Age),"]"),"beta_global","sigmasq_global",paste0("beta_Cohort[",1:ncol(mysim$beta_Cohort),"]"),"mu_m","sigmasq_m","omega")
}else if (my.model==1.54){
    all.estimates=cbind(c(best.draw$beta_Age,best.draw$mu_m,best.draw$sigmasq_m,best.draw$omega),
        c(mean(mysim$beta_Age),apply(mysim$mu_m,2,mean),apply(mysim$sigmasq_m,2,mean),mean(mysim$omega)),
        c(quantile(mysim$beta_Age,probs=0.025),apply(mysim$mu_m,2,quantile,probs=0.025),apply(mysim$sigmasq_m,2,quantile,probs=0.025),quantile(mysim$omega,probs=0.025)),
        c(quantile(mysim$beta_Age,probs=0.975),apply(mysim$mu_m,2,quantile,probs=0.975),apply(mysim$sigmasq_m,2,quantile,probs=0.975),quantile(mysim$omega,probs=0.975)))
    rownames(all.estimates)=c("beta_Age",paste0("mu_m[",1:ncol(mysim$mu_m),"]"),paste0("sigmasq_m[",1:ncol(mysim$sigmasq_m),"]"),"omega")
}else if (my.model==1.6){
    all.estimates=cbind(        c(best.draw$beta_Age, best.draw$tausq,best.draw$sigmasq_m,best.draw$mu_m,best.draw$beta_global,best.draw$sigmasq_global),
        c(apply(mysim$beta_Age,2,mean),mean(mysim$tausq),apply(mysim$sigmasq_m,2,mean),apply(mysim$mu_m,2,mean),mean(mysim$beta_global),mean(mysim$sigmasq_global)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.025),quantile(mysim$tausq,probs=0.025),apply(mysim$sigmasq_m,2,quantile,probs=0.025),apply(mysim$mu_m,2,quantile,probs=0.025),quantile(mysim$beta_global,probs=0.025),quantile(mysim$sigmasq_global,probs=0.025)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.975),quantile(mysim$tausq,probs=0.975),apply(mysim$sigmasq_m,2,quantile,probs=0.975),apply(mysim$mu_m,2,quantile,probs=0.975),quantile(mysim$beta_global,probs=0.975),quantile(mysim$sigmasq_global,probs=0.975)))
    rownames(all.estimates)=c(paste0("beta_Age[",1:ncol(mysim$beta_Age),"]"),"tausq",paste0("sigmasq_m[",1:ncol(mysim$sigmasq_m),"]"),paste0("mu_m[",1:ncol(mysim$mu_m),"]"),"beta_global","sigmasq_global")
}else if (my.model==1.62){
    all.estimates=cbind(c(best.draw$beta_Age, best.draw$tausq,best.draw$sigmasq_m,best.draw$mu_m),
        c(mean(mysim$beta_Age),mean(mysim$tausq),apply(mysim$sigmasq_m,2,mean),apply(mysim$mu_m,2,mean)),
        c(quantile(mysim$beta_Age,probs=0.025),quantile(mysim$tausq,probs=0.025),apply(mysim$sigmasq_m,2,quantile,probs=0.025),apply(mysim$mu_m,2,quantile,probs=0.025)),
        c(quantile(mysim$beta_Age,probs=0.975),quantile(mysim$tausq,probs=0.975),apply(mysim$sigmasq_m,2,quantile,probs=0.975),apply(mysim$mu_m,2,quantile,probs=0.975)))
    rownames(all.estimates)=c("beta_Age","tausq",paste0("sigmasq_[m",1:ncol(mysim$sigmasq_m),"]"),paste0("mu_m[",1:ncol(mysim$mu_m),"]"))
} else if(my.model==2){
    all.estimates=cbind(        c(best.draw$beta_Age,best.draw$beta_Cohort, best.draw$tausq,best.draw$mu_m),
        c(apply(mysim$beta_Age,2,mean),apply(mysim$beta_Cohort,2,mean),mean(mysim$tausq),mean(mysim$mu_m)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.025),apply(mysim$beta_Cohort,2,quantile,probs=0.025),quantile(mysim$tausq,probs=0.025),quantile(mysim$mu_m,probs=0.025)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.975),apply(mysim$beta_Cohort,2,quantile,probs=0.975),quantile(mysim$tausq,probs=0.975),quantile(mysim$mu_m,probs=0.975)))
    rownames(all.estimates)=c(paste0("beta_Age[",1:ncol(mysim$beta_Age),"]"),paste0("beta_Cohort[",1:ncol(mysim$beta_Cohort),"]"),"tausq","mu_m")
}else if(my.model==2.2){
    all.estimates=cbind(        c(best.draw$beta_Age,best.draw$beta_Cohort, best.draw$tausq,best.draw$mu_m),
        c(mean(mysim$beta_Age),apply(mysim$beta_Cohort,2,mean),mean(mysim$tausq),mean(mysim$mu_m)),
        c(quantile(mysim$beta_Age,probs=0.025),apply(mysim$beta_Cohort,2,quantile,probs=0.025),quantile(mysim$tausq,probs=0.025),quantile(mysim$mu_m,probs=0.025)),
        c(quantile(mysim$beta_Age,probs=0.975),apply(mysim$beta_Cohort,2,quantile,probs=0.975),quantile(mysim$tausq,probs=0.975),quantile(mysim$mu_m,probs=0.975)))
    rownames(all.estimates)=c("beta_Age",paste0("beta_Cohort[",1:ncol(mysim$beta_Cohort),"]"),"tausq","mu_m")
}else if(my.model==3){
    all.estimates=cbind(        c(best.draw$beta_Age,best.draw$beta_global,best.draw$sigmasq_global,best.draw$p_by_cohort,best.draw$mu_m,best.draw$sigmasq_m,best.draw$omega),
        c(apply(mysim$beta_Age,2,mean),mean(mysim$beta_global),mean(mysim$sigmasq_global),apply(mysim$p_by_cohort,2,mean),mean(mysim$mu_m),mean(mysim$sigmasq_m),mean(mysim$omega)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.025),quantile(mysim$beta_global,probs=0.025),quantile(mysim$sigmasq_global,probs=0.025),apply(mysim$p_by_cohort,2,quantile,probs=0.025),quantile(mysim$mu_m,probs=0.025),quantile(mysim$sigmasq_m,probs=0.025),
          quantile(mysim$omega,probs=0.025)),
        c(apply(mysim$beta_Age,2,quantile,probs=0.975),quantile(mysim$beta_global,probs=0.975),quantile(mysim$sigmasq_global,probs=0.975),apply(mysim$p_by_cohort,2,quantile,probs=0.975),quantile(mysim$mu_m,probs=0.975),quantile(mysim$sigmasq_m,probs=0.975),
          quantile(mysim$omega,probs=0.975)))
    rownames(all.estimates)=c(paste0("beta_Age[",1:ncol(mysim$beta_Age),"]"),"beta_global","sigmasq_global",paste0("p_by_cohort[",1:ncol(mysim$p_by_cohort),"]"),"mu_m","sigmasq_m","omega")
}else if(my.model==3.2){
    all.estimates=cbind(        c(best.draw$beta_Age,best.draw$p_by_cohort,best.draw$mu_m,best.draw$sigmasq_m,best.draw$omega),
        c(mean(mysim$beta_Age),apply(mysim$p_by_cohort,2,mean),mean(mysim$mu_m),mean(mysim$sigmasq_m),mean(mysim$omega)),
        c(quantile(mysim$beta_Age,probs=0.025),apply(mysim$p_by_cohort,2,quantile,probs=0.025),quantile(mysim$mu_m,probs=0.025),quantile(mysim$sigmasq_m,probs=0.025),
          quantile(mysim$omega,probs=0.025)),
        c(quantile(mysim$beta_Age,probs=0.975),apply(mysim$p_by_cohort,2,quantile,probs=0.975),quantile(mysim$mu_m,probs=0.975),quantile(mysim$sigmasq_m,probs=0.975),
          quantile(mysim$omega,probs=0.975)))
    rownames(all.estimates)=c("beta_Age",paste0("p_by_cohort[",1:ncol(mysim$p_by_cohort),"]"),"mu_m","sigmasq_m","omega")
}
colnames(all.estimates)=c("highest_lp","mean_of_posterior","2.5th_percentile_posterior","97.5th_percentile_posterior")
write.table(all.estimates,paste0(mydir,"/parameter_estimates.",description,".txt"),quote=F,sep="\t")

cat("plot histograms of observed data and data simulated using parameters from draw with highest -log likelihood\n")
best.draw.sims=NULL
for(i in 1:100){
    if(my.model==1|my.model==1.3|my.model==1.32|my.model==1.33|my.model==1.7){
        y.sim <- simulate_y(N=J,a0=best.draw$a0[mother],beta_Age=best.draw$beta_Age,age=Age,tausq=best.draw$tausq,cohort=cohort,beta_Cohort=best.draw$beta_Cohort,model=my.model)
    }else if(my.model==1.2){
        y.sim <- simulate_y(N=J,a0=best.draw$a0[mother],beta_Age=best.draw$beta_Age,age=Age,tausq=best.draw$tausq,cohort=cohort,beta_Cohort=best.draw$beta_Cohort,model=my.model)
    }else if(my.model==1.4|my.model==1.6|my.model==1.62){
        y.sim <- simulate_y(N=J,a0=best.draw$a0[mother],beta_Age=best.draw$beta_Age,age=Age,tausq=best.draw$tausq,cohort=cohort,model=my.model)
    }else if(my.model==1.5|my.model==1.52|my.model==1.53){
        y.sim <- simulate_y(N=J,a0=best.draw$exp_a0[mother],beta_Age=best.draw$beta_Age,age=Age,omega=best.draw$omega,cohort=cohort,beta_Cohort=best.draw$beta_Cohort,model=my.model)
    }else if(my.model==1.54){
        y.sim <- simulate_y(N=J,a0=best.draw$exp_a0[mother],beta_Age=best.draw$beta_Age,age=Age,omega=best.draw$omega,cohort=cohort,model=my.model)
    }else if(my.model==2){
        y.sim <- simulate_y(N=J,mu_m=best.draw$mu_m,beta_Age=best.draw$beta_Age,age=Age,tausq=best.draw$tausq,cohort=cohort,beta_Cohort=best.draw$beta_Cohort,model=my.model)
    }else if(my.model==2.2){
        y.sim <- simulate_y(N=J,mu_m=best.draw$mu_m,beta_Age=best.draw$beta_Age,age=Age,tausq=best.draw$tausq,cohort=cohort,beta_Cohort=best.draw$beta_Cohort,model=my.model)
    }else if(my.model==3){
        y.sim <- simulate_y(N=J,a0=best.draw$exp_a0[mother],beta_Age=best.draw$beta_Age,age=Age,omega=best.draw$omega,cohort=cohort,p_by_cohort=best.draw$p_by_cohort,model=my.model)
    }else if(my.model==3.2){
        y.sim <- simulate_y(N=J,a0=best.draw$exp_a0[mother],beta_Age=best.draw$beta_Age,age=Age,omega=best.draw$omega,cohort=cohort,p_by_cohort=best.draw$p_by_cohort,model=my.model)
    }
    best.draw.sims=rbind(best.draw.sims,y.sim)
}

pdf(paste(mydir,"/","ecdfs_of_observed_and_simulated_y_from_highest_lp_estimates.",description,".pdf",sep=""),height=5,width=5)
for(i in 1:nrow(best.draw.sims)){
    if(i==1){
        plot(ecdf(best.draw.sims[i,]),xlab="crossover counts",col="grey",verticals=T,main="ecdf of Y",pch=NA)
    } else {
        lines(ecdf(best.draw.sims[i,]),col="grey",verticals=T,pch=NA)
    }
}
    lines(ecdf(y),col="red",verticals=T,pch=NA)
legend("topleft",c("observed","simulated"),col=c("red","grey"),lty=1)
dev.off()


if(my.model<3){
    pdf(paste(mydir,"/","ecdfs_of_observed_and_simulated_y_from_highest_lp_estimates.by_cohort.",description,".pdf",sep=""),height=9,width=9)
    par(mfrow=c(3,3))
} else {
    pdf(paste(mydir,"/","ecdfs_of_observed_and_simulated_y_from_highest_lp_estimates.by_cohort.",description,".pdf",sep=""),height=12,width=21)
    par(mfrow=c(4,7))
}
for(c in 1:length(unique(cohort))){
for(i in 1:nrow(best.draw.sims)){
    if(i==1){
        plot(ecdf(best.draw.sims[i,cohort==c]),xlab="crossover counts",col="grey",verticals=T,main=paste0(cohort.codes[cohort.codes$Code==c,"Cohort"],", n=",sum(cohort==c)),pch=NA,cex.main=1,cex.lab=2,cex.axis=2)
    } else {
        lines(ecdf(best.draw.sims[i,cohort==c]),col="grey",verticals=T,pch=NA)
    }
}
lines(ecdf(y[cohort==c]),col="red",verticals=T,pch=NA)
if(c==1){
legend("topleft",c("observed","simulated"),col=c("red","grey"),lty=1,cex=0.8)
}
}
dev.off()

pdf(paste(mydir,"/","histograms_of_observed_and_simulated_y_from_highest_lp_estimates.",description,".pdf",sep=""),height=20,width=20)
obs.hist=hist(y,plot=F)
par(mfrow=c(4,5))
for(i in 1:20){
    sim.hist=hist(best.draw.sims[i,],plot=F)
    both.breaks=obs.hist$breaks
    if(max(best.draw.sims[i,])>max(both.breaks)){
        gaps=both.breaks[2]-both.breaks[1]
        while(both.breaks[length(both.breaks)] <max(best.draw.sims[i,])){
            both.breaks=c(both.breaks,both.breaks[length(both.breaks)]+gaps)
        }
    }
    if(min(best.draw.sims[i,])<min(both.breaks)){
        gaps=both.breaks[2]-both.breaks[1]
        while(both.breaks[1] >min(best.draw.sims[i,])){
            both.breaks=c(both.breaks[1]-gaps,both.breaks)
        }
    }
    sim.hist=hist(best.draw.sims[i,],plot=F,breaks=both.breaks)
    hist(y,main=paste0("simulation ",i),xlab="crossover counts",xlim=range(c(obs.hist$breaks,sim.hist$breaks)),cex.lab=2,cex.main=2,cex.axis=2,col=alpha("red",0.5),ylim=range(c(obs.hist$counts,sim.hist$counts)))
    hist(best.draw.sims[i,],col=alpha("blue",0.5),add=T,breaks=both.breaks)
    if(i==1){
    legend("topleft",c("observed","simulated"),fill=c(alpha("red",0.5),alpha("blue",0.5)),cex=1.5)
}
}
dev.off()




### plot histograms of summary statistics
cat("plot histograms of summary statistics\n")
mean_y <- mean(y)
median_y <- median(y)
max_y <- max(y)
min_y <- min(y)
test_10pc_90pc_y <-test_10pc_90pc(y)

mean_best.draw.sim <- apply(best.draw.sims,1,mean)
median_best.draw.sim <- apply(best.draw.sims,1,median)
max_best.draw.sim <- apply(best.draw.sims,1,max)
min_best.draw.sim <- apply(best.draw.sims,1,min)
test_10pc_90pc_best.draw.sim <-apply(best.draw.sims,1,test_10pc_90pc)

pdf(paste(mydir,"/","histogram_of_summary_statistics.highest_lp_estimates.",description,".pdf",sep=""),height=10,width=15)
par(mfrow=c(2,3))
hist(mean_best.draw.sim,xlim=range(c(mean_y,mean_best.draw.sim)),xlab="mean predicted y",main="Mean",cex.main=2,cex.axis=2,cex.lab=2)
abline(v=mean_y,col="red")
legend("topright",paste("p=",round(mean(mean_best.draw.sim>=mean(y)),3)),cex=2)

hist(median_best.draw.sim,xlim=range(c(median_y,median_best.draw.sim)),xlab="median predicted y",main="Median",cex.main=2,cex.axis=2,cex.lab=2)
abline(v=median_y,col="red")
legend("topright",paste("p=",round(mean(median_best.draw.sim>=median(y)),3)),cex=2)

hist(max_best.draw.sim,xlim=range(c(max_y,max_best.draw.sim)),xlab="maximum predicted y",main="Maximum",cex.main=2,cex.axis=2,cex.lab=2)
abline(v=max_y,col="red")
legend("topright",paste("p=",round(mean(max_best.draw.sim>=max(y)),3)),cex=2)

hist(min_best.draw.sim,xlim=range(c(min_y,min_best.draw.sim)),xlab="minimum predicted y",main="Minimum",cex.main=2,cex.axis=2,cex.lab=2)
abline(v=min_y,col="red")
legend("topleft",paste("p=",round(mean(min_best.draw.sim>=min(y)),3)),cex=2)

hist(test_10pc_90pc_best.draw.sim,xlim=range(c(test_10pc_90pc_y,test_10pc_90pc_best.draw.sim)),xlab="|90pc_y - 50pc_y|-|10pc_y - 50pc_y|",main="Asymmetry",cex.main=2,cex.axis=2,cex.lab=2)
abline(v=test_10pc_90pc_y,col="red")
legend("topright",paste("p=",round(mean(test_10pc_90pc_best.draw.sim>=test_10pc_90pc_y),3)),cex=2)

dev.off()

if(FALSE){
#plot Bayesian residuals for a few draws
pdf(paste(mydir,"/","Bayesian_residual_plot.",description,".pdf",sep=""),height=20,width=20)
par(mfrow=c(4,5))
for(i in 1:20){
s=random_draws[i]
y_predict_s= y.predictions[as.character(s),]
residuals=y-y_predict_s
plot(y_predict_s,residuals,ylab="residual",xlab="expected count",main=paste("rep",s,sep=""),cex.lab=2,cex.main=2,cex.axis=2)
abline(h=0,col="red")
}
dev.off()

#plot Bayesian residuals for a few draws,binned
pdf(paste(mydir,"/","Bayesian_binned_residual_plot.",description,".pdf",sep=""),height=20,width=20)
par(mfrow=c(4,5))
for(i in 1:20){
    s=random_draws[i]
    y_predict_s= y.predictions[as.character(s),]
    residuals=y-y_predict_s
    names(y_predict_s)=1:length(y_predict_s)
    y_predict_sorted=sort(y_predict_s)
    names(residuals)=1:length(y_predict_s)
    residuals_sorted=residuals[names(y_predict_sorted)]
    y_predict_binned=rep(NA,20)
    residuals_binned=rep(NA,20)
    increment=round(length(y_predict_s)/20)
    c=1
    for(i in 1:20){
        y_predict_binned[i]=mean(y_predict_sorted[c:min((c+increment),length(y_predict_s))])
        residuals_binned[i]=mean(residuals_sorted[c:min((c+increment),length(y_predict_s))])
        c=c+increment
    }
    plot(y_predict_binned,residuals_binned,ylab="estimated residual, binned",xlab="expected count, binned",main=paste("rep",s,sep=""),cex.lab=2,cex.main=2,cex.axis=2)
    abline(h=0,col="red")
}
dev.off()


}

#### predict using highest lp estimate
cat("predict using highest lp estimate\n")
if(my.model==1|my.model==1.3|my.model==1.32|my.model==1.33|my.model==1.7){
    y.predict.best <- predict_y(a0=best.draw$a0[mother],beta_Age=best.draw$beta_Age,age=Age,cohort=cohort,beta_Cohort=best.draw$beta_Cohort,model=my.model)
}else if(my.model==1.2){
    y.predict.best <- predict_y(a0=best.draw$a0[mother],beta_Age=best.draw$beta_Age,age=Age,cohort=cohort,beta_Cohort=best.draw$beta_Cohort,model=my.model)
}else if(my.model==1.4|my.model==1.6|my.model==1.62){
    y.predict.best <- predict_y(a0=best.draw$a0[mother],beta_Age=best.draw$beta_Age,age=Age,cohort=cohort,model=my.model)
}else if(my.model==1.5|my.model==1.52|my.model==1.53){
    y.predict.best <- predict_y(a0=best.draw$exp_a0[mother],beta_Age=best.draw$beta_Age,age=Age,cohort=cohort,beta_Cohort=best.draw$beta_Cohort,model=my.model)
}else if(my.model==1.54){
    y.predict.best <- predict_y(a0=best.draw$exp_a0[mother],beta_Age=best.draw$beta_Age,age=Age,cohort=cohort,model=my.model)
}else if(my.model==2){
    y.predict.best <- predict_y(mu_m=best.draw$mu_m,beta_Age=best.draw$beta_Age,age=Age,cohort=cohort,beta_Cohort=best.draw$beta_Cohort,model=my.model)
}else if(my.model==2.2){
    y.predict.best <- predict_y(mu_m=best.draw$mu_m,beta_Age=best.draw$beta_Age,age=Age,cohort=cohort,beta_Cohort=best.draw$beta_Cohort,model=my.model)
}else if(my.model==3){
    y.predict.best <- predict_y(a0=best.draw$exp_a0[mother],beta_Age=best.draw$beta_Age,age=Age,cohort=cohort,p_by_cohort=best.draw$p_by_cohort,model=my.model)
}else if(my.model==3.2){
    y.predict.best <- predict_y(a0=best.draw$exp_a0[mother],beta_Age=best.draw$beta_Age,age=Age,cohort=cohort,p_by_cohort=best.draw$p_by_cohort,model=my.model)
}
names(y.predict.best)=1:length(y.predict.best)
classic_residuals.best=y-y.predict.best
names(classic_residuals.best)=1:length(y.predict.best)


prediction.list=list(y.predict.best)
residual.list=list(classic_residuals.best)
mydescriptions=c(paste0("highest_lp_estimates.",description))

for(x in 1:length(prediction.list)){
y_predict=prediction.list[[x]]
classic_residuals=residual.list[[x]]
description2=mydescriptions[x]
####plot residuals against age
cat("plot residuals against age\n")
pdf(paste(mydir,"/","residual_plot_vs_age.",description2,".pdf",sep=""),height=10,width=10)
plot(Age,classic_residuals,ylab="estimated residual",xlab="Age",main="Estimated residuals vs. age",cex.lab=2,cex.main=2,cex.axis=2)
abline(h=0,col="red")
dev.off()


###bin residual plot by predicted value
cat("bin residual plot by predicted value\n")
pdf(paste(mydir,"/","binned_residual_plot.",description2,".pdf",sep=""),height=10,width=10)
y_predict_sorted=sort(y_predict)
classic_residuals_sorted=classic_residuals[names(y_predict_sorted)]
y_predict_binned=rep(NA,20)
classic_residuals_binned=rep(NA,20)
increment=round(length(y_predict)/20)
c=1
for(i in 1:20){
y_predict_binned[i]=mean(y_predict_sorted[c:min((c+increment),length(y_predict))])
classic_residuals_binned[i]=mean(classic_residuals_sorted[c:min((c+increment),length(y_predict))])
c=c+increment
}
plot(y_predict_binned,classic_residuals_binned,ylab="estimated residual, binned",xlab="expected count, binned",main="Binned estimated residuals vs. expected counts",cex.lab=2,cex.main=2,cex.axis=2)
abline(h=0,col="red")
dev.off()

###bin residual plot by predicted value with 95% interval
cat("bin residual plot by predicted value with 95% interval\n")
pdf(paste(mydir,"/","binned_residual_plot.",description2,".with_95pc_interval.pdf",sep=""),height=10,width=10)
y_predict_sorted=sort(y_predict)
classic_residuals_sorted=classic_residuals[names(y_predict_sorted)]
y_predict_binned=rep(NA,20)
classic_residuals_binned=rep(NA,20)
increment=round(length(y_predict)/20)
c=1
for(i in 1:20){
y_predict_binned[i]=mean(y_predict_sorted[c:min((c+increment),length(y_predict))])
classic_residuals_binned[i]=mean(classic_residuals_sorted[c:min((c+increment),length(y_predict))])
c=c+increment
}
plot(y_predict_binned,classic_residuals_binned,ylab="estimated residual, binned",xlab="expected count, binned",main="Binned estimated residuals vs. expected counts",cex.lab=2,cex.main=2,cex.axis=2,
     ylim=range(classic_residuals_sorted),xlim=range(y_predict),col="white")
abline(h=0,col="red")
c=1
for(i in 1:20){
    my.repictions=y_predict_sorted[c:min((c+increment),length(y_predict))]
    myresiduals=classic_residuals_sorted[c:min((c+increment),length(y_predict))]
    segments(y_predict_binned[i],quantile(myresiduals,0.025),y_predict_binned[i],quantile(myresiduals,0.975))
    segments(quantile(my.repictions,0.025),classic_residuals_binned[i],quantile(my.repictions,0.975),classic_residuals_binned[i])
    c=c+increment
}
dev.off()
names(Age)=1:length(Age)

####bin residual plot by age
cat("bin residual plot by age\n")
pdf(paste(mydir,"/","binned_residual_plot_by_age_with_50pc_interval.",description2,".pdf",sep=""),height=10,width=10)
Age_sorted=sort(Age)
classic_residuals_sorted=classic_residuals[names(Age_sorted)]
Age_binned=rep(NA,20)
classic_residuals_binned=rep(NA,20)
increment=round(length(y_predict)/20)
c=1
for(i in 1:20){
Age_binned[i]=mean(Age_sorted[c:min((c+increment),length(y_predict))])
classic_residuals_binned[i]=mean(classic_residuals_sorted[c:min((c+increment),length(y_predict))])
c=c+increment
}
plot(Age_binned,classic_residuals_binned,ylab="estimated residual, binned",xlab="Age, binned",main="Binned estimated residuals vs. Age",cex.lab=2,cex.main=2,cex.axis=2,
     ylim=c(-10,10),xlim=range(Age),col="white")
abline(h=0,col="red")
c=1
range.25th.to.75th = c()
for(i in 1:20){
    my.repictions=Age_sorted[c:min((c+increment),length(y_predict))]
    myresiduals=classic_residuals_sorted[c:min((c+increment),length(y_predict))]
    range.25th.to.75th =c(range.25th.to.75th ,quantile(myresiduals,0.75)-quantile(myresiduals,0.25))
    segments(Age_binned[i],quantile(myresiduals,0.25),Age_binned[i],quantile(myresiduals,0.75))
    segments(quantile(my.repictions,0.25),classic_residuals_binned[i],quantile(my.repictions,0.75),classic_residuals_binned[i])
    c=c+increment
}
plot(Age_binned,range.25th.to.75th,ylab="interquartile range of residuals",xlab="Age,binned",main="Interquartile range of residuals vs. Age",cex.lab=2,cex.main=2,cex.axis=2)
dev.off()


####plot residuals against predicted values
cat("plot residuals against predicted values\n")
pdf(paste(mydir,"/","residual_plot.",description2,".pdf",sep=""),height=10,width=10)
plot(y_predict,classic_residuals,ylab="estimated residual",xlab="expected count",main="Estimated residuals vs. expected counts",cex.lab=2,cex.main=2,cex.axis=2)
abline(h=0,col="red")
dev.off()

###Distribution of residuals - density and qq plot
cat("Distribution of residuals - density and qq plot\n")
pdf(paste(mydir,"/","distribution_of_residuals.",description2,".pdf",sep=""),height=5,width=10)
par(mfrow=c(1,2))
#density
plot(density(classic_residuals),xlab="estimated residuals",main="Distribution of estimated residuals")
#qqplot
my.sd=sd(classic_residuals)
myquantiles=qnorm(seq(0,1,l=length(classic_residuals)),mean=0,sd=my.sd)
qqplot(myquantiles,classic_residuals,xlab="Expected",ylab="Observed",main="QQ plot of residuals")
abline(a=0,b=1,col="grey")
dev.off()

####plot residuals against predicted values, by cohort
cat("plot residuals against predicted values, by cohort\n")
if(my.model<3){
    pdf(paste(mydir,"/","residual_plot.by_cohort.",description2,".pdf",sep=""),height=9,width=9)
    par(mfrow=c(3,3))
} else {
    pdf(paste(mydir,"/","residual_plot.by_cohort.",description2,".pdf",sep=""),height=12,width=21)
    par(mfrow=c(4,7))
}
for(c in 1:length(unique(cohort))){
plot(y_predict[cohort==c],classic_residuals[cohort==c],ylab="estimated residual",xlab="expected count",main=paste0(cohort.codes[cohort.codes$Code==c,"Cohort"],", n=",sum(cohort==c)),cex.main=1,cex.lab=2,cex.axis=2)
abline(h=0,col="red")
}
dev.off()


###Distribution of residuals - density and qq plot, by cohort
cat("Distribution of residuals - density and qq plot, by cohort\n")
if(my.model<3){
    pdf(paste(mydir,"/","distribution_of_residuals.by_cohort.",description2,".pdf",sep=""),height=5*ceiling(max(cohort)/2),width=20)
    par(mfrow=c(ceiling(max(cohort)/2),4))
} else {
    pdf(paste(mydir,"/","distribution_of_residuals.by_cohort.",description2,".pdf",sep=""),height=5*ceiling(max(cohort)/2),width=20)
    par(mfrow=c(ceiling(max(cohort)/2),4))
}
for(c in 1:length(unique(cohort))){
    plot(density(classic_residuals[cohort==c]),xlab="estimated residuals",main=paste0(cohort.codes[cohort.codes$Code==c,"Cohort"],", n=",sum(cohort==c)),cex.main=1,cex.lab=2,cex.axis=2)
    qqnorm(classic_residuals[cohort==c],xlab="Expected",ylab="Observed",main=paste0("QQ plot of residuals for ",cohort.codes[cohort.codes$Code==c,"Cohort"]),cex.main=1,cex.lab=2,cex.axis=2,cex=2)
    qqline(classic_residuals[cohort==c],col="grey")
}
dev.off()



}
