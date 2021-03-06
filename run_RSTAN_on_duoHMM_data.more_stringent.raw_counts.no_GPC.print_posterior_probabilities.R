library(rstan)
argv <- commandArgs(trailingOnly = TRUE)

cohort.specific=0
#Model 1
if(argv[1] ==43 | argv[1]==44){
    if(argv[1]==43){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_41_100.sigmasq_m_IG_2_40.RData"))
        mysim<-extract(model1.6.mat,permuted=T)
        model="Model1.mat"
    } else {
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_27_64.sigmasq_m_IG_2_15.RData"))
        mysim<-extract(model1.6.pat,permuted=T)
                model="Model1.pat"
    }
}

#Model 1 with cohort-specific tausq                 
if(argv[1] ==3 | argv[1]==4){
    if(argv[1]==3){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.6.2.cohort_specific_tausq.all_chains.mat.including_alpha.mu_m_N_41_100.sigmasq_m_IG_2_40.RData"))
        mysim<-extract(model1.6.mat,permuted=T)
        model="Model1.cohort_specific_tausq.mat"
    } else {
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.6.2.cohort_specific_tausq.all_chains.pat.including_alpha.mu_m_N_27_64.sigmasq_m_IG_2_15.RData"))
        mysim<-extract(model1.6.pat,permuted=T)
        model="Model1.cohort_specific_tausq.pat"
    }
}


#Model 1.2
if(argv[1] ==29 | argv[1]==30){
    cohort.specific=1
    if(argv[1]==29){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.6.all_chains.mat.including_alpha.mu_m_N_41_100.sigmasq_m_IG_2_40.RData"))
        mysim<-extract(model1.6.mat,permuted=T)
        model="Model1.2.mat"
    } else {
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.6.all_chains.pat.including_alpha.mu_m_N_27_64.sigmasq_m_IG_2_15.RData"))
        mysim<-extract(model1.6.pat,permuted=T)
        model="Model1.2.pat"
    }
}

#Model 2
if(argv[1] ==45 | argv[1]==46){
    if(argv[1]==45){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.5.4.all_chains.mat.including_alpha.sigmasq_m_IG_5_5.RData"))
        mysim<-extract(model1.5.mat,permuted=T)
        model="Model2.mat"
    } else {
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.5.4.all_chains.pat.including_alpha.sigmasq_m_IG_5_5.RData"))
        mysim<-extract(model1.5.pat,permuted=T)
        model="Model2.pat"
    }
}


    #Model 2.2
if(argv[1] ==33 | argv[1]==34){
        cohort.specific=1
    if(argv[1]==33){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.5.all_chains.maternal.mu_m_N_37_6.sigmasq_m_IG_5_5.including_alpha.RData"))
        mysim<-extract(model1.5.mat,permuted=T)
        model="Model2.2.mat"
    } else {
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.5.all_chains.paternal.mu_m_N_32_6.sigmasq_m_IG_5_5.including_alpha.RData"))
        mysim<-extract(model1.5.pat,permuted=T)
        model="Model2.2.pat"
    }
}

  
#Model 2 with cohort-specific omega
if(argv[1] ==5 | argv[1]==6){
    if(argv[1]==5){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.5.4.cohort_specific_omega.all_chains.mat.including_alpha.sigmasq_m_IG_5_5.RData"))
        mysim<-extract(model1.5.mat,permuted=T)
        model="Model2.cohort_specific_omega.mat"
    } else{
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.modell1.5.4.cohort_specific_omega.all_chains.pat.including_alpha.sigmasq_m_IG_5_5.RData"))
        mysim<-extract(model1.5.pat,permuted=T)
        model="Model2.cohort_specific_omega.pat"
    }
}



#Model 3
if(argv[1] %in% c(27,28)){
    if(argv[1] ==27){
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.model3.2.adjusted_priors.all_chains.mat.including_alpha.RData")
        mysim<-extract(model3.2.mat,permuted=T)
        model="Model3.mat"
    } else if(argv[1] ==28){
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.model3.2.adjusted_priors.all_chains.pat.including_alpha.RData")
        mysim<-extract(model3.2.pat,permuted=T)
        model="Model3.pat"
    }
}

    
#Model 3.2
if(argv[1]==1|argv[1]==2){
        cohort.specific=1
    if(argv[1] ==1){
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.model3.beta_age_by_cohort.adjusted_priors.all_chains.mat.including_alpha.RData")
        mysim<-extract(model3.mat,permuted=T)
        model="Model3.2.mat"
    } else if(argv[1] ==2){
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC.raw_counts/RSTAN.model3.beta_age_by_cohort.adjusted_priors.all_chains.pat.including_alpha.RData")
        mysim<-extract(model3.pat,permuted=T)
        model="Model3.2.pat"
    }
}


if(cohort.specific==0){
    cat(model,"\t",sum(mysim$beta_Age > 0)/length(mysim$beta_Age))
} else {
    cat(model,"\t",colSums(mysim$beta_Age > 0)/nrow(mysim$beta_Age))
}
