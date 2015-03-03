data {
    int<lower=0> J; // number of children
    int<lower=0> I; // number of families
    int<lower=0> C; // number of cohorts
    int y[J]; // number of recombinations
    int<lower=0> family[J]; // indicates which family is with each child
    vector[J] Age; //maternal age at birth
    int cohort_by_family[I];//which cohort (0=QTR,1=NTR)
    int cohort[J];//which cohort
    real mean_alpha_prior;	
    int<lower=0> sigmasq_m_alpha;
    int<lower=0> sigmasq_m_beta;	
   }
  parameters {
    vector[I] exp_a0; // random effect of mothers
    real beta_Age;
    real<lower=0,upper=1> inv_omega[C]; //overdispersion - this is distributed Uniform[0,1] so omega is constrained on the range (1, infinity)
    vector[C] mu_m; //mean for intercept
    real<lower=0> sigmasq_m[C]; //variance for intercept
  }
  transformed parameters {
    real<lower=1> omega[C];
    vector[J] exp_a0_bychild;
    vector[J] mu;
    for(c in 1:C){
    omega[c] <- 1 / inv_omega[c];
}
    for(j in 1:J){
        exp_a0_bychild[j] <- exp_a0[family[j]];
        mu[j] <-  exp(0.1*(exp_a0_bychild[j]+beta_Age*Age[j])) / (omega[cohort[j]] - 1); ///needed to keep in exp(0.1*...) to ensure that mu > 0, since shape parameter of neg binomial must be >0
    }
  }
  model {
  //Priors
    beta_Age ~ normal(0,0.05); //coefficient for age

    mu_m ~ normal(mean_alpha_prior,sqrt(6)); //smaller since baseline is now exp(0.1*mu_m)
    for(c in 1:C){    
	  sigmasq_m[c] ~ inv_gamma(sigmasq_m_alpha,sigmasq_m_beta);
    inv_omega[c] ~ uniform(0,1);
    }


    for(i in 1:I){ //random effects of the mothers
      exp_a0[i] ~ normal(mu_m[cohort_by_family[i]],sqrt(sigmasq_m[cohort_by_family[i]]));
    }

    for( j in 1:J){
      y[j] ~ neg_binomial(mu[j], 1.0 / (omega[cohort[j]] - 1));
    }

}

