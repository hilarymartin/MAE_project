data {
    int<lower=0> J; // number of children
    int<lower=0> I; // number of families
    int<lower=0> C; // number of cohorts
    int y[J]; // number of recombinations
    int<lower=0> family[J]; // indicates which family is with each child
    vector[J] Age; //maternal age at birth
    int cohort_by_family[I];//which cohort (0=QTR,1=NTR)
    int cohort[J];//which cohort
    
   }
  parameters {
    vector[I] exp_a0; // random effect of mothers
    real beta_Age;
    real<lower=0,upper=1> inv_omega; //overdispersion - this is distributed Uniform[0,1] so omega is constrained on the range (1, infinity)
    vector[C] mu_m; //mean for intercept
//    real<lower=0> inv_sigmasq_m[C]; //variance for intercept
//    real log_sigmasq_m[C]; //variance for intercept
      real sigma_m[C];	
  }
  transformed parameters {
    vector[J] exp_a0_bychild;
    real<lower=0> sigmasq_m[C]; //variance for intercept
    vector[J] mu;
    real<lower=1> omega;
    omega <- 1 / inv_omega;
    for(c in 1:C){
//    	  sigmasq_m[c] <- inv_sigmasq_m[c];
    //sigmasq_m[c] <- exp(log_sigmasq_m[c]);	
    sigmasq_m[c] <- sigma_m[c]*sigma_m[c];	
    }
    for(j in 1:J){
        exp_a0_bychild[j] <- exp_a0[family[j]];
        mu[j] <-  exp(0.1*(exp_a0_bychild[j]+beta_Age*Age[j])) / (omega - 1); ///needed to keep in exp(0.1*...) to ensure that mu > 0, since shape parameter of neg binomial must be >0
    }
  }
  model {
  //Priors
    for(i in 1:I){ //random effects of the mothers
      exp_a0[i] ~ normal(mu_m[cohort_by_family[i]],sqrt(sigmasq_m[cohort_by_family[i]]));
    }

    for( j in 1:J){
      y[j] ~ neg_binomial(mu[j], 1.0 / (omega - 1));
    }

}

