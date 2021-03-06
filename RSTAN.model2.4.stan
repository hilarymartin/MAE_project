  data {
    int<lower=0> J; // number of children
    int<lower=0> C; // number of cohorts
    int cohort[J];
    real y[J]; // number of recombinations
    real Age[J]; //maternal age at birth
    real<lower=0> sigmasq_beta_nkids;	 
    int<lower=0> tausq_alpha_prior;
    int<lower=0> tausq_beta_prior;
    real<lower=0> nkids[J];
   }
  parameters {
    real beta_Age;
    real beta_nkids;
    real<lower=0> tausq; //variance of number of scrossovers
  }
  transformed parameters {
    vector[J] mu;
    for(j in 1:J){
       mu[j] <- beta_Age*Age[j] + beta_nkids*nkids[j];
    }
  }
    model {	
  //Priors
    beta_Age ~ normal(0,1); //coefficient for age
    beta_nkids ~ normal(0,sqrt(sigmasq_beta_nkids));
    tausq ~ inv_gamma(tausq_alpha_prior,tausq_beta_prior);   //variance in number of crossovers

    y ~ normal(mu,sqrt(tausq)) ;
}

