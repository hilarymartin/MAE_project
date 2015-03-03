data {
    int<lower=0> J; // number of children
    int<lower=0> I; // number of families
    int<lower=0> C; // number of cohorts
    real<lower=0> y[J]; // number of recombinations
    int<lower=0> family[J]; // indicates which family is with each child
    real<lower=0> Age[J]; //maternal age at birth
    int cohort[J];//which cohort (0=QTR,1=NTR)
    int cohort_by_family[I];//which cohort (0=QTR,1=NTR)	
//    int<lower=0> df_alpha;
      real mean_alpha_prior;	
      real<lower=0> sigmasq_mu_m_prior;		
      int<lower=0> sigmasq_m_alpha;
      int<lower=0> sigmasq_m_beta;
   }
  parameters {
    real<lower=0> a0[I]; // random effect of families
    real beta_Age;
    real<lower=0> tausq[C]; //variance of number of scrossovers
    vector[C] mu_m; //mean for intercept
    real<lower=0> sigmasq_m[C] ; //variance for intercept
  }
  transformed parameters {
    real<lower=0> a0_bychild[J];
    vector[J] mu;
    for(j in 1:J){
       a0_bychild[j]<-a0[family[j]];
       mu[j] <- a0_bychild[j]+ beta_Age*Age[j];
    }
  }
  model {
  //Priors
    beta_Age ~ normal(0,1); //coefficient for age

    mu_m ~ normal(mean_alpha_prior,sigmasq_mu_m_prior);
    for(c in 1:C){    
       sigmasq_m[c] ~ inv_gamma(sigmasq_m_alpha,sigmasq_m_beta);//this is the variance of a0
       tausq[c] ~ inv_gamma(2,70);   //variance in number of crossovers
    }
   for(i in 1:I){ //random effects of the families
      a0[i] ~ normal(mu_m[cohort_by_family[i]],sqrt(sigmasq_m[cohort_by_family[i]]));
   }
   for(j in 1:J){
     y[j] ~ normal(mu[j],sqrt(tausq[cohort[j]])) ;
     }
}

