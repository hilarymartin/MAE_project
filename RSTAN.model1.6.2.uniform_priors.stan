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
   }
  parameters {
    real<lower=0> a0[I]; // random effect of families
    real beta_Age;
    vector[C] mu_m; //mean for intercept
//    real<lower=0> inv_tausq;
//    real log_tausq;
	real tau;
//    real<lower=0> inv_sigmasq_m[C];
//    real log_sigmasq_m[C];
    real sigma_m[C];
  }
  transformed parameters {
    real<lower=0> tausq; //variance of number of scrossovers
    real<lower=0> sigmasq_m[C] ; //variance for intercept	
    real<lower=0> a0_bychild[J];
    vector[J] mu;
      tausq <- tau*tau;
    for(c in 1:C){
    sigmasq_m[c] <- sigma_m[c]*sigma_m[c];		
    }   

    for(j in 1:J){
       a0_bychild[j]<-a0[family[j]];
       mu[j] <- a0_bychild[j]+ beta_Age*Age[j];
    }
  }
  model {
   for(i in 1:I){ //random effects of the families
      a0[i] ~ normal(mu_m[cohort_by_family[i]],sqrt(sigmasq_m[cohort_by_family[i]]));
   }
   y ~ normal(mu,sqrt(tausq)) ;
}

