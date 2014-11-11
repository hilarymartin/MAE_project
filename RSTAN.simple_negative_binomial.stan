  data { 
    int<lower=0> J; // number of children 
    int<lower=0> I; // number of mothers 
    int<lower=0> C; // number of cohorts
    int y[J]; // number of recombinations
    int family[J]; // indicates which mother is with each child
    vector[J] Age; //maternal age at birth
    int cohort[J];//which cohort (0=QTR,1=NTR)
  }
  parameters {
    vector[I] exp_a0; // random effect of mothers 
    vector[C] beta_Age;  //coefficient for age
    vector[C-1] beta_Cohort;
    real<lower=0,upper=1> inv_omega; //overdispersion - this is distributed Uniform[0,1] so omega is constrained on the range (1, infinity)
    real mu_m; //mean  of a0
    real<lower=0> sigmasq_m; //variance of a0
  }
  transformed parameters {
      real<lower=1> omega;
      vector[J] exp_a0_bychild;
      vector[J] mu;
      vector[J] my_beta_Cohort;
      omega <- 1 / inv_omega;
      for(j in 1:J){
          exp_a0_bychild[j] <- exp_a0[family[j]];
         if(cohort[j]==1){
              my_beta_Cohort[j] <- 0;
         } else {
              my_beta_Cohort[j] <- beta_Cohort[cohort[j]-1];
         }
//        mu[j] <-  (exp_a0_bychild[j]+beta_Age[cohort[j]]*Age[j]+my_beta_Cohort[j]) / (omega - 1);
        mu[j] <-  exp(0.1*(exp_a0_bychild[j]+beta_Age[cohort[j]]*Age[j]+my_beta_Cohort[j])) / (omega - 1);
      }
  }
  model {
    beta_Age ~ normal(0,1); //coefficient for age
    beta_Cohort ~ normal(0,4);

    mu_m ~ normal(38,sqrt(3));//pick something sensible - this is the mean of a0; I chose the intercept from lm(mat.counts$Count~mat.counts$Maternal.age.at.birth)
    sigmasq_m ~ inv_gamma(5,10);//pick something sensible - this is the variance of a0
    inv_omega ~ uniform(0,1);

    for(i in 1:I){ //random effects of the mothers
      exp_a0[i] ~ normal(mu_m,sqrt(sigmasq_m));
    }

    y ~ neg_binomial(mu, 1.0 / (omega - 1));
  }