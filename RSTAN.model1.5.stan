data {
    int<lower=0> J; // number of children
    int<lower=0> I; // number of families
    int<lower=0> C; // number of cohorts
    int y[J]; // number of recombinations
    int<lower=0> family[J]; // indicates which family is with each child
 	vector[J] Age; //maternal age at birth
    int cohort[J];//which cohort
   }
  parameters {
    vector[I] exp_a0; // random effect of mothers
    real beta_global; //mean of the distribution from which the age coefficients are drawn
    real sigmasq_global; //variance of the distribution from which the age coefficients are drawn
    vector[C] beta_Age;
    vector[C-1] beta_Cohort;
    real<lower=0,upper=1> inv_omega; //overdispersion - this is distributed Uniform[0,1] so omega is constrained on the range (1, infinity)
    real mu_m; //mean for intercept
    real<lower=0> sigmasq_m; //variance for intercept
  }
  transformed parameters {
    real<lower=1> omega;
    vector[J] my_beta_Cohort;
    vector[J] exp_a0_bychild;
    vector[J] mu;
    omega <- 1 / inv_omega;
    for(j in 1:J){
        exp_a0_bychild[j] <- exp_a0[family[j]];
        if(cohort[j]==1){
             my_beta_Cohort[j] <- 0;
        } else {
             my_beta_Cohort[j] <- beta_Cohort[cohort[j]-1];
        }
        mu[j] <-  exp(0.1*(exp_a0_bychild[j]+beta_Age[cohort[j]]*Age[j] + my_beta_Cohort[cohort[j]] )) / (omega - 1); ///needed to keep in exp(0.1*...) to ensure that mu > 0, since shape parameter of neg binomial must be >0
    }
  }
  model {
  //Priors
    beta_global ~ normal(0,0.05);///smaller variance since the age effect is now multiplicative: baseline * exp(0.1* beta_Age * age)
    sigmasq_global ~ inv_gamma(3,0.1);///smaller variance since the age effect is now multiplicative: baseline * exp(0.1* beta_Age * age)

    beta_Age ~ normal(beta_global,sqrt(sigmasq_global)); //coefficient for age

    mu_m ~ normal(36,sqrt(6)); //smaller since baseline is now exp(0.1*mu_m)
    sigmasq_m ~ inv_gamma(5,5);

    inv_omega ~ uniform(0,1);

    beta_Cohort ~ normal(0,1); //coefficient for cohort

    for(i in 1:I){ //random effects of the mothers
      exp_a0[i] ~ normal(mu_m,sqrt(sigmasq_m));
    }

    for( j in 1:J){
      y[j] ~ neg_binomial(mu[j], 1.0 / (omega - 1));
    }

}

