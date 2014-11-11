  data { 
    int<lower=0> J; // number of children 
    int<lower=0> I; // number of families
    int<lower=0> C; // number of cohorts
    int y[J]; // number of recombinations
    int family[J]; // indicates which family is with each child
    vector[J] Age; //maternal age at birth
    int cohort[J];//which cohort
    real alpha_parameters[C]; //alpha parameter for beta distributions that are the prior for p(observed crossover|true crossovers)
    real beta_parameters[C]; //beta parameter for beta distributions that are the prior for p(observed crossover|true crossovers)
  }
  parameters {
    vector[I] exp_a0; // random effect of familys 
    vector[C] beta_Age;  //coefficient for age
    real<lower=0,upper=1> p_by_cohort[C]; //probability of calling crossovers in each cohort (where cohort is actually a cohort-family type pair)
    real<lower=0,upper=1> inv_omega; //overdispersion - this is distributed Uniform[0,1] so omega is constrained on the range (1, infinity)
    real mu_m; //mean  of a0
    real<lower=0> sigmasq_m; //variance of a0
  }
  transformed parameters {
      real<lower=1> omega;
      vector[J] exp_a0_bychild;
      real<lower=0> mu[J];
      real p[J];
      omega <- 1 / inv_omega;
      for(j in 1:J){
        exp_a0_bychild[j] <- exp_a0[family[j]];
        p[j] <- p_by_cohort[cohort[j]];
//        mu[j] <-  (exp_a0_bychild[j]+beta_Age[cohort[j]]*Age[j]) / (omega - 1);
//        mu[j] <-  (exp_a0_bychild[j]+beta_Age[cohort[j]]*Age[j]) / 0.4;
//        mu[j] <-  (exp_a0_bychild[j]+ beta_Age[cohort[j]] * Age[j]) / (omega - 1); // beta_Age[c] ~ normal(1,0.5) T[0,];  exp_a0[i] ~ normal(38,2) T[0,];
//        mu[j] <-  exp_a0_bychild[j] / (omega - 1);
        mu[j] <-  exp_a0_bychild[j];
//        mu[j] <-  38;
      }
  }
  model {
//    beta_Age ~ normal(0,1); //coefficient for age
//    beta_Age ~ normal(0,0.5); //coefficient for age
    for(c in 1:C){
       p_by_cohort[c] ~ beta(alpha_parameters[c],beta_parameters[c]);
       beta_Age[c] ~ normal(1,0.5) T[0,]; //coefficient for age	
    }       	      
//    mu_m ~ normal(38,sqrt(3));//pick something sensible - this is the mean of a0; I chose the intercept from lm(mat.counts$Count~mat.counts$Maternal.age.at.birth)
//    sigmasq_m ~ inv_gamma(5,10);//pick something sensible - this is the variance of a0
      inv_omega ~ uniform(0,1);

    for(i in 1:I){ //random effects of the families
//      exp_a0[i] ~ normal(mu_m,sqrt(sigmasq_m));
 //     exp_a0[i] ~ normal(mu_m,2);
     	exp_a0[i] ~ normal(38,2) T[0,];
    }
//    for (j in 1:J){
//        y_true[j] ~ neg_binomial(mu[j], 1.0 / (omega - 1));
//        y_true[j] ~ neg_binomial(40.1, 1.1);    
//   }
//    y ~ binomial(y_true,p)
//      y ~ neg_binomial(mu,(omega - 2) / p);
	for( j in 1:J){
      y[j] ~ neg_binomial(mu[j], 1.0 / (p[j] * (omega - 1)));
}
  }

