  data { 
    int<lower=0> J; // number of children 
    int<lower=0> I; // number of mothers 
    int<lower=0> C; // number of cohorts
    int y[J]; // number of recombinations
    int family[J]; // indicates which mother is with each child
    vector[J] Age; //maternal age at birth
    int cohort[J];//which cohort (0=QTR,1=NTR)
    real alpha_parameters[C]; //alpha parameter for beta distributions that are the prior for p(observed crossover|true crossovers)
    real beta_parameters[C]; //beta parameter for beta distributions that are the prior for p(observed crossover|true crossovers)
    real mean_alpha_prior;
    int cohort_by_family[I];//which cohort (0=QTR,1=NTR)

  }
  parameters {
    vector[I] exp_a0; // random effect of mothers 
    vector[C] beta_Age;  //coefficient for age
    real beta_global; //mean of the distribution from which the age coefficients are drawn
    real sigmasq_global; //variance of the distribution from which the age coefficients are drawn
    real<lower=0,upper=1> p_by_cohort[C]; //probability of calling crossovers in each cohort (where cohort is actually a cohort-family type pair)
    real<lower=0,upper=1> inv_omega; //overdispersion - this is distributed Uniform[0,1] so omega is constrained on the range (1, infinity)
    vector[C] mu_m; //mean for intercept
    real<lower=0> sigmasq_m[C]; //variance for intercept
  }
  transformed parameters {
      real<lower=1> omega;
      vector[J] exp_a0_bychild;
      vector[J] mu;
      real p[J];
      omega <- 1 / inv_omega;
      for(j in 1:J){
        exp_a0_bychild[j] <- exp_a0[family[j]];
        p[j] <- p_by_cohort[cohort[j]];
        mu[j] <-  exp(0.1*(exp_a0_bychild[j]+beta_Age[cohort[j]]*Age[j])) / (omega - 1); ///needed to keep in exp(0.1*...) to ensure that mu > 0, since shape parameter of neg binomial must be >0
      }
  }
  model {
    beta_global ~ normal(0,0.05);///smaller variance since the age effect is now multiplicative: baseline * exp(0.1* beta_Age * age)
    sigmasq_global ~ inv_gamma(3,0.1);///smaller variance since the age effect is now multiplicative: baseline * exp(0.1* beta_Age * age)

    beta_Age ~ normal(beta_global,sqrt(sigmasq_global)); //coefficient for age

    for(c in 1:C){
       p_by_cohort[c] ~ beta(alpha_parameters[c],beta_parameters[c]);
    }

//set this differently for males and females
    mu_m ~ normal(mean_alpha_prior,sqrt(6)); //smaller since baseline is now exp(0.1*mu_m)


    for(c in 1:C){
          sigmasq_m[c] ~ inv_gamma(5,5);
    }

    inv_omega ~ uniform(0,1);


//have cohort-specific mean and variance
    for(i in 1:I){ //random effects of the mothers
      exp_a0[i] ~ normal(mu_m[cohort_by_family[i]],sqrt(sigmasq_m[cohort_by_family[i]]));
    }

    for( j in 1:J){
      y[j] ~ neg_binomial(mu[j], 1.0 / (p[j] * (omega - 1)));
    }

  }