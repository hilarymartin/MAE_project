data {
    int<lower=0> J; // number of children
    int<lower=0> I; // number of families
    int<lower=0> C; // number of cohorts
    real<lower=0> y[J]; // number of recombinations
    int<lower=0> family[J]; // indicates which family is with each child
    real<lower=0> Age[J]; //maternal age at birth
    int cohort[J];//which cohort (0=QTR,1=NTR)

   }
  parameters {
    real<lower=0> a0[I]; // random effect of families
    real beta_global; //mean of the distribution from which the age coefficients are drawn
    real sigmasq_global; //variance of the distribution from which the age coefficients are drawn
    vector[C] beta_Age;
    real<lower=0> tausq; //variance of number of scrossovers
    //parameters for distribution of a0
    real mu_m; //mean for intercept
    real<lower=0> sigmasq_m; //variance for intercept
  }
  transformed parameters {
    real<lower=0> a0_bychild[J];
    vector[J] mu;

    for(j in 1:J){
       a0_bychild[j]<-a0[family[j]];
       mu[j] <- a0_bychild[j]+ beta_Age[cohort[j]]*Age[j];
    }
  }
  model {
  //Priors
    beta_global ~ normal(0,1);
    sigmasq_global ~ inv_gamma(3,0.5);

    beta_Age ~ normal(beta_global,sqrt(sigmasq_global)); //coefficient for age

    tausq ~ inv_gamma(2,70);   //variance in number of crossovers

    mu_m ~ normal(38,sqrt(3));//pick something sensible - this is the mean of a0; I chose the intercept from lm(mat.counts$Count~mat.counts$Maternal.age.at.birth)
    sigmasq_m ~ inv_gamma(5,10);//pick something sensible - this is the variance of a0

   for(i in 1:I){ //random effects of the families
      a0[i] ~ normal(mu_m,sqrt(sigmasq_m));
   }
   y ~ normal(mu,sqrt(tausq)) ;
}

