
data {
    int<lower=0> J; // number of children
    int<lower=0> I; // number of families
    int<lower=0> C; // number of cohorts
//    real<lower=0> y[J]; // number of recombinations
    vector[J] y; // number of recombinations
    int<lower=0> family[J]; // indicates which family is with each child
    real<lower=0> Age[J]; //maternal age at birth
    int cohort[J];//which cohort (0=QTR,1=NTR)
    int cohort_by_family[I];//which cohort (0=QTR,1=NTR)	
    real mean_alpha_prior;	
    real<lower=0> sigmasq_mu_m_prior;		
    int<lower=0> sigmasq_m_alpha;
    int<lower=0> sigmasq_m_beta;
    matrix[J,J] P;
    matrix[J,J] Tw;
    matrix[J,J] D;
   }
  parameters {
    real<lower=0> a0[I]; // random effect of families
    real beta_Age;
    real<lower=0> tausq; //variance of number of scrossovers
    vector[C] mu_m; //mean for intercept
    real<lower=0> sigmasq_m[C] ; //variance for intercept
    real<lower=0,upper=1> lambda;
  }
  transformed parameters {
    vector[J] mu;
    cov_matrix[J] Sigma;
    for(j in 1:J){
       mu[j] <- mu_m[cohort[j]]+ beta_Age*Age[j];
//build covariance matrix
	for(i in 1:J){
       	      Sigma[i,j] <- (D[i,j]+lambda*Tw[i,j])*tausq+P[i,j]*sigmasq_m[cohort[j]];
	}      
   }

//     print("testSigma=",testSigma);
  }
  model {
  //Priors
    beta_Age ~ normal(0,1); //coefficient for age
    tausq ~ inv_gamma(2,70);   //variance in number of crossovers

//prior on lambda is uniform on (0,1) so no need to specify

    mu_m ~ normal(mean_alpha_prior,sigmasq_mu_m_prior);
    for(c in 1:C){    
       sigmasq_m[c] ~ inv_gamma(sigmasq_m_alpha,sigmasq_m_beta);//this is the variance of a0
    }
   for(i in 1:I){ //random effects of the families
      a0[i] ~ normal(mu_m[cohort_by_family[i]],sqrt(sigmasq_m[cohort_by_family[i]]));
   }
//make this multivariate normal
   y ~ multi_normal(mu,Sigma) ;
}

