data{
  int<lower=1> N;//number of annual samples (time-series length)
  real TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower = 0> log_a0;// average productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this
  vector[N] z_yr; //annual deviations in slope
  
 //variance components  
  real<lower=0> sigma_e;
  real<lower=0> gp_tau;
  real<lower=0> gp_rho;

}
transformed parameters{
	vector[N] log_a;
	real b;
	
	matrix[N,N] L_Tmat;
	matrix[N,N] Tmat;
	Tmat = cov_exp_quad(TT, gp_tau, gp_rho);
	for(n in 1:N) Tmat[n, n] += 1e-9;
	L_Tmat = cholesky_decompose(Tmat);

	log_a =log_a0+L_Tmat*z_yr;
	
	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
}
model{
  //priors
  log_a0 ~ normal(0,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(-12,3); //average per capita capacity parameter
  
  //variance terms
  sigma_e ~ gamma(2,5);
  gp_tau ~ gamma(2,5);
  gp_rho ~ inv_gamma(5, 5);
  
  //scaled effects of time
  z_yr ~ std_normal();

  for(i in 1:N) R_S[i] ~ normal(log_a[i] - b*S[i], sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  vector[N] y_rep;
  for (t in 1:N){ log_lik[t] = normal_lpdf(R_S[t]|log_a[t] - S[t]*b, sigma_e);
				  y_rep[t] = normal_rng(log_a[t] - S[t]*b, sigma_e);
  }
  }
   

