// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// not needed
functions {

  real half_cauchy_custom_lp(real hc_prior) {

    real a1;
    real a2;

    a1 ~ inv_gamma(0.5, 1/hc_prior);
    a2 ~ inv_gamma(0.5, 1/a1);
    return (a2);
  }

}

// data
data {
  int<lower=0> nY;
  int<lower=0> nZ;
  int<lower=0> nG;
  real y[nY];
  matrix[nY, nZ] x;
  vector[nG] gSize;
  int gInd[nZ];
  vector[2] pr_sigma;
}

// The parameters accepted by the model.
parameters {
  vector[nZ] theta;
  real<lower=0> sigma_2;

  // shrinkage parameters
  real<lower=0> tau; // global
  real<lower=0> delta[nG]; // group
  real<lower=0> lambda[nZ]; // local
}

transformed parameters {

  vector[nZ] tr_Sigma;

  for (i in 1:nZ) {
    tr_Sigma[i] = sigma_2 * square(lambda[i]) * square(delta[gInd[i]]);
  };

}

// The model to be estimated.
model {

  // global
  tau ~ cauchy(0, 1);

  // group
  for (i in 1:nG) {
    delta[i] ~ exponential(1);
  };

  // local
  for (i in 1:nZ) {
    lambda[i] ~ exponential(1);
  };


  sigma_2 ~ inv_gamma(pr_sigma[1], pr_sigma[2]);

  theta ~ normal(0, tr_Sigma);

  y ~ normal(x * theta, sigma_2);

}


generated quantities {
  real yrep[nY];
  for (i in 1:nY) {
    yrep[i] = normal_rng(x[i] * theta, sigma_2);
  };
}
