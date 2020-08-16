// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

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

  real<lower=0> delta_prior[nG]; // additional prior
  real<lower=0> lambda_prior[nZ];

}

transformed parameters {

  vector[nZ] tr_Sigma;

  for (i in 1:nZ) {
    tr_Sigma[i] = sigma_2 * lambda[i] * delta[gInd[i]];
  };

}

// The model to be estimated.
model {

  // global
  tau ~ cauchy(0, 1);

  // group
  for (i in 1:nG) {
    delta_prior[i] ~ cauchy(0, 1);
    delta[i] ~ cauchy(0, delta_prior[i]);
  };

  // local
  for (i in 1:nZ) {
    lambda_prior[i] ~ cauchy(0, 1);
    lambda[i] ~ cauchy(0, lambda_prior[i]);
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
