//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> nY;
  int<lower=0> nK;
  int<lower=0> nG;
  int<lower=0> nMaxGroup;
  real y[nY];
  matrix[nY, nK] x;
  vector[nG] gSize;
  int gInd[nK];
  vector[2] pr_sigma;
  real<lower=0> pi_0;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[nK] beta;
  real<lower=0> sigma;
  real<lower=0> tau[nG];
  real<lower=0> lambda;
}

transformed parameters {

  vector[nK] tr_Sigma;
  vector[nK] beta_ss;


  for (i in 1:nK) {
    tr_Sigma[i] = sigma * square(tau[gInd[i]]);
  };

  beta_ss = (1 - pi_0) * beta + pi_0 * 0;

}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  pi_0 ~ beta(0.5, 0.5);
  lambda ~ exponential(1);
  // print(lambda);
  tau ~ gamma((gSize + 1)/2, lambda/2);

  sigma ~ inv_gamma(pr_sigma[1], pr_sigma[2]);

  beta ~ normal(0, tr_Sigma);

  y ~ normal(x * beta_ss, sigma);

}


generated quantities {
  real yrep[nY];
  for (i in 1:nY) {
    yrep[i] = normal_rng(x[i] * beta, sigma);
  };
}
