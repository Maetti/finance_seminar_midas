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
  int N;
  int K;
  int N_samples;
  matrix[N, K] x_test;
  matrix[N_samples, K] theta;
  vector[N_samples] sigma_2;

}
parameters {
}
model {
}

generated quantities {
  matrix[N_samples, N] y_test;
  for(n in 1:N) {
    for(i in 1:N_samples) {
      y_test[i, n] = normal_rng(x_test[n]*theta[i]', sigma_2[i]);
    }
  }
}
