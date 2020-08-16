// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> G;
  int<lower=0> g1;
  int<lower=0> g2;
  int<lower=0> g3;
  int<lower=0> g4;
  vector[N] y;
  matrix[N, K] x;
  vector<lower=0>[G] g_full;
  vector[2] a1;
  vector[2] a2;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[g1] theta_1;
  vector[g2] theta_2;
  vector[g3] theta_3;
  vector[g4] theta_4;
  vector<lower=0>[G] tau;
  real<lower=0> sigma;
  vector<lower=0>[G] lambda;
}


transformed parameters {
  vector[K] theta_all;
  theta_all = append_row(theta_1, append_row(theta_2, append_row(theta_3, theta_4)));        // school treatment effects
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y ~ normal(x * theta_all, sigma);

  for (i in 1:G)  {
    tau[i] ~ gamma((g_full[i] + 1)/2, pow(lambda[i], 2)/2);
    lambda[i] ~ inv_gamma(a2[1], a2[1]);
  }

  theta_1 ~ normal(0, sigma * tau[1]);
  theta_2 ~ normal(0, sigma * tau[2]);
  theta_3 ~ normal(0, sigma * tau[3]);
  theta_4 ~ normal(0, sigma * tau[4]);

  sigma ~ inv_gamma(a1[1], a1[1]);

}


generated quantities {

  real yrep[N];
  for (i in 1:N) {
    yrep[i] = normal_rng(x[i] * theta_all, sigma);
  };

  // vector[N] yhat;                  // linear predictor
  // real<lower=0> rss;             // residual sum of squares
  // real<lower=0> totalss;         // total SS
  // real Rsq;                      // Rsq
  //
  // yhat = x * theta_all;
  // rss = dot_self(y-yhat);
  // totalss = dot_self(y-mean(y));
  // Rsq = 1 - rss/totalss;

}

