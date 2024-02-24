data {
  int<lower=0> N;          // Number of observations
  array[N] int<lower=0, upper=1> y;  // Binary outcome variable
  matrix[N, 2] X;          // Matrix of predictors
}

parameters {
  vector[2] beta;          // Coefficients for predictors
}

model {
  // Prior distribution for coefficients
  beta ~ normal(0, 1);
  
  // Likelihood function
  for (i in 1:N) {
    y[i] ~ bernoulli_logit(dot_product(X[i], beta));
  }
}
