data {
  int<lower=0> N;
  real y[N];
}
transformed data {
  real sum_log_y;
  sum_log_y <- 0.0;
  for(i in 1:N)
    sum_log_y <- sum_log_y + log(y[i]);
}
parameters {
  real<lower=0> mu;
  real<lower=0> lambda;
}
model {
  lp__ <- lp__ + 0.5 * N * log(lambda / (2.0 * pi())) - 1.5 * sum_log_y;
  for(i in 1:N)
    lp__ <- lp__ - lambda * pow((y[i] - mu) / mu, 2) / (2.0 * y[i]);
}