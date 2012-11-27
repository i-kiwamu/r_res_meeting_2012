data {
  int<lower=0> N;
  real<lower=0> Nile[N];
}
parameters {
  real theta[N];
  real F;
  real G;
  //real theta_zero;
  real<lower=0> sigma_y;
  real<lower=0> sigma_theta;
}
model {
  theta[1] ~ normal(Nile[1], sigma_theta);
  Nile[1] ~ normal(F * theta[1], sigma_y);
  for(i in 2:N){
    theta[i] ~ normal(G * theta[i-1], sigma_theta);
    Nile[i] ~ normal(F * theta[i], sigma_y);
  }
}