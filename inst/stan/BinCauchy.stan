data {
  int<lower=0> nCT;
  int<lower=0> nCC;
  int<lower=0> nEC;
  int<lower=0> p;
  int<lower=0,upper=1> yCT[nCT];
  int<lower=0,upper=1> yCC[nCC];
  int<lower=0,upper=1> yEC[nEC];
  row_vector[p] xCT[nCT];
  row_vector[p] xCC[nCC];
  row_vector[p] xEC[nEC];
  real<lower=0> scale;
}
parameters {
  real theta;
  real gammaCC;
  real gammaEC;
  real<lower=0> tau;
  vector[p] beta;
}
model {
  tau ~ cauchy(0,scale);
  gammaCC ~ normal(gammaEC,tau);

  for (i in 1:nCT)
    yCT[i] ~ bernoulli(inv_logit(gammaCC+theta+xCT[i]*beta));
  for (i in 1:nCC)
    yCC[i] ~ bernoulli(inv_logit(gammaCC      +xCC[i]*beta));
  for (i in 1:nEC)
    yEC[i] ~ bernoulli(inv_logit(gammaEC      +xEC[i]*beta));
}
