data {
  int<lower=1> nCT;
  int<lower=1> nCC;
  int<lower=1> nEC;
  int<lower=1> p;
  int<lower=0,upper=1> yCT[nCT];
  int<lower=0,upper=1> yCC[nCC];
  int<lower=0,upper=1> yEC[nEC];
  row_vector[p] xCT[nCT];
  row_vector[p] xCC[nCC];
  row_vector[p] xEC[nEC];
}
parameters {
  real theta;
  real gammaCC;
  vector[p] beta;
}
model {
  for (i in 1:nCT)
    yCT[i] ~ bernoulli(inv_logit(gammaCC+theta+xCT[i]*beta));
  for (i in 1:nCC)
    yCC[i] ~ bernoulli(inv_logit(gammaCC      +xCC[i]*beta));
  for (i in 1:nEC)
    yEC[i] ~ bernoulli(inv_logit(gammaCC      +xEC[i]*beta));
}
