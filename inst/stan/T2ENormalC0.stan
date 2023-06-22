data {
  int<lower=1> nCT_o;
  int<lower=1> nCT_c;
  int<lower=1> nCC_o;
  int<lower=1> nEC_o;
  int<lower=1> nEC_c;
  int<lower=1> p;
  vector[nCT_o] yCT_o;
  vector[nCT_c] yCT_c;
  vector[nCC_o] yCC_o;
  vector[nEC_o] yEC_o;
  vector[nEC_c] yEC_c;
  row_vector[p] xCT_o[nCT_o];
  row_vector[p] xCT_c[nCT_c];
  row_vector[p] xCC_o[nCC_o];
  row_vector[p] xEC_o[nEC_o];
  row_vector[p] xEC_c[nEC_c];
  real<lower=0> scale;
}
parameters {
  real theta;
  real gammaCC;
  real gammaEC;
  real<lower=0> tau;
  vector[p] beta;
  real<lower=0> alpha;
}
model {
  tau ~ normal(0,scale);
  gammaCC ~ normal(gammaEC,tau);

  for (i in 1:nCT_o)
    yCT_o[i] ~ weibull(alpha,exp(-(gammaCC+theta+xCT_o[i]*beta)/alpha));
  for (i in 1:nCC_o)
    yCC_o[i] ~ weibull(alpha,exp(-(gammaCC      +xCC_o[i]*beta)/alpha));
  for (i in 1:nEC_o)
    yEC_o[i] ~ weibull(alpha,exp(-(gammaEC      +xEC_o[i]*beta)/alpha));

  for (i in 1:nCT_c)
    target += weibull_lccdf(yCT_c[i]|alpha,exp(-(gammaCC+theta+xCT_c[i]*beta)/alpha));
  for (i in 1:nEC_c)
    target += weibull_lccdf(yEC_c[i]|alpha,exp(-(gammaEC      +xEC_c[i]*beta)/alpha));
}
