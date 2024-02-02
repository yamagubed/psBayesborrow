data {
  int<lower=0> nCT_o;
  int<lower=0> nCT_c;
  int<lower=0> nCC_o;
  int<lower=0> nCC_c;
  int<lower=0> p;
  vector[nCT_o] yCT_o;
  vector[nCT_c] yCT_c;
  vector[nCC_o] yCC_o;
  vector[nCC_c] yCC_c;
  row_vector[p] xCT_o[nCT_o];
  row_vector[p] xCT_c[nCT_c];
  row_vector[p] xCC_o[nCC_o];
  row_vector[p] xCC_c[nCC_c];
}
parameters {
  real theta;
  real gammaCC;
  vector[p] beta;
  real<lower=0> alpha;
}
model {
  for (i in 1:nCT_o)
    yCT_o[i] ~ weibull(alpha,exp(-(gammaCC+theta+xCT_o[i]*beta)/alpha));
  for (i in 1:nCC_o)
    yCC_o[i] ~ weibull(alpha,exp(-(gammaCC      +xCC_o[i]*beta)/alpha));

  if(nCT_c>0){
    for (i in 1:nCT_c)
      target += weibull_lccdf(yCT_c[i]|alpha,exp(-(gammaCC+theta+xCT_c[i]*beta)/alpha));
  }
  if(nCC_c>0){
    for (i in 1:nCC_c)
      target += weibull_lccdf(yCC_c[i]|alpha,exp(-(gammaCC      +xCC_c[i]*beta)/alpha));
  }
}
