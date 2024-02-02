data {
  int<lower=0> nCT;
  int<lower=0> nCC;
  int<lower=0> p;
  vector[nCT] yCT;
  vector[nCC] yCC;
  row_vector[p] xCT[nCT];
  row_vector[p] xCC[nCC];
}
parameters {
  real theta;
  real gammaCC;
  vector[p] beta;
  real<lower=0> sigma;
}
model {
  for (i in 1:nCT)
    yCT[i] ~ normal(gammaCC+theta+xCT[i]*beta,sigma);
  for (i in 1:nCC)
    yCC[i] ~ normal(gammaCC      +xCC[i]*beta,sigma);
}
