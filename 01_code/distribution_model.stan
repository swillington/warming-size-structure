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
  int<lower=1>              N;      // number of sets
  int<lower=1>              R;      // number of locations
  int<lower=1>              J;      // number of size classes
  int<lower=1>              X;      // age-classes
  real                  lnEps;      // ln(epsilon)
  real<lower=0.001,upper=1> S[X];   // probability alive
  real<lower=0>          Lmid[J];   // length boundaries
  real<lower=0>             L[J+1]; // length boundaries
  int<lower=1,upper=R>      r[N];   // location index
  int<lower=0,upper=1>    mpa[N];   // mpa status: 0 = no, 1 = yes
  real                     zT[N];   // z-transformed sea-surface temperatures
  int<lower=1>             nT[N];   // total number of fish per set
  int<lower=0>              y[N,J]; // fish counts [set, sizeclass]
}

parameters {
  real<lower= 1.0,upper=50.00> phi;        // multinomial overdispersion parameter
  real<lower=30.00,upper=55.00> ELInf;     // expected maximum size 
  real<lower= 0.50,upper= 2.50> MpK;       // typical M / K across sites and years
  real<lower=-0.50,upper= 0.50> beta_z;    // effect of temperature on M / K
  real<lower=-0.50,upper= 0.50> beta_mpa;  // effect of mpa on M / K
  real<lower= 0.00,upper=25.00> L50;       // length when 50% observed                     
  real<lower= 1.00,upper=20.00> DL95;      // length when 95% observed= L50+DL95                   
  real<lower= 0.05,upper= 0.30> cv;        // coefficient of variation in LInf
  real<lower=0.1,upper=1.00>    sigma_Loc; // between location variation                 
  real<lower=0.1,upper=1.00>    sigma_N;   // between year:location variation                 
  real MPK_Loc_RE[R];   // random deviations around overall mean for all locations
  real MPK_N_RE[N];   // random deviations around overall mean for all locations
}

transformed parameters{
}

model {
  vector[J] g;     // probability in length bin j
  vector[J] alpha; // multinomial parameter for length bin j
  vector[J] s;     // probability observe
  vector[X] Lbar;  // expected length at age
  
  array[J, X] real z; // z-transformed length

  real sum1;
  real sum2;
  real alpha0;
  real MpKi; // M / K for specific i
  
  // set priors 
  phi      ~ exponential(0.1);
  ELInf    ~ normal(40.0, 10.0);
  MpK      ~ normal( 1.5,  0.5);
  beta_z   ~ normal( 0.0,  0.2);
  beta_mpa ~ normal( 0.0,  0.2);
  L50      ~ normal(15.0,  5.0);
  DL95     ~ normal(10.0,  5.0);
  cv       ~ normal( 0.1,  0.05);
  sigma_Loc ~ normal(0.0,  0.5);
  sigma_N   ~ normal(0.0,  0.5);

  MPK_Loc_RE ~ normal(0.0, sigma_Loc); // random effect deviations (vectorised) 
  MPK_N_RE   ~ normal(0.0, sigma_N);   
  
  alpha0 = 1.0 / phi; // alpha sum
  
  for (j in 1:J) { // probability observed when in length bin j
    s[j] = 1.0 / (1 + exp(-(log(19)/DL95)*(Lmid[j]-L50)));
  }

  for (i in 1:N) { // for each set of counts
    // set specific M /K
    MpKi = exp(log(MpK) + beta_z*zT[i] + beta_mpa*mpa[i] + 
      MPK_Loc_RE[r[i]] +  MPK_N_RE[i]); 
  
    for (x in 1:X) { // expected length at age
      Lbar[x] = ELInf*(1.0 - exp(lnEps*x/(X*MpKi)));
    }
    
    for (x in 1:X) {
      for (j in 1:J) { // z-transformed length-at-age
        z[j,x] = (L[j+1] - Lbar[x])/(cv*Lbar[x]);
      }
    }
    
    // calculate g_j for this observation i
    
    sum2 = 0.0; // first length box
    for (x in 1:X) {
      sum2 += S[x]*normal_cdf(z[1,x], 0, 1);
    }
    g[1] = s[1]*sum2;
    sum1 = g[1];
    
    for (j in 2:(J-1)) { // intermediate length boxes
      sum2 = 0.0; 
      for (x in 1:X) {
        sum2 += S[x]*(normal_cdf(z[j,x], 0, 1) - normal_cdf(z[j-1,x], 0, 1));
      }
      g[j] = s[j]*sum2;
      sum1 += g[j];
    }
  
    sum2 = 0.0; // last length box
    for (x in 1:X) {
      sum2 += S[x]*(1.0 - normal_cdf(z[J,x], 0, 1));
    }
    g[J] = s[J]*sum2;
    sum1 += g[J];
    
    // calculate alpha
    for (j in 1:J) {
       g[j] = g[j] / sum1; // normalise
       alpha[j] = g[j] / alpha0; // scale probabilities to multinomial params
    }
    
    // calculate log-likelihood
    sum1 = lgamma(alpha0) + lgamma(nT[i] + 1) - lgamma(nT[i] + alpha0);
    for (j in 1:J) { // for each sub-period
      sum1 += lgamma(y[i,j] + alpha[j]) - 
        lgamma(alpha[j]) - lgamma(y[i,j] + 1); 
    }
    
    target += sum1; // add log-likelihood term associated with set of obs
  }
}

generated quantities {
}
