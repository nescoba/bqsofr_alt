data{
    int<lower=1> n;
    vector[n] x;
    vector[n] y;
    real<lower=0.0, upper=1.0> p_0;
    real<lower=0.0> c;
    real beta_0;
    real<lower=0.0> B_0;
    real<lower=0.0> n_0;
    real d_0;
    real delta_0;
    real<lower=0.0> D_0;
    real<lower=0.0> shape_1;
    real<lower=0.0> shape_2;
}
transformed data {
    real L;
    real U;
    L = -sqrt(1 - p_0);
    U = sqrt(p_0);
}
parameters {
   real beta_par;
   real gamma_par;
   real<lower=0.0> sigma_par;
   real<lower=0.0> delta_par;
   vector[n] h_raw;
   vector[n] nu_rvar;
}
transformed parameters {
   real<lower=0.0, upper=1.0> p;
   real indic;
   real A;
   real B;
   real C;
   real gamma_trans;
   real xi_3;
   real <lower=0.0> sigma_par_2;
   vector[n] h_rvar;
   indic = 0;
   if (gamma_par < 0){
    indic = 1;
   }
   p = indic + (p_0 - indic)/(2*Phi(-fabs(gamma_par))*exp(gamma_par^2/2));
   A = (1 - 2*p)/(p*(1 - p));
   B = 2/(p*(1 - p));
   C = (1 - indic - p)^(-1);
   gamma_trans = (gamma_par - L )/(U - L);
   xi_3 = c + exp(delta_par);
   sigma_par_2 = sigma_par^2;
   h_rvar = fabs(h_raw);
}
model{
    vector[n] mean;
    vector[n] sd;
    h_raw ~ normal(rep_vector(0.0, n), sigma_par);
    nu_rvar ~ exponential(rep_vector(sigma_par, n));
    for(i in 1:n){
        mean[i] = x[i]*beta_par + A*nu_rvar[i] + C*fabs(gamma_par)*h_rvar[i];
        sd[i] = sqrt(sigma_par*B*nu_rvar[i]);
        target += -2*log(sd[i]);
        if (y[i] == 1){
            target += log(Phi(-mean[i]/sd[i]));
        } else if (y[i] == 2){
            target += log(Phi((c - mean[i])/sd[i]) - Phi(-mean[i]/sd[i]));
        } else if (y[i] == 3){
            target += log(Phi((xi_3 - mean[i])/sd[i]) - Phi(c-mean[i]/sd[i]));
        } else {
            target += log(1 - Phi((xi_3 - mean[i])/sd[i]));
        }
    }
    beta_par ~ normal(beta_0, B_0);
    sigma_par_2 ~ inv_gamma(n_0/2, d_0/2);
    gamma_trans ~ beta(shape_1, shape_2);
    delta_par ~ normal(delta_0, D_0);
}

