#---- This will fit the model using STAN
#-------------------------------------------------------------------------------
#-- Implementing Model in Rstan
library(brms)
library(rstan)
library(splines)

#-------------------------------------------------------------------------------
#--- Import Data
#-------------------------------------------------------------------------------
source("~/Library/CloudStorage/OneDrive-IndianaUniversity/Join_folder_Annie_DrZoh/NHANES Dataset/2013-2014/BQReg_Nhanes2012_13Application.R")
load("/Users/rszoh/Library/CloudStorage/OneDrive-IndianaUniversity/Join_folder_Annie_DrZoh/NHANES Dataset/2013-2014/FinalDataV2.RData")

#-------------------------------------------------------------------------------
#--- Set up the data
#-------------------------------------------------------------------------------
pn = 10
a <- seq(0,1, length = dim(Warray)[2])
bs2 <- bs(a, df = pn, intercept = T)


#Xt <- apply(Warray, c(1,2), mean)
Xt =  as.matrix(FUI(model="gaussian", smooth=T, Warray ,silent = FALSE)) #} #- (abs(range(W)[1]) +1)
#Xt <- as.matrix(FUI(model="gaussian",smooth=T, Warray,silent = FALSE))
Wt <- (Xt%*%bs2)/length(a)


#---- Transform the data
# df_sample$AgeYRTru <- df_sample$AgeYR
# #df_sample$AgeYR <- scale(df_sample$AgeYR, center = T, scale = T)
# df_sampleOld <- df_sample
# df_sample$AgeYR <- df_sample$AgeYR/100
# df_sample$MnASTP <- df_sample$MnASTP*10
# df_sample$MnMVPA <-  df_sample$MnMVPA/100


df_sampleWt <- data.frame(df_sample, W = Wt)[1:200,]



#********************************************************************************
# trial code (Cases we need to look at)
# 
#-------------------------------------------------------------------------------
#case1 (Fit the model in brms ignoring MVPA and ASTP)
#case2 (Fit the model in STAN with smoothing prior and ignoring MVPA and ASTP)
#case3 (Fit the model in brms with MVPA and ASTP)
#Case4 (Fit the model in STAN with MVPA and ASTP and smoothing prior)

#-------------------------------------------------------------------------------
#---- Data set prep. 
#-------------------------------------------------------------------------------
tau0 = .1

form <- paste0("Y ~ ", paste(c("Gender","Race","HealthCondt2","AgeYR", colnames(df_sampleWt)[grepl("W",colnames(df_sampleWt))]), collapse = " + "))
Dt0 <- make_standata(bf(as.formula(form), quantile = tau0), data = df_sampleWt[c(1:nrow(df_sampleWt)),], family = asym_laplace())
#, stanvars = stanvars2, prior = bprior,chains = 2,iter = 500, control = list(adapt_delta = 0.97), data2 = list(M = rep(0, 29), V = SigbPrio))


fitCase1  <- list()
fitCase1b <- list()
fitCase2  <- list()
fitCase2b <- list()
fitCase3  <- list()
fitCase4  <- list()
fitCase5  <- list()
tauVal <- c(.1,.5,.9)
Niter = 100
Nchain = 2


for(l in 1:length(tauVal)){ 
  
  tau0 <- tauVal[l]
  Bd = GamBnd(tau0)[1:2]
  
  Dt0$pn = pn
  Dt0$Xold <- Dt0$X
  Dt0$X <- Dt0$Xold[,!grepl("W.[0-9]", colnames(Dt0$Xold))]
  Dt0$Wn <- Dt0$Xold[,grepl("W.[0-9]", colnames(Dt0$Xold))]; dim(Dt0$Wn)
  Dt0$M <- numeric(Dt0$pn)
  Dt0$V <- as.positive.definite(as.inverse(P.mat(pn) + diag(rep(.01, pn))))*.1
  Dt0$K <- ncol(Dt0$X)
  Bd <- GamBnd(tau0)[1:2]
  Dt0$Bd <- Bd*.99
  Dt0$tau0 <- tau0
  
  varX <- c("MnASTP", "MnMVPA")
  
  #-------------------------------------------------------------------------------
  # Case 1 (no MPVA and ASTP)
  #-------------------------------------------------------------------------------
  #@' This function is the stan_funs2 for the GAL distribution
  
  Bd = GamBnd(tau0)[1:2]
  
  {
    GAL2 <- custom_family(
      "GAL2", dpars = c("mu", "sigma","gam", "tau"), links=c("identity","log","identity","identity"),
      lb=c(NA,0, Bd[1]*.9,0), ub=c(NA,NA, Bd[2]*.9,1), type="real") #, vars = "vint1[n]"
    
    
    stan_funs2 <- "
  /*
  A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + log1m_exp(fabs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2)))); 
  gam = (gamU - gamL) * ligam + gamL;
  real gam = (gamU - gamL) * ligam + gamL;
  real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau, real gamL, real gamU){
  real GAL2_rng(real mu, real sigma, real ligam, real tau, real gamL, real gamU){
   */
     /* helper function for asym_laplace_lpdf
  * Args:
    *   y: the response value
  *   tau: quantile parameter in (0, 1)
  */
    real rho_quantile(real y, real tau) {
      if (y < 0) {
        return y * (tau - 1);
      } else {
        return y * tau;
      }
    }
  /* asymmetric laplace log-PDF for a single response
  * Args:
    *   y: the response value
  *   mu: location parameter
  *   sigma: positive scale parameter
  *   tau: quantile parameter in (0, 1)
  * Returns:
    *   a scalar to be added to the log posterior
  */
    real asym_laplace_lpdf(real y, real mu, real sigma, real tau) {
      return log(tau * (1 - tau)) -
        log(sigma) -
        rho_quantile((y - mu) / sigma, tau);
    }
  /* asymmetric laplace log-CDF for a single quantile
  * Args:
    *   y: a quantile
  *   mu: location parameter
  *   sigma: positive scale parameter
  *   tau: quantile parameter in (0, 1)
  * Returns:
    *   a scalar to be added to the log posterior
  */
    real asym_laplace_lcdf(real y, real mu, real sigma, real tau) {
      if (y < mu) {
        return log(tau) + (1 - tau) * (y - mu) / sigma;
      } else {
        return log1m((1 - tau) * exp(-tau * (y - mu) / sigma));
      }
    }
  /* asymmetric laplace log-CCDF for a single quantile
  * Args:
    *   y: a quantile
  *   mu: location parameter
  *   sigma: positive scale parameter
  *   tau: quantile parameter in (0, 1)
  * Returns:
    *   a scalar to be added to the log posterior
  */
    real asym_laplace_lccdf(real y, real mu, real sigma, real tau) {
      if (y < mu) {
        return log1m(tau * exp((1 - tau) * (y - mu) / sigma));
      } else {
        return log1m(tau) - tau * (y - mu) / sigma;
      }
    }
   
   real GAL2_lpdf(real y, real mu, real sigma, real gam, real tau){
   
   real p_pos;
   real p_neg;
   real a3;
   real a2;
   real p;
   real est;
   real A;
   real B;
   real Res = 0;
   //real gam = ligam;
    p = 1 * (gam < 0) + (tau - 1 * (gam < 0))/(2*Phi(-fabs(gam)) * exp(.5 * pow(gam, 2)));
    p_pos = p - 1 * (gam > 0);
    p_neg = p -  1* (gam < 0);  
    est = (y - mu) / sigma; 
    
    if(fabs(gam) > 0){
    a3 = p_pos * (est / fabs(gam));
    a2 = fabs(gam) * (p_neg / p_pos);
    
    
    if(est/gam > 0){
      A =  0.5 * pow(gam, 2) * pow(p_neg/p_pos, 2) - est * p_neg + log_diff_exp(log(Phi_approx(a2-a3)), log(Phi_approx(a2)) ); 
      B =  0.5 * pow(gam, 2) - p_pos * est + log(Phi_approx(-fabs(gam) + a3));
      Res = log(2*p*(1-p)) - log(sigma) +  log_sum_exp(A, B);
    }else{
      Res =  log(2*p*(1-p)) - log(sigma) - p_pos * est + 0.5 * pow(gam, 2) + log(Phi_approx(-fabs(gam) ));
    }
    }else{
    Res = asym_laplace_lpdf( y | mu, sigma, tau); 
    }
     
    return Res;
   }
  
  real GAL2_rng(real mu, real sigma, real ligam, real tau){
  
     real A;
     real B;
     real C;
     real p;
     real hi;
     real nui;
     real mui=0;
     real Up = uniform_rng(.5, 1.0);
     
     real gam = ligam;
     p = (gam < 0) + (tau - (gam < 0))/(2*Phi_approx(-fabs(gam))*exp(.5*pow(gam, 2)));
     A = (1 - 2*p)/(p - pow(p,2));
     B = 2/(p - pow(p,2));
     C = 1/((gam > 0) - p);
     
      hi = sigma * inv_Phi(Up);
     nui = sigma * exponential_rng(1);
     mui += mu + A * nui + C * fabs(gam) * hi;
  
     return normal_rng(mui, sqrt(sigma*B*nui));
  }
  "
  }  
  
  #--- Now define all of these here
  form <- paste0("Y ~ ", paste(c("Gender","Race","HealthCondt2","AgeYR", colnames(df_sampleWt)[grepl("W",colnames(df_sampleWt))]), collapse = " + "))
  
  stanvars2 <- stanvar(scode = stan_funs2, block = "functions")
  fitCase1[[l]] <- brms::brm(bf(as.formula(form), tau = tau0), data = df_sampleWt, family = GAL2, stanvars = stanvars2,
                             chains = Nchain,iter = Niter, control = list(adapt_delta = 0.99), cores = 2, seed = 1123) # init = 0.1,
  
  
  #-------------------------------------------------------------------------------
  # Case 1b (with MPVA and ASTP) - Wt
  #-------------------------------------------------------------------------------
  #@' This function is the stan_funs2 for the GAL distribution
  
  Bd = GamBnd(tau0)[1:2]
  
  {
    GAL2 <- custom_family(
      "GAL2", dpars = c("mu", "sigma","gam", "tau"), links=c("identity","log","identity","identity"),
      lb=c(NA,0, Bd[1]*.9,0), ub=c(NA,NA, Bd[2]*.9,1), type="real") #, vars = "vint1[n]"
    
    
    stan_funs2 <- "
  /*
  A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + log1m_exp(fabs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2)))); 
  gam = (gamU - gamL) * ligam + gamL;
  real gam = (gamU - gamL) * ligam + gamL;
  real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau, real gamL, real gamU){
  real GAL2_rng(real mu, real sigma, real ligam, real tau, real gamL, real gamU){
   */
     /* helper function for asym_laplace_lpdf
  * Args:
    *   y: the response value
  *   tau: quantile parameter in (0, 1)
  */
    real rho_quantile(real y, real tau) {
      if (y < 0) {
        return y * (tau - 1);
      } else {
        return y * tau;
      }
    }
  /* asymmetric laplace log-PDF for a single response
  * Args:
    *   y: the response value
  *   mu: location parameter
  *   sigma: positive scale parameter
  *   tau: quantile parameter in (0, 1)
  * Returns:
    *   a scalar to be added to the log posterior
  */
    real asym_laplace_lpdf(real y, real mu, real sigma, real tau) {
      return log(tau * (1 - tau)) -
        log(sigma) -
        rho_quantile((y - mu) / sigma, tau);
    }
  /* asymmetric laplace log-CDF for a single quantile
  * Args:
    *   y: a quantile
  *   mu: location parameter
  *   sigma: positive scale parameter
  *   tau: quantile parameter in (0, 1)
  * Returns:
    *   a scalar to be added to the log posterior
  */
    real asym_laplace_lcdf(real y, real mu, real sigma, real tau) {
      if (y < mu) {
        return log(tau) + (1 - tau) * (y - mu) / sigma;
      } else {
        return log1m((1 - tau) * exp(-tau * (y - mu) / sigma));
      }
    }
  /* asymmetric laplace log-CCDF for a single quantile
  * Args:
    *   y: a quantile
  *   mu: location parameter
  *   sigma: positive scale parameter
  *   tau: quantile parameter in (0, 1)
  * Returns:
    *   a scalar to be added to the log posterior
  */
    real asym_laplace_lccdf(real y, real mu, real sigma, real tau) {
      if (y < mu) {
        return log1m(tau * exp((1 - tau) * (y - mu) / sigma));
      } else {
        return log1m(tau) - tau * (y - mu) / sigma;
      }
    }
   
   real GAL2_lpdf(real y, real mu, real sigma, real gam, real tau){
   
   real p_pos;
   real p_neg;
   real a3;
   real a2;
   real p;
   real est;
   real A;
   real B;
   real Res = 0;
   //real gam = ligam;
    p = 1 * (gam < 0) + (tau - 1 * (gam < 0))/(2*Phi(-fabs(gam)) * exp(.5 * pow(gam, 2)));
    p_pos = p - 1 * (gam > 0);
    p_neg = p -  1* (gam < 0);  
    est = (y - mu) / sigma; 
    
    if(fabs(gam) > 0){
    a3 = p_pos * (est / fabs(gam));
    a2 = fabs(gam) * (p_neg / p_pos);
    
    
    if(est/gam > 0){
      A =  0.5 * pow(gam, 2) * pow(p_neg/p_pos, 2) - est * p_neg + log_diff_exp(log(Phi_approx(a2-a3)), log(Phi_approx(a2)) ); 
      B =  0.5 * pow(gam, 2) - p_pos * est + log(Phi_approx(-fabs(gam) + a3));
      Res = log(2*p*(1-p)) - log(sigma) +  log_sum_exp(A, B);
    }else{
      Res =  log(2*p*(1-p)) - log(sigma) - p_pos * est + 0.5 * pow(gam, 2) + log(Phi_approx(-fabs(gam) ));
    }
    }else{
    Res = asym_laplace_lpdf( y | mu, sigma, tau); 
    }
     
    return Res;
   }
  
  real GAL2_rng(real mu, real sigma, real ligam, real tau){
  
     real A;
     real B;
     real C;
     real p;
     real hi;
     real nui;
     real mui=0;
     real Up = uniform_rng(.5, 1.0);
     
     real gam = ligam;
     p = (gam < 0) + (tau - (gam < 0))/(2*Phi_approx(-fabs(gam))*exp(.5*pow(gam, 2)));
     A = (1 - 2*p)/(p - pow(p,2));
     B = 2/(p - pow(p,2));
     C = 1/((gam > 0) - p);
     
      hi = sigma * inv_Phi(Up);
     nui = sigma * exponential_rng(1);
     mui += mu + A * nui + C * fabs(gam) * hi;
  
     return normal_rng(mui, sqrt(sigma*B*nui));
  }
  "
  }  
  
  #--- Now define all of these here
  formCase1b <- paste0("Y ~ ", paste(c("Gender","Race","HealthCondt2","AgeYR"), collapse = " + "))
  for(pa in varX){
    formCase1b <- paste0(formCase1b, " + s(", pa, ", bs='ps', k=11)")
  }
  
  stanvars2 <- stanvar(scode = stan_funs2, block = "functions")
  fitCase1b[[l]] <- brms::brm(bf(as.formula(formCase1b), tau = tau0), data = df_sampleWt, family = GAL2, stanvars = stanvars2,
                              chains = Nchain,iter = Niter, control = list(adapt_delta = 0.99), cores = 2,seed = 1123) # init = 0.1,
  
  
  #-------------------------------------------------------------------------------
  # Case 2
  #-------------------------------------------------------------------------------
  {
    StandCodeVGalLin <- "
      // > fit_fakeqGal3a.5$model
      // generated with brms 2.18.0
      functions {
        
        /*
        A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + log1m_exp(fabs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2)))); 
        gam = (gamU - gamL) * ligam + gamL;
        real gam = (gamU - gamL) * ligam + gamL;
        real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau, real gamL, real gamU){
        real GAL2_rng(real mu, real sigma, real ligam, real tau, real gamL, real gamU){
         */
           /* helper function for asym_laplace_lpdf
        * Args:
          *   y: the response value
        *   tau: quantile parameter in (0, 1)
        */
          real rho_quantile(real y, real tau) {
            if (y < 0) {
              return y * (tau - 1);
            } else {
              return y * tau;
            }
          }
        /* asymmetric laplace log-PDF for a single response
        * Args:
          *   y: the response value
        *   mu: location parameter
        *   sigma: positive scale parameter
        *   tau: quantile parameter in (0, 1)
        * Returns:
          *   a scalar to be added to the log posterior
        */
          real asym_laplace_lpdf(real y, real mu, real sigma, real tau) {
            return log(tau * (1 - tau)) -
              log(sigma) -
              rho_quantile((y - mu) / sigma, tau);
          }
        /* asymmetric laplace log-CDF for a single quantile
        * Args:
          *   y: a quantile
        *   mu: location parameter
        *   sigma: positive scale parameter
        *   tau: quantile parameter in (0, 1)
        * Returns:
          *   a scalar to be added to the log posterior
        */
          real asym_laplace_lcdf(real y, real mu, real sigma, real tau) {
            if (y < mu) {
              return log(tau) + (1 - tau) * (y - mu) / sigma;
            } else {
              return log1m((1 - tau) * exp(-tau * (y - mu) / sigma));
            }
          }
        /* asymmetric laplace log-CCDF for a single quantile
        * Args:
          *   y: a quantile
        *   mu: location parameter
        *   sigma: positive scale parameter
        *   tau: quantile parameter in (0, 1)
        * Returns:
          *   a scalar to be added to the log posterior
        */
          real asym_laplace_lccdf(real y, real mu, real sigma, real tau) {
            if (y < mu) {
              return log1m(tau * exp((1 - tau) * (y - mu) / sigma));
            } else {
              return log1m(tau) - tau * (y - mu) / sigma;
            }
          }
         
         real GAL2_lpdf(real y, real mu, real sigma, real gam, real tau){
         
         real p_pos;
         real p_neg;
         real a3;
         real a2;
         real p;
         real est;
         real A;
         real B;
         real Res = 0;
         //real gam = ligam;
          p = 1 * (gam < 0) + (tau - 1 * (gam < 0))/(2*Phi(-fabs(gam))*exp(.5*pow(gam, 2)));
          p_pos = p -  1 * (gam > 0);
          p_neg = p -  1 * (gam < 0);  
          est = (y - mu) / sigma; 
          
          if(fabs(gam) > 0){
          a3 = p_pos * (est / fabs(gam));
          a2 = fabs(gam) * (p_neg / p_pos);
          
          
          if(est/gam > 0){
            A =  0.5*pow(gam, 2)*pow(p_neg/p_pos, 2) - est*p_neg + log_diff_exp(log(Phi_approx(a2-a3)), log(Phi_approx(a2)) ); 
            B =  0.5*pow(gam, 2) - p_pos*est + log(Phi_approx(-fabs(gam) + a3));
            Res = log(2*p*(1-p)) - log(sigma) +  log_sum_exp(A, B);
          }else{
            Res =  log(2*p*(1-p)) - log(sigma) - p_pos * est + 0.5 * pow(gam, 2) + log(Phi_approx(-fabs(gam) ));
          }
          }else{
          Res = asym_laplace_lpdf( y | mu, sigma, tau); 
          }
           
          return Res;
         }
        
        real GAL2_rng(real mu, real sigma, real gam, real tau){
        
           real A;
           real B;
           real C;
           real p;
           real hi;
           real nui;
           real mui=0;
           real Up = uniform_rng(.5, 1.0);
           
          // real gam = ligam;
           p = (gam < 0) + (tau - (gam < 0))/(2*Phi_approx(-fabs(gam))*exp(.5*pow(gam, 2)));
           A = (1 - 2*p)/(p - pow(p,2));
           B = 2/(p - pow(p,2));
           C = 1/((gam > 0) - p);
           
            hi = sigma * inv_Phi(Up);
           nui = sigma * exponential_rng(1);
           mui += mu + A * nui + C * fabs(gam) * hi;
        
           return normal_rng(mui, sqrt(sigma*B*nui));
        }
        
      }
      data {
        int<lower=1> pn; //Number of basis;
        int<lower=1> N;  // total number of observations
        vector[N] Y;  // response variable
        int<lower=1> K;  // number of population-level effects
        matrix[N, K] X;  // population-level design matrix
        int prior_only;  // should the likelihood be ignored?
        matrix[N,pn] Wn;
          vector[pn] M;
         matrix[pn, pn] V;
         vector[2] Bd;
         real<lower=0,upper=1> tau0;
      }
      transformed data {
        int Kc = K - 1;
        matrix[N, Kc] Xc;  // centered version of X without an intercept
        vector[Kc] means_X;  // column means of X before centering
        for (i in 2:K) {
          means_X[i - 1] = mean(X[, i]);
          Xc[, i - 1] = X[, i] - means_X[i - 1];
        }
      }
      parameters {
        vector[Kc] b;  // population-level effects
        real Intercept;  // temporary intercept for centered predictors
        vector[pn] bw;
        real<lower=0.01,upper=.99> psi;
         
        // standarized spline coefficients
        real<lower=0> sigma;  // dispersion parameter
        real<lower=Bd[1],upper=Bd[2]> gam;
        real<lower=0> teta;
      }
      transformed parameters {
         //bwa ar(1) stationary prior
         vector[pn] bwa;
        real tau = tau0;
        real lprior = 0;  // prior contributions to the log posterior
        // compute bwa
        bwa[1] = bw[1];
        for(i in 2:pn)
          bwa[i] = psi*bwa[i-1] + bw[i]*teta; 
          
        lprior += normal_lpdf(b | 0.0, 10);
        lprior += std_normal_lpdf(bw);
        lprior += uniform_lpdf(psi | 0.01, .99);
        lprior += uniform_lpdf(gam | Bd[1], Bd[2]);
        //lprior += cauchy_lpdf(teta| 0, 1) - 1* cauchy_lccdf(0| 0, 1);
        //lprior += gamma_lpdf(teta| 1, 1);
        //lprior += weibull_lpdf(teta| 1, 2.5);
        lprior += std_normal_lpdf(teta) - 1* log(0.5);
        lprior += student_t_lpdf(Intercept | 3, 28.7, 5.6);
        //lprior += normal_lpdf(bs | 0, 5);
        // lprior += student_t_lpdf(sds_1_1 | 3, 0, 5.6)
        //   - 1 * student_t_lccdf(0 | 3, 0, 5.6);
        // lprior += student_t_lpdf(sds_2_1 | 3, 0, 5.6)
        //   - 1 * student_t_lccdf(0 | 3, 0, 5.6);
        lprior += student_t_lpdf(sigma | 3, 0, 5.6)
          - 1 * student_t_lccdf(0 | 3, 0, 5.6);
      }
      model {
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term + Xs * bs
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b + Wn * bwa;
          for (n in 1:N) {
          //  target += GAL2_lpdf(Y[n] | mu[n], sigma, ligam, tau0);
              target += GAL2_lpdf(Y[n] | mu[n], sigma, gam, tau0);
          }
        }
        // priors including constants
        target += lprior;
        //target += std_normal_lpdf(zs_1_1);
        //target += std_normal_lpdf(zs_2_1);
      }
      generated quantities {
        // actual population-level intercept
        real b_Intercept = Intercept - dot_product(means_X, b);
        vector[N] mu = rep_vector(0.0, N);
        vector[N] log_lik;
          mu += Intercept + Xc * b + Wn * bwa;
        for(n in 1:N)
              log_lik[n] +=  GAL2_lpdf(Y[n] | mu[n], sigma, gam, tau0);

        //Vector[N] ASTP = Xs[1:N,1] * bs[1] + Zs_1_1 * s_1_1;
        //vector[N] MVPA = Xs[1:N,2] * bs[2] + Zs_2_1 * s_2_1;
      }
      "
  }
  df2 <- df_sampleWt
  
  form2 <- paste0("Y ~ ", paste(c("Gender","Race","HealthCondt2","AgeYR", colnames(df2)[grepl("W",colnames(df2))]), collapse = " + "))
  Dt2 <- make_standata(bf(as.formula(form2), quantile = tau0), data = df2[c(1:nrow(df2)),], family = asym_laplace())
  
  Dt2$pn = pn
  Dt2$Xold <- Dt2$X
  Dt2$X <- Dt2$Xold[,!grepl("W.[0-9]", colnames(Dt2$Xold))]
  Dt2$Wn <- Dt2$Xold[,grepl("W.[0-9]", colnames(Dt2$Xold))]; dim(Dt2$Wn)
  Dt2$M <- numeric(Dt2$pn)
  Dt2$V <- as.positive.definite(as.inverse(P.mat(pn) + diag(rep(.01, pn))))*.1
  Dt2$K <- ncol(Dt2$X)
  Bd <- GamBnd(tau0)[1:2]
  Dt2$Bd <- Bd*.9
  Dt2$tau0 <- tau0
  Dt2$X[,ncol(Dt2$X)] <- Dt2$X[,ncol(Dt2$X)]/100
  
  
  fitCase2[[l]] <- stan(model_code = StandCodeVGalLin,
                        data = Dt2, iter = Niter, chains = Nchain, verbose = TRUE,control = list(adapt_delta = 0.99,max_treedepth=11), cores = 2, seed = 1123)
  
  
  
  #-------------------------------------------------------------------------------
  # Case 2b  (with mixture distribution)
  #-------------------------------------------------------------------------------
  {
    StandCodeVGalLin <- "
      // > fit_fakeqGal3a.5$model
      // generated with brms 2.18.0
      functions {
        
        /*
        A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + log1m_exp(fabs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2)))); 
        gam = (gamU - gamL) * ligam + gamL;
        real gam = (gamU - gamL) * ligam + gamL;
        real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau, real gamL, real gamU){
        real GAL2_rng(real mu, real sigma, real ligam, real tau, real gamL, real gamU){
         */
           /* helper function for asym_laplace_lpdf
        * Args:
          *   y: the response value
        *   tau: quantile parameter in (0, 1)
        */
          real rho_quantile(real y, real tau) {
            if (y < 0) {
              return y * (tau - 1);
            } else {
              return y * tau;
            }
          }
        /* asymmetric laplace log-PDF for a single response
        * Args:
          *   y: the response value
        *   mu: location parameter
        *   sigma: positive scale parameter
        *   tau: quantile parameter in (0, 1)
        * Returns:
          *   a scalar to be added to the log posterior
        */
          real asym_laplace_lpdf(real y, real mu, real sigma, real tau) {
            return log(tau * (1 - tau)) -
              log(sigma) -
              rho_quantile((y - mu) / sigma, tau);
          }
        /* asymmetric laplace log-CDF for a single quantile
        * Args:
          *   y: a quantile
        *   mu: location parameter
        *   sigma: positive scale parameter
        *   tau: quantile parameter in (0, 1)
        * Returns:
          *   a scalar to be added to the log posterior
        */
          real asym_laplace_lcdf(real y, real mu, real sigma, real tau) {
            if (y < mu) {
              return log(tau) + (1 - tau) * (y - mu) / sigma;
            } else {
              return log1m((1 - tau) * exp(-tau * (y - mu) / sigma));
            }
          }
        /* asymmetric laplace log-CCDF for a single quantile
        * Args:
          *   y: a quantile
        *   mu: location parameter
        *   sigma: positive scale parameter
        *   tau: quantile parameter in (0, 1)
        * Returns:
          *   a scalar to be added to the log posterior
        */
          real asym_laplace_lccdf(real y, real mu, real sigma, real tau) {
            if (y < mu) {
              return log1m(tau * exp((1 - tau) * (y - mu) / sigma));
            } else {
              return log1m(tau) - tau * (y - mu) / sigma;
            }
          }
         
         real GAL2_lpdf(real y, real mu, real sigma, real gam, real tau){
         
         real p_pos;
         real p_neg;
         real a3;
         real a2;
         real p;
         real est;
         real A;
         real B;
         real Res = 0;
         //real gam = ligam;
          p = 1 * (gam < 0) + (tau - 1 * (gam < 0))/(2*Phi(-fabs(gam))*exp(.5*pow(gam, 2)));
          p_pos = p -  1 * (gam > 0);
          p_neg = p -  1 * (gam < 0);  
          est = (y - mu) / sigma; 
          
          if(fabs(gam) > 0){
          a3 = p_pos * (est / fabs(gam));
          a2 = fabs(gam) * (p_neg / p_pos);
          
          if(est/gam > 0){
            A =  0.5*pow(gam, 2)*pow(p_neg/p_pos, 2) - est*p_neg + log_diff_exp(log(Phi_approx(a2-a3)), log(Phi_approx(a2)) ); 
            B =  0.5*pow(gam, 2) - p_pos*est + log(Phi_approx(-fabs(gam) + a3));
            Res = log(2*p*(1-p)) - log(sigma) +  log_sum_exp(A, B);
          }else{
            Res =  log(2*p*(1-p)) - log(sigma) - p_pos * est + 0.5 * pow(gam, 2) + log(Phi_approx(-fabs(gam) ));
          }
          }else{
          Res = asym_laplace_lpdf( y | mu, sigma, tau); 
          }
           
          return Res;
         }
        
        real GAL2_rng(real mu, real sigma, real gam, real tau){
        
           real A;
           real B;
           real C;
           real p;
           real hi;
           real nui;
           real mui=0;
           real Up = uniform_rng(.5, 1.0);
           
          // real gam = ligam;
           p = (gam < 0) + (tau - (gam < 0))/(2*Phi_approx(-fabs(gam))*exp(.5*pow(gam, 2)));
           A = (1 - 2*p)/(p - pow(p,2));
           B = 2/(p - pow(p,2));
           C = 1/((gam > 0) - p);
           
            hi = sigma * inv_Phi(Up);
           nui = sigma * exponential_rng(1);
           mui += mu + A * nui + C * fabs(gam) * hi;
        
           return normal_rng(mui, sqrt(sigma*B*nui));
        }
        
      }
      data {
        int<lower=1> pn; //Number of basis;
        int<lower=1> N;  // total number of observations
        vector[N] Y;  // response variable
        int<lower=1> K;  // number of population-level effects
        matrix[N, K] X;  // population-level design matrix
        int prior_only;  // should the likelihood be ignored?
        matrix[N,pn] Wn;
          vector[pn] M;
         matrix[pn, pn] V;
         vector[2] Bd;
         real<lower=0,upper=1> tau0;
         int<lower=1> Kep;
      }
      transformed data {
        int Kc = K - 1;
        matrix[N, Kc] Xc;  // centered version of X without an intercept
        vector[Kc] means_X;  // column means of X before centering
        for (i in 2:K) {
          means_X[i - 1] = mean(X[, i]);
          Xc[, i - 1] = X[, i] - means_X[i - 1];
        }
      }
      parameters {
        simplex[Kep] pis;
        vector[Kc] b;  // population-level effects
        real Intercept;  // temporary intercept for centered predictors
        vector[pn] bw;
        real<lower=0.01,upper=.99> psi;
         
        // standarized spline coefficients
        vector<lower=0>[Kep] sigma;  // dispersion parameter
        //ordered<lower=Bd[1],upper=Bd[2]>[Kep] gam;
        vector<lower=Bd[1],upper=Bd[2]>[Kep] gam;
        real<lower=0> teta;
      }
      transformed parameters {
        //bwa ar(1) stationary prior
         vector[pn] bwa;
        real tau = tau0;
        real lprior = 0;  // prior contributions to the log posterior
        // compute bwa
        bwa[1] = bw[1];
        for(i in 2:pn)
          bwa[i] = psi*bwa[i-1] + bw[i]*teta; 
          
        lprior += normal_lpdf(b | 0.0, 10);
        lprior += std_normal_lpdf(bw);
        lprior += uniform_lpdf(psi | 0.01, .99);
        lprior += uniform_lpdf(gam | Bd[1], Bd[2]);
        //lprior += cauchy_lpdf(teta| 0, 1) - 1* cauchy_lccdf(0| 0, 1);
        //lprior += gamma_lpdf(teta| 1, 1);
        //lprior += weibull_lpdf(teta| 1, 2.5);
        lprior += std_normal_lpdf(teta) - 1* log(0.5);
        lprior += student_t_lpdf(Intercept | 3, 28.7, 5.6);
        //lprior += normal_lpdf(bs | 0, 5);
        // lprior += student_t_lpdf(sds_1_1 | 3, 0, 5.6)
        //   - 1 * student_t_lccdf(0 | 3, 0, 5.6);
        // lprior += student_t_lpdf(sds_2_1 | 3, 0, 5.6)
        //   - 1 * student_t_lccdf(0 | 3, 0, 5.6);
        lprior += student_t_lpdf(sigma | 3, 0, 5.6)
          - 1 * student_t_lccdf(0 | 3, 0, 5.6);
      }
      model {
        // likelihood including constants
        if (!prior_only) {
         // pis
         
          vector[Kep] log_pis = log(pis);
          // initialize linear predictor term + Xs * bs
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b + Wn * bwa;
          for (n in 1:N) {
          //  target += GAL2_lpdf(Y[n] | mu[n], sigma, ligam, tau0);
            vector[Kep] lps = log_pis;
            for(l in 1:Kep)
              lps[l] +=  GAL2_lpdf(Y[n] | mu[n], sigma[l], gam[l], tau0);
              
          target += log_sum_exp(lps);
          //log_lik[n] = log_sum_exp(lps);
          }
        }
        // priors including constants
        target += lprior;
        //target += std_normal_lpdf(zs_1_1);
        //target += std_normal_lpdf(zs_2_1);
      }
      generated quantities {
        // actual population-level intercept
        real b_Intercept = Intercept - dot_product(means_X, b);
        vector[N] mu = rep_vector(0.0, N);
        vector[N] log_lik;
        vector[Kep] lps0 = log(pis);
        vector[Kep] lps1 = log(pis);
          mu += Intercept + Xc * b + Wn * bwa;
         
        for(n in 1:N) {
             lps1 = log(pis);
            for(l in 1:Kep)
              lps1[l] +=  GAL2_lpdf(Y[n] | mu[n], sigma[l], gam[l], tau0);
              
          log_lik[n] = log_sum_exp(lps1);
          }
        //Vector[N] ASTP = Xs[1:N,1] * bs[1] + Zs_1_1 * s_1_1;
        //vector[N] MVPA = Xs[1:N,2] * bs[2] + Zs_2_1 * s_2_1;
      }
      "
  }
  df2b <- df_sampleWt
  
  form1 <- paste0("Y ~ ", paste(c("Gender","Race","HealthCondt2","AgeYR", colnames(df2b)[grepl("W",colnames(df2b))]), collapse = " + "))
  Dt2b <- make_standata(bf(as.formula(form1), quantile = tau0), data = df2b[c(1:nrow(df2b)),], family = asym_laplace())
  
  Dt2b$pn = pn
  Dt2b$Xold <- Dt2b$X
  Dt2b$X <- Dt2b$Xold[,!grepl("W.[0-9]", colnames(Dt2b$Xold))]
  Dt2b$Wn <- Dt2b$Xold[,grepl("W.[0-9]", colnames(Dt2b$Xold))]; dim(Dt2b$Wn)
  Dt2b$M <- numeric(Dt2b$pn)
  Dt2b$V <- as.positive.definite(as.inverse(P.mat(pn) + diag(rep(.01, pn))))*.1
  Dt2b$K <- ncol(Dt2b$X)
  Bd <- GamBnd(tau0)[1:2]
  Dt2b$Bd <- Bd*.9
  Dt2b$tau0 <- tau0
  Dt2b$X[,ncol(Dt2b$X)] <- Dt2b$X[,ncol(Dt2b$X)]/100
  Dt2b$Kep <- 2
  
  fitCase2b[[l]] <- stan(model_code = StandCodeVGalLin,
                         data = Dt2b, iter = Niter, chains = Nchain, verbose = TRUE,control = list(adapt_delta = 0.99,max_treedepth=11), cores = 2, seed = 1123)
  
  #-------------------------------------------------------------------------------
  # Case 3
  #-------------------------------------------------------------------------------
  formCase3 <- paste0("Y ~ ", paste(c("Gender","Race","HealthCondt2","AgeYR", colnames(df_sampleWt)[grepl("W",colnames(df_sampleWt))]), collapse = " + "))
  for(pa in varX){
    formCase3 <- paste0(formCase3, " + s(", pa, ", bs='ps', k=11)")
  }
  
  df_sampleWtSd <- df_sampleWt
  df_sampleWtSd$AgeYR <- df_sampleWtSd$AgeYR/sd(df_sampleWtSd$AgeYR)
  
  stanvars2 <- stanvar(scode = stan_funs2, block = "functions")
  fitCase3[[l]] <- brms::brm(bf(as.formula(formCase3), tau = tau0), data = df_sampleWt, family = GAL2, 
                             stanvars = stanvars2,chains = Nchain,iter = Niter, control = list(adapt_delta = 0.99),  cores = 2, seed = 1123) # init = 0.1,
  
  
  
  
  #-------------------------------------------------------------------------------
  # Case 4
  #-------------------------------------------------------------------------------
  df4 <- df_sample
  df4$AgeYRTru <- df4$AgeYR
  #df4$AgeYR <- scale(df4$AgeYR, center = T, scale = T)
  
  df4$AgeYR <- df_sample$AgeYR/100
  df4$MnASTP <- df_sample$MnASTP*6
  df4$MnMVPA <-  df_sample$MnMVPA/100
  
  df_sampleWt4 <- data.frame(df4, W = Wt)
  #tau0 = .1
  
  form4 <- paste0("Y ~ ", paste(c("Gender","Race","HealthCondt2","AgeYR", colnames(df_sampleWt4)[grepl("W",colnames(df_sampleWt4))]), collapse = " + "))
  for(pa in varX){
    form4 <- paste0(form4, " + s(", pa, ", bs='ps', k=11)")
  }
  
  
  Dt4 <- brms::make_standata(bf(as.formula(form4), quantile = tau0), data = df_sampleWt4, family = asym_laplace())
  #, stanvars = stanvars2, prior = bprior,chains = 2,iter = 500, control = list(adapt_delta = 0.97), data2 = list(M = rep(0, 29), V = SigbPrio))
  Bd = GamBnd(tau0)[1:2]
  
  Dt4$pn = pn
  Dt4$Xold <- Dt4$X
  Dt4$X <- Dt4$Xold[,!grepl("W.[0-9]", colnames(Dt4$Xold))]
  Dt4$Wn <- Dt4$Xold[,grepl("W.[0-9]", colnames(Dt4$Xold))]; dim(Dt4$Wn)
  Dt4$M <- numeric(Dt4$pn)
  Dt4$V <- as.positive.definite(as.inverse(P.mat(pn) + diag(rep(.01, pn))))*.1
  Dt4$K <- ncol(Dt4$X)
  Bd <- GamBnd(tau0)[1:2]
  Dt4$Bd <- Bd*.9
  Dt4$tau0 <- tau0
  
  #-- Option 2 GAL
  #fitXt[[l]]
  {
    StandCodeVGal <- "
  // > fit_fakeqGal3a.5$model
  // generated with brms 2.18.0
  functions {
    
    /*
    A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + log1m_exp(fabs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2)))); 
    gam = (gamU - gamL) * ligam + gamL;
    real gam = (gamU - gamL) * ligam + gamL;
    real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau, real gamL, real gamU){
    real GAL2_rng(real mu, real sigma, real ligam, real tau, real gamL, real gamU){
     */
       /* helper function for asym_laplace_lpdf
    * Args:
      *   y: the response value
    *   tau: quantile parameter in (0, 1)
    */
      real rho_quantile(real y, real tau) {
        if (y < 0) {
          return y * (tau - 1);
        } else {
          return y * tau;
        }
      }
    /* asymmetric laplace log-PDF for a single response
    * Args:
      *   y: the response value
    *   mu: location parameter
    *   sigma: positive scale parameter
    *   tau: quantile parameter in (0, 1)
    * Returns:
      *   a scalar to be added to the log posterior
    */
      real asym_laplace_lpdf(real y, real mu, real sigma, real tau) {
        return log(tau * (1 - tau)) -
          log(sigma) -
          rho_quantile((y - mu) / sigma, tau);
      }
    /* asymmetric laplace log-CDF for a single quantile
    * Args:
      *   y: a quantile
    *   mu: location parameter
    *   sigma: positive scale parameter
    *   tau: quantile parameter in (0, 1)
    * Returns:
      *   a scalar to be added to the log posterior
    */
      real asym_laplace_lcdf(real y, real mu, real sigma, real tau) {
        if (y < mu) {
          return log(tau) + (1 - tau) * (y - mu) / sigma;
        } else {
          return log1m((1 - tau) * exp(-tau * (y - mu) / sigma));
        }
      }
    /* asymmetric laplace log-CCDF for a single quantile
    * Args:
      *   y: a quantile
    *   mu: location parameter
    *   sigma: positive scale parameter
    *   tau: quantile parameter in (0, 1)
    * Returns:
      *   a scalar to be added to the log posterior
    */
      real asym_laplace_lccdf(real y, real mu, real sigma, real tau) {
        if (y < mu) {
          return log1m(tau * exp((1 - tau) * (y - mu) / sigma));
        } else {
          return log1m(tau) - tau * (y - mu) / sigma;
        }
      }
     
     real GAL2_lpdf(real y, real mu, real sigma, real gam, real tau){
     
     real p_pos;
     real p_neg;
     real a3;
     real a2;
     real p;
     real est;
     real A;
     real B;
     real Res = 0;
     //real gam = ligam;
      p = 1 * (gam < 0) + (tau - 1 * (gam < 0))/(2*Phi(-fabs(gam))*exp(.5*pow(gam, 2)));
      p_pos = p -  1 * (gam > 0);
      p_neg = p -  1 * (gam < 0);  
      est = (y - mu) / sigma; 
      
      if(fabs(gam) > 0){
      a3 = p_pos * (est / fabs(gam));
      a2 = fabs(gam) * (p_neg / p_pos);
      
      
      if(est/gam > 0){
        A =  0.5*pow(gam, 2)*pow(p_neg/p_pos, 2) - est*p_neg + log_diff_exp(log(Phi_approx(a2-a3)), log(Phi_approx(a2)) ); 
        B =  0.5*pow(gam, 2) - p_pos*est + log(Phi_approx(-fabs(gam) + a3));
        Res = log(2*p*(1-p)) - log(sigma) +  log_sum_exp(A, B);
      }else{
        Res =  log(2*p*(1-p)) - log(sigma) - p_pos * est + 0.5 * pow(gam, 2) + log(Phi_approx(-fabs(gam) ));
      }
      }else{
      Res = asym_laplace_lpdf( y | mu, sigma, tau); 
      }
       
      return Res;
     }
    
    real GAL2_rng(real mu, real sigma, real gam, real tau){
    
       real A;
       real B;
       real C;
       real p;
       real hi;
       real nui;
       real mui=0;
       real Up = uniform_rng(.5, 1.0);
       
      // real gam = ligam;
       p = (gam < 0) + (tau - (gam < 0))/(2*Phi_approx(-fabs(gam))*exp(.5*pow(gam, 2)));
       A = (1 - 2*p)/(p - pow(p,2));
       B = 2/(p - pow(p,2));
       C = 1/((gam > 0) - p);
       
        hi = sigma * inv_Phi(Up);
       nui = sigma * exponential_rng(1);
       mui += mu + A * nui + C * fabs(gam) * hi;
    
       return normal_rng(mui, sqrt(sigma*B*nui));
    }
    
  }
  data {
  int<lower=1> pn; //Number of basis;
    int<lower=1> N;  // total number of observations
    vector[N] Y;  // response variable
    int<lower=1> K;  // number of population-level effects
    matrix[N, K] X;  // population-level design matrix
    // data for splines
    int Ks;  // number of linear effects
    matrix[N, Ks] Xs;  // design matrix for the linear effects
     
    int nb_1;  // number of bases
    int knots_1[nb_1];  // number of knots
    // basis function matrices
    matrix[N, knots_1[1]] Zs_1_1;
    
    int nb_2;  // number of bases
    int knots_2[nb_2];  // number of knots
    // basis function matrices
    matrix[N, knots_2[1]] Zs_2_1;
    int prior_only;  // should the likelihood be ignored?
    matrix[N,pn] Wn;
      vector[pn] M;
     matrix[pn, pn] V;
     vector[2] Bd;
     real<lower=0,upper=1> tau0;
  }
  transformed data {
    int Kc = K - 1;
    matrix[N, Kc] Xc;  // centered version of X without an intercept
    vector[Kc] means_X;  // column means of X before centering
    for (i in 2:K) {
      means_X[i - 1] = mean(X[, i]);
      Xc[, i - 1] = X[, i] - means_X[i - 1];
    }
  }
  parameters {
    vector[Kc] b;  // population-level effects
    real Intercept;  // temporary intercept for centered predictors
    vector[Ks] bs;  // spline coefficients
    vector[pn] bw;
    real<lower=0.01,upper=.99> psi;
     
    // standarized spline coefficients
    vector[knots_1[1]] zs_1_1;
    real<lower=0> sds_1_1;  // standard deviations of spline coefficients
     
    // standarized spline coefficients
    vector[knots_2[1]] zs_2_1;
    real<lower=0> sds_2_1;  // standard deviations of spline coefficients
    real<lower=0> sigma;  // dispersion parameter
    real<lower=Bd[1],upper=Bd[2]> gam;
    real<lower=0> teta;
  }
  transformed parameters {
     //bwa ar(1) stationary prior
     vector[pn] bwa;
    // actual spline coefficients
    vector[knots_1[1]] s_1_1;
    // actual spline coefficients
    vector[knots_2[1]] s_2_1;
    //matrix[pn,pn] V1;
    real tau = tau0;
    real lprior = 0;  // prior contributions to the log posterior
    // compute actual spline coefficients
    s_1_1 = sds_1_1 * zs_1_1;
    //V1 = teta*V;
    // compute actual spline coefficients
    s_2_1 = sds_2_1 * zs_2_1;
    // compute bwa
    bwa[1] = bw[1];
    for(i in 2:pn)
      bwa[i] = psi*bwa[i-1] + bw[i]*teta; 
      
    lprior += normal_lpdf(b | 0.0, 10);
    lprior += std_normal_lpdf(bw);
    lprior += uniform_lpdf(psi | 0.01, .99);
    lprior += uniform_lpdf(gam | Bd[1], Bd[2]);
    //lprior += cauchy_lpdf(teta| 0, 1) - 1* cauchy_lccdf(0| 0, 1);
    //lprior += gamma_lpdf(teta| 1, 1);
    //lprior += weibull_lpdf(teta| 1, 2.5);
    lprior += std_normal_lpdf(teta) - 1* log(0.5);
    lprior += student_t_lpdf(Intercept | 3, 28.7, 5.6);
    lprior += normal_lpdf(bs | 0, 5);
    lprior += student_t_lpdf(sds_1_1 | 3, 0, 5.6)
      - 1 * student_t_lccdf(0 | 3, 0, 5.6);
    lprior += student_t_lpdf(sds_2_1 | 3, 0, 5.6)
      - 1 * student_t_lccdf(0 | 3, 0, 5.6);
    lprior += student_t_lpdf(sigma | 3, 0, 5.6)
      - 1 * student_t_lccdf(0 | 3, 0, 5.6);
  }
  model {
    // likelihood including constants
    if (!prior_only) {
      // initialize linear predictor term
      vector[N] mu = rep_vector(0.0, N);
      mu += Intercept + Xc * b + Wn * bwa + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;
      for (n in 1:N) {
      //  target += GAL2_lpdf(Y[n] | mu[n], sigma, ligam, tau0);
          target += GAL2_lpdf(Y[n] | mu[n], sigma, gam, tau0);
      }
    }
    // priors including constants
    target += lprior;
    target += std_normal_lpdf(zs_1_1);
    target += std_normal_lpdf(zs_2_1);
  }
  generated quantities {
    // actual population-level intercept
    real b_Intercept = Intercept - dot_product(means_X, b);
    vector[N] mu = rep_vector(0.0, N);
    vector[N] log_lik;
          mu += Intercept + Xc * b + Wn * bwa;
        for(n in 1:N)
              log_lik[n] +=  GAL2_lpdf(Y[n] | mu[n], sigma, gam, tau0);

    //Vector[N] ASTP = Xs[1:N,1] * bs[1] + Zs_1_1 * s_1_1;
    //vector[N] MVPA = Xs[1:N,2] * bs[2] + Zs_2_1 * s_2_1;
  }
  "
  }
  
  fitCase4[[l]] <- stan(model_code = StandCodeVGal, 
                        data = Dt4, iter = Niter, chains = Nchain, verbose = TRUE,control = list(adapt_delta = 0.99,max_treedepth=11), cores = 2, seed = 1123) 
  
  
  
  #-------------------------------------------------------------------------------
  # Case 5
  #-------------------------------------------------------------------------------
  df5 <- df_sample
  df5$AgeYRTru <- df5$AgeYR
  #df5$AgeYR <- scale(df5$AgeYR, center = T, scale = T)
  
  df5$AgeYR <- df_sample$AgeYR/100
  df5$MnASTP <- df_sample$MnASTP*6
  df5$MnMVPA <-  df_sample$MnMVPA/100
  
  df_sampleWt5 <- data.frame(df5, W = Wt)
  #tau0 = .1
  
  form5 <- paste0("Y ~ ", paste(c("Gender","Race","HealthCondt2","AgeYR", colnames(df_sampleWt5)[grepl("W",colnames(df_sampleWt4))]), collapse = " + "))
  for(pa in varX){
    form5 <- paste0(form5, " + s(", pa, ", bs='ps', k=11)")
  }
  
  
  Dt5 <- brms::make_standata(bf(as.formula(form4), quantile = tau0), data = df_sampleWt5, family = asym_laplace())
  #, stanvars = stanvars2, prior = bprior,chains = 2,iter = 500, control = list(adapt_delta = 0.97), data2 = list(M = rep(0, 29), V = SigbPrio))
  Bd = GamBnd(tau0)[1:2]
  
  Dt5$pn = pn
  Dt5$Xold <- Dt5$X
  Dt5$X <- Dt5$Xold[,!grepl("W.[0-9]", colnames(Dt5$Xold))]
  Dt5$Wn <- Dt5$Xold[,grepl("W.[0-9]", colnames(Dt5$Xold))]; dim(Dt5$Wn)
  Dt5$M <- numeric(Dt5$pn)
  Dt5$V <- as.positive.definite(as.inverse(P.mat(pn) + diag(rep(.01, pn))))*.1
  Dt5$K <- ncol(Dt5$X)
  Bd <- GamBnd(tau0)[1:2]
  Dt5$Bd <- Bd*.9
  Dt5$tau0 <- tau0
  
  #-- Option 2 GAL
  #fitXt[[l]]
  {
    StandCodeVGal5 <- "
    // > fit_fakeqGal3a.5$model
    // generated with brms 2.18.0
    functions {
      
      /*
      A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + log1m_exp(fabs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2)))); 
      gam = (gamU - gamL) * ligam + gamL;
      real gam = (gamU - gamL) * ligam + gamL;
      real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau, real gamL, real gamU){
      real GAL2_rng(real mu, real sigma, real ligam, real tau, real gamL, real gamU){
       */
         /* helper function for asym_laplace_lpdf
      * Args:
        *   y: the response value
      *   tau: quantile parameter in (0, 1)
      */
        real rho_quantile(real y, real tau) {
          if (y < 0) {
            return y * (tau - 1);
          } else {
            return y * tau;
          }
        }
      /* asymmetric laplace log-PDF for a single response
      * Args:
        *   y: the response value
      *   mu: location parameter
      *   sigma: positive scale parameter
      *   tau: quantile parameter in (0, 1)
      * Returns:
        *   a scalar to be added to the log posterior
      */
        real asym_laplace_lpdf(real y, real mu, real sigma, real tau) {
          return log(tau * (1 - tau)) -
            log(sigma) -
            rho_quantile((y - mu) / sigma, tau);
        }
      /* asymmetric laplace log-CDF for a single quantile
      * Args:
        *   y: a quantile
      *   mu: location parameter
      *   sigma: positive scale parameter
      *   tau: quantile parameter in (0, 1)
      * Returns:
        *   a scalar to be added to the log posterior
      */
        real asym_laplace_lcdf(real y, real mu, real sigma, real tau) {
          if (y < mu) {
            return log(tau) + (1 - tau) * (y - mu) / sigma;
          } else {
            return log1m((1 - tau) * exp(-tau * (y - mu) / sigma));
          }
        }
      /* asymmetric laplace log-CCDF for a single quantile
      * Args:
        *   y: a quantile
      *   mu: location parameter
      *   sigma: positive scale parameter
      *   tau: quantile parameter in (0, 1)
      * Returns:
        *   a scalar to be added to the log posterior
      */
        real asym_laplace_lccdf(real y, real mu, real sigma, real tau) {
          if (y < mu) {
            return log1m(tau * exp((1 - tau) * (y - mu) / sigma));
          } else {
            return log1m(tau) - tau * (y - mu) / sigma;
          }
        }
       
       real GAL2_lpdf(real y, real mu, real sigma, real gam, real tau){
       
       real p_pos;
       real p_neg;
       real a3;
       real a2;
       real p;
       real est;
       real A;
       real B;
       real Res = 0;
       //real gam = ligam;
        p = 1 * (gam < 0) + (tau - 1 * (gam < 0))/(2*Phi(-fabs(gam))*exp(.5*pow(gam, 2)));
        p_pos = p -  1 * (gam > 0);
        p_neg = p -  1 * (gam < 0);  
        est = (y - mu) / sigma; 
        
        if(fabs(gam) > 0){
        a3 = p_pos * (est / fabs(gam));
        a2 = fabs(gam) * (p_neg / p_pos);
        
        
        if(est/gam > 0){
          A =  0.5*pow(gam, 2)*pow(p_neg/p_pos, 2) - est*p_neg + log_diff_exp(log(Phi_approx(a2-a3)), log(Phi_approx(a2)) ); 
          B =  0.5*pow(gam, 2) - p_pos*est + log(Phi_approx(-fabs(gam) + a3));
          Res = log(2*p*(1-p)) - log(sigma) +  log_sum_exp(A, B);
        }else{
          Res =  log(2*p*(1-p)) - log(sigma) - p_pos * est + 0.5 * pow(gam, 2) + log(Phi_approx(-fabs(gam) ));
        }
        }else{
        Res = asym_laplace_lpdf( y | mu, sigma, tau); 
        }
         
        return Res;
       }
      
      real GAL2_rng(real mu, real sigma, real gam, real tau){
      
         real A;
         real B;
         real C;
         real p;
         real hi;
         real nui;
         real mui=0;
         real Up = uniform_rng(.5, 1.0);
         
        // real gam = ligam;
         p = (gam < 0) + (tau - (gam < 0))/(2*Phi_approx(-fabs(gam))*exp(.5*pow(gam, 2)));
         A = (1 - 2*p)/(p - pow(p,2));
         B = 2/(p - pow(p,2));
         C = 1/((gam > 0) - p);
         
          hi = sigma * inv_Phi(Up);
         nui = sigma * exponential_rng(1);
         mui += mu + A * nui + C * fabs(gam) * hi;
      
         return normal_rng(mui, sqrt(sigma*B*nui));
      }
      
    }
    data {
    int<lower=1> pn; //Number of basis;
      int<lower=1> N;  // total number of observations
      vector[N] Y;  // response variable
      int<lower=1> K;  // number of population-level effects
      matrix[N, K] X;  // population-level design matrix
      // data for splines
      int Ks;  // number of linear effects
      matrix[N, Ks] Xs;  // design matrix for the linear effects
       
      int nb_1;  // number of bases
      int knots_1[nb_1];  // number of knots
      // basis function matrices
      matrix[N, knots_1[1]] Zs_1_1;
      
      int nb_2;  // number of bases
      int knots_2[nb_2];  // number of knots
      // basis function matrices
      matrix[N, knots_2[1]] Zs_2_1;
      int prior_only;  // should the likelihood be ignored?
      matrix[N,pn] Wn;
        vector[pn] M;
       matrix[pn, pn] V;
       vector[2] Bd;
       real<lower=0,upper=1> tau0;
    }
    transformed data {
      int Kc = K - 1;
      matrix[N, Kc] Xc;  // centered version of X without an intercept
      vector[Kc] means_X;  // column means of X before centering
      for (i in 2:K) {
        means_X[i - 1] = mean(X[, i]);
        Xc[, i - 1] = X[, i] - means_X[i - 1];
      }
    }
    parameters {
      vector[Kc] b;  // population-level effects
      real Intercept;  // temporary intercept for centered predictors
      vector[Ks] bs;  // spline coefficients
      vector<upper=0>[pn] bw;
      real<lower=0.01,upper=.99> psi;
       
      // standarized spline coefficients
      vector[knots_1[1]] zs_1_1;
      real<lower=0> sds_1_1;  // standard deviations of spline coefficients
       
      // standarized spline coefficients
      vector[knots_2[1]] zs_2_1;
      real<lower=0> sds_2_1;  // standard deviations of spline coefficients
      real<lower=0> sigma;  // dispersion parameter
      real<lower=Bd[1],upper=Bd[2]> gam;
      real<lower=0> teta;
    }
    transformed parameters {
       //bwa ar(1) stationary prior
       vector<upper=0>[pn] bwa;
      // actual spline coefficients
      vector[knots_1[1]] s_1_1;
      // actual spline coefficients
      vector[knots_2[1]] s_2_1;
      //matrix[pn,pn] V1;
      real tau = tau0;
      real lprior = 0;  // prior contributions to the log posterior
      // compute actual spline coefficients
      s_1_1 = sds_1_1 * zs_1_1;
      //V1 = teta*V;
      // compute actual spline coefficients
      s_2_1 = sds_2_1 * zs_2_1;
      // compute bwa
      bwa[1] = bw[1];
      for(i in 2:pn)
        bwa[i] = psi*bwa[i-1] + bw[i]*teta; 
        
      lprior += normal_lpdf(b | 0.0, 10);
      lprior += std_normal_lpdf(bw) - 1* log(0.5); // imposing bw < 0
      lprior += uniform_lpdf(psi | 0.01, .99);
      lprior += uniform_lpdf(gam | Bd[1], Bd[2]);
      //lprior += cauchy_lpdf(teta| 0, 1) - 1* cauchy_lccdf(0| 0, 1);
      //lprior += gamma_lpdf(teta| 1, 1);
      //lprior += weibull_lpdf(teta| 1, 2.5);
      lprior += std_normal_lpdf(teta) - 1* log(0.5);
      lprior += student_t_lpdf(Intercept | 3, 28.7, 5.6);
      lprior += normal_lpdf(bs | 0, 5);
      lprior += student_t_lpdf(sds_1_1 | 3, 0, 5.6)
        - 1 * student_t_lccdf(0 | 3, 0, 5.6);
      lprior += student_t_lpdf(sds_2_1 | 3, 0, 5.6)
        - 1 * student_t_lccdf(0 | 3, 0, 5.6);
      lprior += student_t_lpdf(sigma | 3, 0, 5.6)
        - 1 * student_t_lccdf(0 | 3, 0, 5.6);
    }
    model {
      // likelihood including constants
      if (!prior_only) {
        // initialize linear predictor term
        vector[N] mu = rep_vector(0.0, N);
        mu += Intercept + Xc * b + Wn * bwa + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;
        for (n in 1:N) {
        //  target += GAL2_lpdf(Y[n] | mu[n], sigma, ligam, tau0);
            target += GAL2_lpdf(Y[n] | mu[n], sigma, gam, tau0);
        }
      }
      // priors including constants
      target += lprior;
      target += std_normal_lpdf(zs_1_1);
      target += std_normal_lpdf(zs_2_1);
    }
    generated quantities {
      // actual population-level intercept
      real b_Intercept = Intercept - dot_product(means_X, b);
      vector[N] mu = rep_vector(0.0, N);
      vector[N] log_lik;
          mu += Intercept + Xc * b + Wn * bwa;
        for(n in 1:N)
              log_lik[n] +=  GAL2_lpdf(Y[n] | mu[n], sigma, gam, tau0);

      //Vector[N] ASTP = Xs[1:N,1] * bs[1] + Zs_1_1 * s_1_1;
      //vector[N] MVPA = Xs[1:N,2] * bs[2] + Zs_2_1 * s_2_1;
    }
    "
  }
  
  fitCase5[[l]] <- stan(model_code = StandCodeVGal5, 
                        data = Dt5, iter = Niter, chains = Nchain, verbose = TRUE,control = list(adapt_delta = 0.99,max_treedepth=11), cores = 2, seed = 1123) 
  
  
  
  
  
}
