
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


Xt <- apply(Warray, c(1,2), mean)
#Xt =  as.matrix(FUI(model="gaussian", smooth=T, Warray ,silent = FALSE)) #} #- (abs(range(W)[1]) +1)
#Xt <- as.matrix(FUI(model="gaussian",smooth=T, Warray,silent = FALSE))
Wt <- (Xt%*%bs2)/length(a)


#---- Transform the data
# df_sample$AgeYRTru <- df_sample$AgeYR
# #df_sample$AgeYR <- scale(df_sample$AgeYR, center = T, scale = T)
# df_sampleOld <- df_sample
# df_sample$AgeYR <- df_sample$AgeYR/100
# df_sample$MnASTP <- df_sample$MnASTP*10
# df_sample$MnMVPA <-  df_sample$MnMVPA/100

Nsub = 300 

df_sampleWt <- data.frame(df_sample[c(1:Nsub),], W = Wt[c(1:Nsub),])



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
fitCase3  <- list()
fitCase4  <- list()
fitCase5  <- list()
fitCase2b <- list()
tauVal <- c(.1,.5,.9)

Niter = 5000
Nchain = 2

#for(l in 1:length(tauVal)){ 
l = 3  
tau0 <- tauVal[l]
Bd = GamBnd(tau0)[1:2]




#-------------------------------------------------------------------------------
# Case 2b  (with mizture distribution)
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

library(loo)
library(StanDDM)

?  waic()

LogLik2b <- extract_log_lik(fitCase2b[[l]], parameter_name = "log_lik")  
loo::waic(LogLik2b) 
LooVal <- loo::loo(LogLik2b, save_psis = TRUE)
print(LooVal)
