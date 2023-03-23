#-- Implementing Model in Rstan
library(brms)
# library(rstan)
# library(splines)

#-------------------------------------------------------------------------------
#--- Import Data
#-------------------------------------------------------------------------------
# source("~/Library/CloudStorage/OneDrive-IndianaUniversity/Join_folder_Annie_DrZoh/NHANES Dataset/2013-2014/BQReg_Nhanes2012_13Application.R")
# load("/Users/rszoh/Library/CloudStorage/OneDrive-IndianaUniversity/Join_folder_Annie_DrZoh/NHANES Dataset/2013-2014/FinalDataV2.RData")

#-------------------------------------------------------------------------------
#--- Set up the data
#-------------------------------------------------------------------------------
# a <- seq(0,1, length = dim(Warray)[2])
# bs2 <- bs(a, df = 20, intercept = T)


# Xt <- apply(Warray, c(1,2), mean)
# #Xt <- as.matrix(FUI(model="gaussian",smooth=T, Warray,silent = FALSE))
# Wt <- (Xt%*%bs2)/length(a)

# df_sampleWt <- data.frame(df_sample, W = Wt)
# form <- paste0("Y ~ ", paste(c("Gender","Race","HealthCondt2","AgeYR","SmokStat", colnames(df_sampleWt)[grepl("W",colnames(df_sampleWt))]), collapse = " + "))
# for(pa in varX){
#   form <- paste0(form, " + s(", pa, ", bs='ps', k=11)")
# }


# parameters

# iterations of the simulation procedure
N <- 100
# mean number of data points per value of x sampled
mean_obs <- 10
# true value of beta, as in mu_i = beta * x_i
beta_mean <- 1
# multiplier for heteroskedasticity
hetero_fact <- 1


# simulation procedure
x <- c()
y <- c()
for (i in 1:N) {
  # generate x value
  x_i <- runif(n = 1)
  # generate how many samples will have that x value
  n_i <- rpois(n = 1, lambda = mean_obs)
  if (n_i > 0) {
    # generate the y values of each data point
    y_i <- rnorm(
      n = n_i, mean = beta_mean * x_i,
      sd = hetero_fact * x_i
      # sd = .1
    )
    # append them to data
    for (j in 1:n_i) {
      x <- append(x, x_i)
      y <- append(y, y_i[j])
    }
  }
}

simulated_data <- data.frame(x = x, y = y)



#-------------------------------------------------------------------------------
#--- Useful functions
#-- find the root
#-------------------------------------------------------------------------------
GamF <- function(gam, p0) {
  (2 * pnorm(-abs(gam)) * exp(.5 * gam^2) - p0)^2
}

# optimize(GamF,interval = c(-30,30), p0=1 - 0.05)
GamBnd <- function(p0) {
  Re1 <- optimize(GamF, interval = c(-30, 30), p0 = 1 - p0)
  Re2 <- optimize(GamF, interval = c(-30, 30), p0 = p0)

  c(-abs(Re1$minimum), abs(Re2$minimum), Re1$objective, Re2$objective)
}
GamBnd <- Vectorize(GamBnd, "p0")

#-------------------------------------------------------------------------------
#--- Function
#-------------------------------------------------------------------------------
#--- Stan function (case 2)
tau0 <- .1

{
  Bd <- GamBnd(tau0)[1:2]
  # GAL2 <- custom_family(
  #  "GAL2", dpars = c("mu", "sigma","ligam", "tau","gamL","gamU"), links=c("identity","log","identity","identity", rep("identity",2)),
  #  lb=c(NA,0, -1.087,0,NA,NA), ub=c(NA,NA, 1.087,1,NA,NA), type="real") #, vars = "vint1[n]"
  # GAL2 <- custom_family(
  #   "GAL2", dpars = c("mu", "sigma","ligam"), links= c("identity","log","identity"),
  #   lb=c(NA,0, -1.087), ub=c(NA,NA, 1.087), type="real", vars = c("tau","gamL","gamU")) #, vars = "vint1[n]"
  GAL2 <- custom_family(
    "GAL2",
    dpars = c("mu", "sigma", "ligam", "tau"), links = c("identity", "log", "identity", "identity"),
    lb = c(NA, 0, -1.086, 0), ub = c(NA, NA, 1.086, 1), type = "real"
  ) # , vars = c("tau")) #, vars = "vint1[n]"

  GAL2 <- custom_family(
    "GAL2",
    dpars = c("mu", "sigma", "ligam", "tau"), links = c("identity", "log", "logit", "identity"),
    lb = c(NA, 0, -15, 0), ub = c(NA, NA, 15, 1), type = "real", vars = c("gamL", "gamU")
  ) # , vars = c("tau")) #, vars = "vint1[n]"

  GAL2 <- custom_family(
    "GAL2",
    dpars = c("mu", "sigma", "ligam", "tau", "gamL", "gamU"), links = c("identity", "log", "identity", "identity", rep("identity", 2)),
    lb = c(NA, 0, -15, 0, NA, NA), ub = c(NA, NA, 15, 1, NA, NA), type = "real"
  ) # , vars = "vint1[n]"

  GAL2 <- custom_family(
    "GAL2",
    dpars = c("mu", "sigma", "ligam", "tau"), links = c("identity", "log", "identity", "identity"),
    lb = c(NA, 0, Bd[1] * .99, 0), ub = c(NA, NA, Bd[2] * .99, 1), type = "real"
  ) # , vars = "vint1[n]"

  GAL2 <- custom_family(
    "GAL2",
    dpars = c("mu", "sigma", "ligam", "tau"), links = c("identity", "log", "identity", "identity"),
    lb = c(NA, 0, Bd[1] * .99, 0), ub = c(NA, NA, Bd[2] * .99, 1), type = "real"
  )

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

   real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau){

   real p_pos;
   real p_neg;
   real a3;
   real a2;
   real p;
   real est;
   real A;
   real B;
   real Res = 0;
   real gam = ligam;
    p = (gam < 0) + (tau - (gam < 0))/(2*Phi(-fabs(gam))*exp(.5*pow(gam, 2)));
    p_pos = p -  (gam > 0);
    p_neg = p -  (gam < 0);
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
# stanvars <- stanvar(scode = stan_funs, block = "functions")
stanvars2 <- stanvar(scode = stan_funs2, block = "functions")

# priorval = c(prior(normal(0, 15), class = "Intercept"),
#           # Prior guess of 80% of the terms are non-zero
#           prior(horseshoe(par_ratio =0.8), class = "b"))

#--- Fit the quantile regression
#-- quantile

# fit_fakeq.5 <- brms::brm(bf(as.formula(form), quantile = .5), data = df_sampleWt, family = asym_laplace(), prior = set_prior(horseshoe(df=1,par_ratio = 0.1), class="b"))
# fit_fakeq.5 <- brms::brm(bf(as.formula(form), quantile = .5), data = df_sampleWt, family = asym_laplace())

# fit_fakeqGal.5 <- brms::brm(bf(as.formula(form), tau = .5, gamL = Bd[1], gamU = Bd[2]), data = df_sampleWt, family = GAL,
#                             stanvars = stanvars,chains = 3,iter = 10000, control = list(adapt_delta = 0.97, max_treedepth = 10))


# fit_fakeqGal2.5 <- brms::brm(bf(as.formula(form), tau = .5, gamL = Bd[1], gamU = Bd[2]), data = df_sampleWt, family = GAL2, stanvars = stanvars2,chains = 3,iter = 1000)
# fit_fakeqGal2.5 <- brms::brm(bf(as.formula(form), tau = .5,gamL = Bd[1], gamU = Bd[2]), data = df_sampleWt, family = GAL2, stanvars = stanvars2,chains = 2,iter = 2000, init = 0.1, control = list(adapt_delta = 0.97))

# fit_fakeqGal2.5 <- brms::brm(bf(as.formula(form), tau = tau0), data = df_sampleWt, family = GAL2, stanvars = stanvars2,chains = 2,iter = 5000, control = list(adapt_delta = 0.97)) # init = 0.1,
# fit_fakeqGal2.9 <- brms::brm(bf(as.formula(form), tau = tau0), data = df_sampleWt, family = GAL2, stanvars = stanvars2,chains = 2,iter = 5000, control = list(adapt_delta = 0.97)) # init = 0.1,
# fit_fakeqGal2.1 <- brms::brm(bf(as.formula(form), tau = tau0), data = df_sampleWt, family = GAL2, stanvars = stanvars2,chains = 2,iter = 5000, control = list(adapt_delta = 0.97)) # init = 0.1,


fit_fakeqGal2.1 <- brms::brm(bf(y ~ x, sigma ~ log(x), tau = tau0), data = simulated_data, family = GAL2, stanvars = stanvars2, chains = 2, iter = 5000, control = list(adapt_delta = 0.97)) # init = 0.1,

summary(fit_fakeqGal2.1)
# 
# plot(x, y)
# lines(x, (1 + qnorm(.9)) * x, col = "blue")
# lines(x, 0.01 + 2.13 * x)
# 
stancode(fit_fakeqGal2.1)

# #brms::conditional_smooths(fit_fakeqGal2.5)
# #brms::conditional_smooths(fit_fakeqGal2.9)
# brms::conditional_smooths(fit_fakeqGal2.1)



# #********************************************************************************
#   # trial code
#   #
#  tau0 = 0.5

#   {
#     Bd = GamBnd(tau0)[1:2]
#     #GAL2 <- custom_family(
#     #  "GAL2", dpars = c("mu", "sigma","ligam", "tau","gamL","gamU"), links=c("identity","log","identity","identity", rep("identity",2)),
#     #  lb=c(NA,0, -1.087,0,NA,NA), ub=c(NA,NA, 1.087,1,NA,NA), type="real") #, vars = "vint1[n]"
#     # GAL2 <- custom_family(
#     #   "GAL2", dpars = c("mu", "sigma","ligam"), links= c("identity","log","identity"),
#     #   lb=c(NA,0, -1.087), ub=c(NA,NA, 1.087), type="real", vars = c("tau","gamL","gamU")) #, vars = "vint1[n]"
#     GAL2 <- custom_family(
#       "GAL2", dpars = c("mu", "sigma","ligam", "tau"), links= c("identity","log","identity","identity"),
#       lb=c(NA,0, -1.086,0), ub=c(NA,NA, 1.086,1), type="real",loop= FALSE) #, vars = c("tau")) #, vars = "vint1[n]"

#     GAL2 <- custom_family(
#       "GAL2", dpars = c("mu", "sigma","ligam", "tau"), links= c("identity","log","logit","identity"),
#       lb=c(NA,0, -15,0), ub=c(NA,NA, 15,1), type="real", vars=c("gamL", "gamU")) #, vars = c("tau")) #, vars = "vint1[n]"

#     GAL2 <- custom_family(
#       "GAL2", dpars = c("mu", "sigma","ligam", "tau","gamL","gamU"), links=c("identity","log","identity","identity", rep("identity",2)),
#       lb=c(NA,0, -15,0,NA,NA), ub=c(NA,NA, 15,1,NA,NA), type="real") #, vars = "vint1[n]"

#     GAL2 <- custom_family(
#       "GAL2", dpars = c("mu", "sigma","ligam", "tau"), links=c("identity","log","identity","identity"),
#       lb=c(NA,0, Bd[1]*.99,0), ub=c(NA,NA, Bd[2]*.99,1), type="real",loop= TRUE) #, vars = "vint1[n]"


#     stan_funs2 <- "
#   /*
#   A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + log1m_exp(fabs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2))));
#   gam = (gamU - gamL) * ligam + gamL;
#   real gam = (gamU - gamL) * ligam + gamL;
#   real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau, real gamL, real gamU){
#   real GAL2_rng(real mu, real sigma, real ligam, real tau, real gamL, real gamU){
#    */
#      /* helper function for asym_laplace_lpdf
#   * Args:
#     *   y: the response value
#   *   tau: quantile parameter in (0, 1)
#   */
#     real rho_quantile(real y, real tau) {
#       if (y < 0) {
#         return y * (tau - 1);
#       } else {
#         return y * tau;
#       }
#     }
#   /* asymmetric laplace log-PDF for a single response
#   * Args:
#     *   y: the response value
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   tau: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lpdf(real y, real mu, real sigma, real tau) {
#       return log(tau * (1 - tau)) -
#         log(sigma) -
#         rho_quantile((y - mu) / sigma, tau);
#     }
#   /* asymmetric laplace log-CDF for a single quantile
#   * Args:
#     *   y: a quantile
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   tau: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lcdf(real y, real mu, real sigma, real tau) {
#       if (y < mu) {
#         return log(tau) + (1 - tau) * (y - mu) / sigma;
#       } else {
#         return log1m((1 - tau) * exp(-tau * (y - mu) / sigma));
#       }
#     }
#   /* asymmetric laplace log-CCDF for a single quantile
#   * Args:
#     *   y: a quantile
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   tau: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lccdf(real y, real mu, real sigma, real tau) {
#       if (y < mu) {
#         return log1m(tau * exp((1 - tau) * (y - mu) / sigma));
#       } else {
#         return log1m(tau) - tau * (y - mu) / sigma;
#       }
#     }

#    real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau){

#    real p_pos;
#    real p_neg;
#    real a3;
#    real a2;
#    real p;
#    real est;
#    real A;
#    real B;
#    real Res = 0;
#    real gam = ligam;
#     p = (gam < 0) + (tau - (gam < 0))/(2*Phi(-fabs(gam))*exp(.5*pow(gam, 2)));
#     p_pos = p -  (gam > 0);
#     p_neg = p -  (gam < 0);
#     est = (y - mu) / sigma;

#     if(fabs(gam) > 0){
#     a3 = p_pos * (est / fabs(gam));
#     a2 = fabs(gam) * (p_neg / p_pos);


#     if(est/gam > 0){
#       A =  0.5*pow(gam, 2)*pow(p_neg/p_pos, 2) - est*p_neg + log_diff_exp(log(Phi_approx(a2-a3)), log(Phi_approx(a2)) );
#       B =  0.5*pow(gam, 2) - p_pos*est + log(Phi_approx(-fabs(gam) + a3));
#       Res = log(2*p*(1-p)) - log(sigma) +  log_sum_exp(A, B);
#     }else{
#       Res =  log(2*p*(1-p)) - log(sigma) - p_pos * est + 0.5 * pow(gam, 2) + log(Phi_approx(-fabs(gam) ));
#     }
#     }else{
#     Res = asym_laplace_lpdf( y | mu, sigma, tau);
#     }

#     return Res;
#    }

#   real GAL2_rng(real mu, real sigma, real ligam, real tau){

#      real A;
#      real B;
#      real C;
#      real p;
#      real hi;
#      real nui;
#      real mui=0;
#      real Up = uniform_rng(.5, 1.0);

#      real gam = ligam;
#      p = (gam < 0) + (tau - (gam < 0))/(2*Phi_approx(-fabs(gam))*exp(.5*pow(gam, 2)));
#      A = (1 - 2*p)/(p - pow(p,2));
#      B = 2/(p - pow(p,2));
#      C = 1/((gam > 0) - p);

#       hi = sigma * inv_Phi(Up);
#      nui = sigma * exponential_rng(1);
#      mui += mu + A * nui + C * fabs(gam) * hi;

#      return normal_rng(mui, sqrt(sigma*B*nui));
#   }
#   "
#   }

# #--- Now define all of these here
#  bprior <- prior("multi_normal(M, V)", class = "b") + set_prior("normal(0, tau)", class = "b") + set_prior("target += normal_lpdf(tau | 0, 10)", check = FALSE)

#  SigbPrio <- bdiag(rep(10, 9), as.positive.definite(as.inverse(P.mat(20) + diag(rep(.01, 20))))*.1 )

#  SigbPrio1 <- diag(rep(10, 9))
#  SigbPrio1 <- as.positive(as.inverse(Pmat(20) + diag(rep(.01, 20))))


#  stanvars <- stanvar(rep(0, 29), "M", scode = "  vector[29] M;") +
#    stanvar(SigbPrio1, "V1", scode = "  matrix[9, 9] V1;") + stanvar(SigbPrio1, "V2", scode = "  matrix[20, 20] V2;") + stanvar(scode = "real<lower=0> tau;",block = "parameters")

#  stanvar(scode = "  matrix[29, 29] V = ;")

#  #--- Assume we know tau
#  bprior <- prior(multi_normal(M, V), class = "b")
#  Tp02 = as.positive.definite(as.inverse(P.mat(20) + diag(rep(.01, 20))))*.1

#  SigbPrio <- as.matrix(1. * bdiag( diag(rep(10, 7)) , Tp02))  # 1.*as.symmetric.matrix

#  dim(SigbPrio)
#  # stanvars <- stanvar(rep(0, 29), "M", scode = "  vector[29] M;") +
# #   stanvar(SigbPrio, "V", scode = "  matrix[29, 29] V;")
#  V <- SigbPrio
#  M <- rep(0, 29)

#  stanvars2 <- stanvar(scode = stan_funs2, block = "functions") + stanvar(rep(0, 27), "M", scode = "  vector[27] M;") + stanvar(1*SigbPrio, "V", scode = " matrix[27, 27] V;")

#  # Tp <- df_sampleWt
#  # Tp$V <- I(V)

#  fit_fakeqGal3a.5 <- brms::brm(bf(as.formula(form), tau = tau0), data = df_sampleWt, family = GAL2, stanvars = stanvars2, prior = bprior,
#                                chains = 2,iter = 500, control = list(adapt_delta = 0.97), data2 = list(M = rep(0, 29), V = SigbPrio)) # init = 0.1,


#  make_stancode(bf(as.formula(form), tau = tau0), data = df_sampleWt, family = GAL2, stanvars = stanvars2, prior = bprior,
#                chains = 2,iter = 500, control = list(adapt_delta = 0.97)) #, data2 = list(M = rep(0, 29), V = SigbPrio))


#  make_stancode(bf(as.formula(form), tau = tau0), data = Dt0, family = GAL2, stanvars = stanvars2, prior = bprior,
#                chains = 2,iter = 500, control = list(adapt_delta = 0.97), data2 = list(M = rep(0, 29), V = SigbPrio))


# Dt0 <- make_standata(bf(as.formula(form), tau = tau0), data = df_sampleWt, family = GAL2, stanvars = stanvars2, prior = bprior,
#                      chains = 2,iter = 500, control = list(adapt_delta = 0.97), data2 = list(M = rep(0, 29), V = SigbPrio))
# Dt0$pn = 20
# Dt0$Xold <- Dt0$X
# Dt0$X <- Dt0$Xold[,!grepl("W.[0-9]", colnames(Dt0$Xold))]
# Dt0$Wn <- Dt0$Xold[,grepl("W.[0-9]", colnames(Dt0$Xold))]; dim(Dt0$Wn)
# Dt0$M <- numeric(Dt0$pn)
# Dt0$V <- as.positive.definite(as.inverse(P.mat(20) + diag(rep(.01, 20))))*.1
# Dt0$K <- ncol(Dt0$X)
# Dt0$Bd <- Bd*.99
# Dt0$tau0 <- tau0

# #--- Fit the model
# library(rstan)

# {
# StandCodeV1 <- "
# // > fit_fakeqGal3a.5$model
# // generated with brms 2.18.0
# functions {

#   /*
#   A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + log1m_exp(fabs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2))));
#   gam = (gamU - gamL) * ligam + gamL;
#   real gam = (gamU - gamL) * ligam + gamL;
#   real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau, real gamL, real gamU){
#   real GAL2_rng(real mu, real sigma, real ligam, real tau, real gamL, real gamU){
#    */
#      /* helper function for asym_laplace_lpdf
#   * Args:
#     *   y: the response value
#   *   tau: quantile parameter in (0, 1)
#   */
#     real rho_quantile(real y, real tau) {
#       if (y < 0) {
#         return y * (tau - 1);
#       } else {
#         return y * tau;
#       }
#     }
#   /* asymmetric laplace log-PDF for a single response
#   * Args:
#     *   y: the response value
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   tau: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lpdf(real y, real mu, real sigma, real tau) {
#       return log(tau * (1 - tau)) -
#         log(sigma) -
#         rho_quantile((y - mu) / sigma, tau);
#     }
#   /* asymmetric laplace log-CDF for a single quantile
#   * Args:
#     *   y: a quantile
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   tau: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lcdf(real y, real mu, real sigma, real tau) {
#       if (y < mu) {
#         return log(tau) + (1 - tau) * (y - mu) / sigma;
#       } else {
#         return log1m((1 - tau) * exp(-tau * (y - mu) / sigma));
#       }
#     }
#   /* asymmetric laplace log-CCDF for a single quantile
#   * Args:
#     *   y: a quantile
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   tau: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lccdf(real y, real mu, real sigma, real tau) {
#       if (y < mu) {
#         return log1m(tau * exp((1 - tau) * (y - mu) / sigma));
#       } else {
#         return log1m(tau) - tau * (y - mu) / sigma;
#       }
#     }

#    real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau){

#    real p_pos;
#    real p_neg;
#    real a3;
#    real a2;
#    real p;
#    real est;
#    real A;
#    real B;
#    real Res = 0;
#    real gam = ligam;
#     p = (gam < 0) + (tau - (gam < 0))/(2*Phi(-fabs(gam))*exp(.5*pow(gam, 2)));
#     p_pos = p -  (gam > 0);
#     p_neg = p -  (gam < 0);
#     est = (y - mu) / sigma;

#     if(fabs(gam) > 0){
#     a3 = p_pos * (est / fabs(gam));
#     a2 = fabs(gam) * (p_neg / p_pos);


#     if(est/gam > 0){
#       A =  0.5*pow(gam, 2)*pow(p_neg/p_pos, 2) - est*p_neg + log_diff_exp(log(Phi_approx(a2-a3)), log(Phi_approx(a2)) );
#       B =  0.5*pow(gam, 2) - p_pos*est + log(Phi_approx(-fabs(gam) + a3));
#       Res = log(2*p*(1-p)) - log(sigma) +  log_sum_exp(A, B);
#     }else{
#       Res =  log(2*p*(1-p)) - log(sigma) - p_pos * est + 0.5 * pow(gam, 2) + log(Phi_approx(-fabs(gam) ));
#     }
#     }else{
#     Res = asym_laplace_lpdf( y | mu, sigma, tau);
#     }

#     return Res;
#    }

#   real GAL2_rng(real mu, real sigma, real ligam, real tau){

#      real A;
#      real B;
#      real C;
#      real p;
#      real hi;
#      real nui;
#      real mui=0;
#      real Up = uniform_rng(.5, 1.0);

#      real gam = ligam;
#      p = (gam < 0) + (tau - (gam < 0))/(2*Phi_approx(-fabs(gam))*exp(.5*pow(gam, 2)));
#      A = (1 - 2*p)/(p - pow(p,2));
#      B = 2/(p - pow(p,2));
#      C = 1/((gam > 0) - p);

#       hi = sigma * inv_Phi(Up);
#      nui = sigma * exponential_rng(1);
#      mui += mu + A * nui + C * fabs(gam) * hi;

#      return normal_rng(mui, sqrt(sigma*B*nui));
#   }

# }
# data {
# int<lower=1> pn; //Number of basis;
#   int<lower=1> N;  // total number of observations
#   vector[N] Y;  // response variable
#   int<lower=1> K;  // number of population-level effects
#   matrix[N, K] X;  // population-level design matrix
#   // data for splines
#   int Ks;  // number of linear effects
#   matrix[N, Ks] Xs;  // design matrix for the linear effects

#   int nb_1;  // number of bases
#   int knots_1[nb_1];  // number of knots
#   // basis function matrices
#   matrix[N, knots_1[1]] Zs_1_1;

#   int nb_2;  // number of bases
#   int knots_2[nb_2];  // number of knots
#   // basis function matrices
#   matrix[N, knots_2[1]] Zs_2_1;
#   int prior_only;  // should the likelihood be ignored?
#   matrix[N,pn] Wn;
#     vector[pn] M;
#    matrix[pn, pn] V;
#    vector[2] Bd;
#    real<lower=0,upper=1> tau0;
# }
# transformed data {
#   int Kc = K - 1;
#   matrix[N, Kc] Xc;  // centered version of X without an intercept
#   vector[Kc] means_X;  // column means of X before centering
#   for (i in 2:K) {
#     means_X[i - 1] = mean(X[, i]);
#     Xc[, i - 1] = X[, i] - means_X[i - 1];
#   }
# }
# parameters {
#   vector[Kc] b;  // population-level effects
#   real Intercept;  // temporary intercept for centered predictors
#   vector[Ks] bs;  // spline coefficients
#   vector[pn] bw;

#   // standarized spline coefficients
#   vector[knots_1[1]] zs_1_1;
#   real<lower=0> sds_1_1;  // standard deviations of spline coefficients

#   // standarized spline coefficients
#   vector[knots_2[1]] zs_2_1;
#   real<lower=0> sds_2_1;  // standard deviations of spline coefficients
#   real<lower=0> sigma;  // dispersion parameter
#   real<lower=Bd[1],upper=Bd[2]> ligam;
#   real<lower=0> teta;
# }
# transformed parameters {
#   // actual spline coefficients
#   vector[knots_1[1]] s_1_1;
#   // actual spline coefficients
#   vector[knots_2[1]] s_2_1;
#   //matrix[pn,pn] V1;
#   real tau = tau0;
#   real lprior = 0;  // prior contributions to the log posterior
#   // compute actual spline coefficients
#   s_1_1 = sds_1_1 * zs_1_1;
#   //V1 = teta*V;
#   // compute actual spline coefficients
#   s_2_1 = sds_2_1 * zs_2_1;
#   lprior += normal_lpdf(b | 0, 10);
#   lprior += multi_normal_lpdf(bw | M, teta*V);
#   //lprior += cauchy_lpdf(teta| 0, 1) - 1* cauchy_lccdf(0| 0, 1);
#   lprior += gamma_lpdf(teta| 1, 1);
#   lprior += student_t_lpdf(Intercept | 3, 28.7, 5.6);
#   lprior += normal_lpdf(bs | 0, 10);
#   lprior += student_t_lpdf(sds_1_1 | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
#   lprior += student_t_lpdf(sds_2_1 | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
#   lprior += student_t_lpdf(sigma | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
# }
# model {
#   // likelihood including constants
#   if (!prior_only) {
#     // initialize linear predictor term
#     vector[N] mu = rep_vector(0.0, N);
#     mu += Intercept + Xc * b + Wn*bw + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;
#     for (n in 1:N) {
#       target += GAL2_lpdf(Y[n] | mu[n], sigma, ligam, tau);
#     }
#   }
#   // priors including constants
#   target += lprior;
#   target += std_normal_lpdf(zs_1_1);
#   target += std_normal_lpdf(zs_2_1);
# }
# generated quantities {
#   // actual population-level intercept
#   real b_Intercept = Intercept - dot_product(means_X, b);
#   //Vector[N] ASTP = Xs[1:N,1] * bs[1] + Zs_1_1 * s_1_1;
#   //vector[N] MVPA = Xs[1:N,2] * bs[2] + Zs_2_1 * s_2_1;
# }
# "
# }

# #--- Translate the Stan model (compile the model)
# #modCopile <- stan_model(model_code = StandCodeV1, verbose = TRUE)


# fit <- stan(model_code = StandCodeV1,
#             data = Dt0, iter = 10000, chains = 2, verbose = TRUE,control = list(adapt_delta = 0.97))

# bs2 = bs(a, df = pn, degree = 3, intercept = T)

# PostAnaly <- function(fit, Bx=bs2, Dt0){

#  PostMat <- as.matrix(fit)
#   b <- PostMat[, c("Intercept","b_Intercept", paste0("b[",1:7,"]"))]
#  bw <- PostMat[, grepl("bw", colnames(PostMat))] #
#  bs <- PostMat[, grepl("bs", colnames(PostMat))] #
#  bs11 <- PostMat[, grepl("bs", colnames(PostMat))] #
#    #---------------------------------------------------------------------------
#    #--- Plot the estimates
#    #---------------------------------------------------------------------------
#     Fw <- bw %*% t(Bx)
#  fAstp <- tcrossprod(bs[,1], Dt0$Xs[,1]) + tcrossprod(PostMat[, grepl("^*s_1_1", colnames(PostMat))], Dt0$Zs_1_1)
#  fMvpa <- tcrossprod(bs[,2], Dt0$Xs[,2]) + tcrossprod(PostMat[, grepl("^*s_2_1", colnames(PostMat))], Dt0$Zs_2_1)

#  list(Fw =Fw, fAstp = fAstp, fMvpa = fMvpa)

#  matplot(Dt0$Xs[,1], t(rbind(apply(fAstp,2,quantile, probs=c(0.025, 0.975)), colMeans(fAstp)) ), type="p", lwd=3)
#  matplot(Dt0$Xs[,2], t(rbind(apply(fMvpa,2,quantile, probs=c(0.025, 0.975)), colMeans(fMvpa)) ), type="p", lwd=3)

#  matplot(seq(0,1,length = 1440), t(rbind(apply(Fw,2,quantile, probs=c(0.025, 0.975)), colMeans(Fw))), type="p", lwd=3)

# }


# #(bf(as.formula(form), quantile = .5), data = df_sampleWt, family = asym_laplace())

# #stanvars <- stanvar(scode = stan_funs, block = "functions")
# #stanvars2 <- stanvar(scode = stan_funs2, block = "functions")

# #--- Fit the quantile regression
# #-- quantile

# #fit_fakeqGal3a.5 <- brms::brm(bf(as.formula(form), tau = tau0), data = df_sampleWt, family = GAL2, stanvars = stanvars2,chains = 2,iter = 500, control = list(adapt_delta = 0.97)) # init = 0.1,




























# #Originial model
# {

#   Code <- "
# // generated with brms 2.18.0
# data {
#   int<lower=1> N;  // total number of observations
#   vector[N] Y;  // response variable
#   int<lower=1> K;  // number of population-level effects
#   matrix[N, K] X;  // population-level design matrix
#   // data for splines
#   int Ks;  // number of linear effects
#   matrix[N, Ks] Xs;  // design matrix for the linear effects
#   // data for spline s(MnASTP, bs = "ps", k = 11)
#   int nb_1;  // number of bases
#   int knots_1[nb_1];  // number of knots
#   // basis function matrices
#   matrix[N, knots_1[1]] Zs_1_1;
#   // data for spline s(MnMVPA, bs = "ps", k = 11)
#   int nb_2;  // number of bases
#   int knots_2[nb_2];  // number of knots
#   // basis function matrices
#   matrix[N, knots_2[1]] Zs_2_1;
#   int prior_only;  // should the likelihood be ignored?
#   // This part is the part we add to the model related to W(t)
#   //matrix[N,]
# }
# transformed data {
#   int Kc = K - 1;
#   matrix[N, Kc] Xc;  // centered version of X without an intercept
#   vector[Kc] means_X;  // column means of X before centering
#   for (i in 2:K) {
#     means_X[i - 1] = mean(X[, i]);
#     Xc[, i - 1] = X[, i] - means_X[i - 1];
#   }
# }
# parameters {
#   vector[Kc] b;  // population-level effects
#   real Intercept;  // temporary intercept for centered predictors
#   vector[Ks] bs;  // spline coefficients
#   // parameters for spline s(MnASTP, bs = "ps", k = 11)
#   // standarized spline coefficients
#   vector[knots_1[1]] zs_1_1;
#   real<lower=0> sds_1_1;  // standard deviations of spline coefficients
#   // parameters for spline s(MnMVPA, bs = "ps", k = 11)
#   // standarized spline coefficients
#   vector[knots_2[1]] zs_2_1;
#   real<lower=0> sds_2_1;  // standard deviations of spline coefficients
#   real<lower=0> sigma;  // dispersion parameter
# }
# transformed parameters {
#   // actual spline coefficients
#   vector[knots_1[1]] s_1_1;
#   // actual spline coefficients
#   vector[knots_2[1]] s_2_1;
#   real lprior = 0;  // prior contributions to the log posterior
#   // compute actual spline coefficients
#   s_1_1 = sds_1_1 * zs_1_1;
#   // compute actual spline coefficients
#   s_2_1 = sds_2_1 * zs_2_1;
#   lprior += student_t_lpdf(Intercept | 3, 28.7, 5.6);
#   lprior += student_t_lpdf(sds_1_1 | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
#   lprior += student_t_lpdf(sds_2_1 | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
#   lprior += student_t_lpdf(sigma | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
# }
# model {
#   // likelihood including constants
#   if (!prior_only) {
#     // initialize linear predictor term
#     vector[N] mu = rep_vector(0.0, N);
#     mu += Intercept + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;
#     target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);
#   }
#   // priors including constants
#   target += lprior;
#   target += std_normal_lpdf(zs_1_1);
#   target += std_normal_lpdf(zs_2_1);
# }
# generated quantities {
#   // actual population-level intercept
#   real b_Intercept = Intercept - dot_product(means_X, b);
# }


# "


# }



# #-- second approach
#   df_sampleWt <- data.frame(df_sample, W = Wt)
#   #colnames(df_sampleWt) <- c(colnames(df_sample), paste0("W", 1:ncol(Wt)))

#   form <- paste0("Y ~ ", paste(c("Gender","Race","HealthCondt2","AgeYR", colnames(df_sampleWt)[grepl("W",colnames(df_sampleWt))]), collapse = " + "))
#   for(pa in varX){
#     form <- paste0(form, " + s(", pa, ", bs='ps', k=11)")
#   }



#   fit_fakeB <- brms::brm(as.formula(form), data = df_sampleWt, family = gaussian())

#   fit_fakeB$model
#   conditional_smooths(fit_fakeB)

#   #-- Change the dataset
#   sdat13 <- make_standata(form, data = df_sampleWt, family = gaussian())
#   sdat13$Kc = sdat13$K - 1;
#   sdat13$Xc = matrix(nrow=sdat13$N,ncol=sdat13$Kc);  #// centered version of X without an intercept
#   sdat13$means_X <- numeric(sdat13$Kc);  #// column means of X before centering
#   for (i in 2:K) {
#     sdat13$means_X[i - 1] = mean(X[, i]);
#     sdat13$Xc[, i - 1] = sdat13$X[, i] - sdat13$means_X[i - 1];
#   }


#   sdat13$gamlB = Bd[1]
#   sdat13$gamuB = Bd[2]
#   sdat13$tau0 = tau0

#  #--- Gaussian model with no smoothing
#   {
#     StanCode <- "

#     >   fit_fakeB$model
# // generated with brms 2.18.0
# functions {
# }
# data {
#   int<lower=1> N;  // total number of observations
#   vector[N] Y;  // response variable
#   int<lower=1> K;  // number of population-level effects
#   matrix[N, K] X;  // population-level design matrix
#   // data for splines
#   int Ks;  // number of linear effects
#   matrix[N, Ks] Xs;  // design matrix for the linear effects
#   // data for spline s(MnASTP, bs = "ps", k = 11)
#   int nb_1;  // number of bases
#   int knots_1[nb_1];  // number of knots
#   // basis function matrices
#   matrix[N, knots_1[1]] Zs_1_1;
#   // data for spline s(MnMVPA, bs = "ps", k = 11)
#   int nb_2;  // number of bases
#   int knots_2[nb_2];  // number of knots
#   // basis function matrices
#   matrix[N, knots_2[1]] Zs_2_1;
#   int prior_only;  // should the likelihood be ignored?
# }
# transformed data {
#   int Kc = K - 1;
#   matrix[N, Kc] Xc;  // centered version of X without an intercept
#   vector[Kc] means_X;  // column means of X before centering
#   for (i in 2:K) {
#     means_X[i - 1] = mean(X[, i]);
#     Xc[, i - 1] = X[, i] - means_X[i - 1];
#   }
# }
# parameters {
#   vector[Kc] b;  // population-level effects
#   real Intercept;  // temporary intercept for centered predictors
#   vector[Ks] bs;  // spline coefficients
#   // parameters for spline s(MnASTP, bs = "ps", k = 11)
#   // standarized spline coefficients
#   vector[knots_1[1]] zs_1_1;
#   real<lower=0> sds_1_1;  // standard deviations of spline coefficients
#   // parameters for spline s(MnMVPA, bs = "ps", k = 11)
#   // standarized spline coefficients
#   vector[knots_2[1]] zs_2_1;
#   real<lower=0> sds_2_1;  // standard deviations of spline coefficients
#   real<lower=0> sigma;  // dispersion parameter
# }
# transformed parameters {
#   // actual spline coefficients
#   vector[knots_1[1]] s_1_1;
#   // actual spline coefficients
#   vector[knots_2[1]] s_2_1;
#   real lprior = 0;  // prior contributions to the log posterior
#   // compute actual spline coefficients
#   s_1_1 = sds_1_1 * zs_1_1;
#   // compute actual spline coefficients
#   s_2_1 = sds_2_1 * zs_2_1;
#   lprior += student_t_lpdf(Intercept | 3, 28.7, 5.6);
#   lprior += student_t_lpdf(sds_1_1 | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
#   lprior += student_t_lpdf(sds_2_1 | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
#   lprior += student_t_lpdf(sigma | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
# }
# model {
#   // likelihood including constants
#   if (!prior_only) {
#     // initialize linear predictor term
#     vector[N] mu = rep_vector(0.0, N);
#     mu += Intercept + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;
#     target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);
#   }
#   // priors including constants
#   target += lprior;
#   target += std_normal_lpdf(zs_1_1);
#   target += std_normal_lpdf(zs_2_1);
# }
# generated quantities {
#   // actual population-level intercept
#   real b_Intercept = Intercept - dot_product(means_X, b);
# }



#     "





#   }

#  #--- Now we do the
#   #--- Gaussian model with no smoothing
#   {
# Code <- '
#  //   >   fit_fakeB$model
# // generated with brms 2.18.0
# // functions {
# //}
# data {
#    real<lower=0> tau0;
#    real gamlB;
#    real gamuB;
#   int<lower=1> N;  // total number of observations
#   vector[N] Y;  // response variable
#   int<lower=1> K;  // number of population-level effects
#   matrix[N, K] X;  // population-level design matrix
#   // data for splines
#   int Ks;  // number of linear effects
#   matrix[N, Ks] Xs;  // design matrix for the linear effects
#   //data for spline s(MnASTP, bs = "ps", k = 11)
#   int nb_1;  // number of bases
#   int knots_1[nb_1];  // number of knots
#   // basis function matrices
#   matrix[N, knots_1[1]] Zs_1_1;
#   //data for spline s(MnMVPA, bs = "ps", k = 11) //
#   int nb_2;  // number of bases
#   int knots_2[nb_2];  // number of knots
#   // basis function matrices
#   matrix[N, knots_2[1]] Zs_2_1;
#   int prior_only;  // should the likelihood be ignored?
# }
# transformed data {
#   int Kc = K - 1;
#   matrix[N, Kc] Xc;  // centered version of X without an intercept
#   vector[Kc] means_X;  // column means of X before centering
#   for (i in 2:K) {
#     means_X[i - 1] = mean(X[, i]);
#     Xc[, i - 1] = X[, i] - means_X[i - 1];
#   }
# }
# parameters {
#   vector[Kc] b;  // population-level effects
#   real Intercept;  // temporary intercept for centered predictors
#   vector[Ks] bs;  // spline coefficients
#   // parameters for spline s(MnASTP, bs = "ps", k = 11)
#   // standarized spline coefficients
#   vector[knots_1[1]] zs_1_1;
#   real<lower=0> sds_1_1;  // standard deviations of spline coefficients
#   // parameters for spline s(MnMVPA, bs = "ps", k = 11)
#   // standarized spline coefficients
#   vector[knots_2[1]] zs_2_1;
#   real<lower=0> sds_2_1;  // standard deviations of spline coefficients
#   real<lower=0> sigma;  // dispersion parameter
#   real<lower=gamlB,upper=gamuB> gam; // Gamma function
#   vector<lower=0>[N] si;
#   vector<lower=0>[N] wi;
# }
# transformed parameters {
#   // actual spline coefficients
#   vector[knots_1[1]] s_1_1;
#   // actual spline coefficients
#   vector[knots_2[1]] s_2_1;
#   // hi
#   vector<lower=0>[N] hi;
#   // nui
#   vector<lower=0>[N] nui;
#   real<lower=0,upper=1> p;
#   real A;
#   real B;
#   real C;

#   real lprior = 0;  // prior contributions to the log posterior

#   // compute actual spline coefficients
#   s_1_1 = sds_1_1 * zs_1_1;
#   // compute actual spline coefficients
#   s_2_1 = sds_2_1 * zs_2_1;
#   // Compute hi
#   hi = sigma * wi;
#   // compute nui
#   nui = sigma * si;
#   // compute p
#   p = (gam < 0) + (tau0 - (gam<0))/(2*Phi_approx(-abs(gam))*exp(.5*pow(gam,2)));
#   A = (1-2*p)/(p*(1-p));
#   B = 2/(p*(1-p));
#   C = 1/((gam > 0) - p);

#   // add the priors
#   lprior += student_t_lpdf(Intercept | 3, 28.7, 5.6);
#   lprior += student_t_lpdf(sds_1_1 | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
#   lprior += student_t_lpdf(sds_2_1 | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
#   lprior += student_t_lpdf(sigma | 3, 0, 5.6)
#     - 1 * student_t_lccdf(0 | 3, 0, 5.6);
#   lprior += std_normal_lpdf(si) - 1 * normal_lccdf(0| 0.0, 1.); //std_normal_lccdf(0.0);
#   lprior += exponential_lpdf(wi | 1.);
#   lprior += uniform_lpdf(gam | gamlB, gamuB);

# }
# model {
#   // likelihood including constants
#   if (!prior_only) {
#     // initialize linear predictor term
#     vector[N] mu = rep_vector(0.0, N);
#     mu += Intercept + Xc * b + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + A*nui + C*abs(gam)*hi; // linear and basis effect and intercept
#     target += normal_lpdf(Y | mu, sqrt(sigma*B*nui));
#     //for( i in 1:N){
#     //target += normal_lpdf(Y | mu[i], sqrt(sigma*B*nui[i]));
#     //target += normal_glm_lpdf(Y | Xc, mu, b, sigma);
#     //}
#   }
#   // priors including constants
#   target += lprior;
#   target += std_normal_lpdf(zs_1_1);
#   target += std_normal_lpdf(zs_2_1);
# }
# generated quantities {
#   // actual population-level intercept
#   real b_Intercept = Intercept - dot_product(means_X, b);
# }
# ' #), file = 'GalVers1.stan')





#   }

#   ols.dso <- rstan::stan_model(model_code = Code)

#   ols.fit <- rstan::sampling(ols.dso, data = sdat13, chains = 3, iter = 5000)
#   ols.fit


# #--- Margina
# msms <- marginal_smooths(fit_fake)
# msms <- conditional_smooths(fit_fake)
# plot(msms)



# sdata2 <- make_standata(count ~ s(zAge) + zBase * Trt + (1|patient),
#                         data = epilepsy, family = "poisson")
# str(sdata2)

# form <- paste0("Y ~ ", paste(colnames(ZMat)[c(2:3,6,9)], collapse = " + "))
# for(pa in varX){
#   form <- paste0(form, " + s(", pa, ", bs='ps', k=11)")
# }

# sdat13 <- make_standata(form, data = df_sample, family = gaussian())




# #--- Quantile regression
# {
# QCode <- '
# //> fit_fakeq.5$model
# // generated with brms 2.18.0
# functions {
#   /* helper function for asym_laplace_lpdf
#   * Args:
#     *   y: the response value
#   *   quantile: quantile parameter in (0, 1)
#   */
#     real rho_quantile(real y, real quantile) {
#       if (y < 0) {
#         return y * (quantile - 1);
#       } else {
#         return y * quantile;
#       }
#     }
#   /* asymmetric laplace log-PDF for a single response
#   * Args:
#     *   y: the response value
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   quantile: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lpdf(real y, real mu, real sigma, real quantile) {
#       return log(quantile * (1 - quantile)) -
#         log(sigma) -
#         rho_quantile((y - mu) / sigma, quantile);
#     }
#   /* asymmetric laplace log-CDF for a single quantile
#   * Args:
#     *   y: a quantile
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   quantile: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lcdf(real y, real mu, real sigma, real quantile) {
#       if (y < mu) {
#         return log(quantile) + (1 - quantile) * (y - mu) / sigma;
#       } else {
#         return log1m((1 - quantile) * exp(-quantile * (y - mu) / sigma));
#       }
#     }
#   /* asymmetric laplace log-CCDF for a single quantile
#   * Args:
#   *   y: a quantile
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   quantile: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lccdf(real y, real mu, real sigma, real quantile) {
#       if (y < mu) {
#         return log1m(quantile * exp((1 - quantile) * (y - mu) / sigma));
#       } else {
#         return log1m(quantile) - quantile * (y - mu) / sigma;
#       }
#     }
# }
# data {
#   int<lower=1> N;  // total number of observations
#   vector[N] Y;  // response variable
#   int<lower=1> K;  // number of population-level effects
#   matrix[N, K] X;  // population-level design matrix
#   // data for splines
#   int Ks;  // number of linear effects
#   matrix[N, Ks] Xs;  // design matrix for the linear effects
#   // data for spline s(MnASTP, bs = "ps", k = 11)
#   int nb_1;  // number of bases
#   int knots_1[nb_1];  // number of knots
#   // basis function matrices
#   matrix[N, knots_1[1]] Zs_1_1;
#   // data for spline s(MnMVPA, bs = "ps", k = 11)
#   int nb_2;  // number of bases
#   int knots_2[nb_2];  // number of knots
#   // basis function matrices
#   matrix[N, knots_2[1]] Zs_2_1;
#   int prior_only;  // should the likelihood be ignored?
# }
# transformed data {
#   int Kc = K - 1;
#   matrix[N, Kc] Xc;  // centered version of X without an intercept
#   vector[Kc] means_X;  // column means of X before centering
#   for (i in 2:K) {
#     means_X[i - 1] = mean(X[, i]);
#     Xc[, i - 1] = X[, i] - means_X[i - 1];
#   }
# }
# parameters {
#   vector[Kc] b;  // population-level effects
#   real Intercept;  // temporary intercept for centered predictors
#   vector[Ks] bs;  // spline coefficients
#   // parameters for spline s(MnASTP, bs = "ps", k = 11)
#   // standarized spline coefficients
#   vector[knots_1[1]] zs_1_1;
#   real<lower=0> sds_1_1;  // standard deviations of spline coefficients
#   // parameters for spline s(MnMVPA, bs = "ps", k = 11)
#   // standarized spline coefficients
#   vector[knots_2[1]] zs_2_1;
#   real<lower=0> sds_2_1;  // standard deviations of spline coefficients
#   real<lower=0> sigma;  // dispersion parameter
# }
# transformed parameters {
#   // actual spline coefficients
#   vector[knots_1[1]] s_1_1;
#   // actual spline coefficients
#   vector[knots_2[1]] s_2_1;
#   real quantile = 0.5;  // quantile parameter
#   real lprior = 0;  // prior contributions to the log posterior
#   // compute actual spline coefficients
#   s_1_1 = sds_1_1 * zs_1_1;
#   // compute actual spline coefficients
#   s_2_1 = sds_2_1 * zs_2_1;
#   lprior += student_t_lpdf(Intercept | 3, 28.7, 5.6);
#   lprior += student_t_lpdf(sds_1_1 | 3, 0, 5.6)
#   - 1 * student_t_lccdf(0 | 3, 0, 5.6);
#   lprior += student_t_lpdf(sds_2_1 | 3, 0, 5.6)
#   - 1 * student_t_lccdf(0 | 3, 0, 5.6);
#   lprior += student_t_lpdf(sigma | 3, 0, 5.6)
#   - 1 * student_t_lccdf(0 | 3, 0, 5.6);
# }
# model {
#   // likelihood including constants
#   if (!prior_only) {
#     // initialize linear predictor term
#     vector[N] mu = rep_vector(0.0, N);
#     mu += Intercept + Xc * b + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;
#     for (n in 1:N) {
#       target += asym_laplace_lpdf(Y[n] | mu[n], sigma, quantile);
#     }
#   }
#   // priors including constants
#   target += lprior;
#   target += std_normal_lpdf(zs_1_1);
#   target += std_normal_lpdf(zs_2_1);
# }
# generated quantities {
#   // actual population-level intercept
#   real b_Intercept = Intercept - dot_product(means_X, b);
# }
# '
# }

# #--- define custom functions
#   GAL <- custom_family(
#     "GAL", dpars = c("mu","sigma","ligam","tau","gamL","gamU"), links=c("identity","log","identity","identity", rep("identity",2)),
#     lb=c(NA,0, NA,0,NA,NA), ub=c(NA,NA,NA,1,NA,NA), type="real") #, vars = "vint1[n]"

# #--- Stan function
#   {
#     GAL <- custom_family(
#       "GAL", dpars = c("mu", "sigma","ligam", "tau","gamL","gamU"), links=c("identity","log","identity","identity", rep("identity",2)),
#       lb=c(NA,0,-10.,0,NA,NA), ub=c(NA,NA,10.,1,NA,NA), type="real") #, vars = "vint1[n]"

#   stan_funs <- "

#   /*
#   A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + log1m_exp(fabs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2))));
#    */
#    real GAL_lpdf(real y, real mu, real sigma, real ligam, real tau, real gamL, real gamU){

#    real p_pos;
#    real p_neg;
#    real a3;
#    real a2;
#    real p;
#    real est;
#    real A;
#    real B;
#    real Res = 0;
#    real gam = (gamU - gamL) * inv_logit(ligam) + gamL;
#     p = (gam < 0) + (tau - (gam < 0))/(2*Phi(-fabs(gam))*exp(.5*pow(gam, 2)));
#     p_pos = p - 1.*(gam > 0);
#     p_neg = p - 1.*(gam < 0);
#     est = (y - mu) / sigma;

#     a3 = p_pos * est / fabs(gam);
#     a2 = fabs(gam) * p_neg/p_pos;
#     A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log_diff_exp( log(Phi_approx(a2-a3)), log(Phi_approx(a2)) );
#     B = -p_pos*est + 0.5*pow(gam, 2) + log(Phi_approx(-abs(gam) + a3));

#     if(est/gam > 0){
#       Res = log(2*p*(1-p)) - log(sigma) + log_sum_exp(A, B);
#     }else{
#       Res =  log(2*p*(1-p)) - log(sigma) - p_pos*est + 0.5*pow(gam, 2) + log(Phi_approx(-fabs(gam) ));
#       }

#     return Res;
#    }

#   real GAL_rng(real mu, real sigma, real ligam, real tau, real gamL, real gamU){

#      real A;
#      real B;
#      real C;
#      real p;
#      real hi;
#      real nui;
#      real mui=0;
#      real Up = uniform_rng(.5, 1.0);

#      real gam = (gamU - gamL) * inv_logit(ligam) + gamL;
#      p = (gam < 0) + (tau - (gam < 0))/(2*Phi_approx(-fabs(gam))*exp(.5*pow(gam, 2)));
#      A = (1-2*p)/(p - pow(p,2));
#      B = 2/(p - pow(p,2));
#      C = 1/((gam > 0) - p);

#       hi = sigma * inv_Phi(Up);
#      nui = sigma * exponential_rng(1);
#      mui += mu + A * nui + C * fabs(gam) * hi;

#      return normal_rng(mui, sqrt(sigma*B*nui));
#   }
#   "
#   }


#   #--- Stan function (case 2)
#   {

#     #GAL2 <- custom_family(
#     #  "GAL2", dpars = c("mu", "sigma","ligam", "tau","gamL","gamU"), links=c("identity","log","identity","identity", rep("identity",2)),
#     #  lb=c(NA,0, -1.087,0,NA,NA), ub=c(NA,NA, 1.087,1,NA,NA), type="real") #, vars = "vint1[n]"
#     # GAL2 <- custom_family(
#     #   "GAL2", dpars = c("mu", "sigma","ligam"), links= c("identity","log","identity"),
#     #   lb=c(NA,0, -1.087), ub=c(NA,NA, 1.087), type="real", vars = c("tau","gamL","gamU")) #, vars = "vint1[n]"
#     GAL2 <- custom_family(
#       "GAL2", dpars = c("mu", "sigma","ligam", "tau"), links= c("identity","log","identity","identity"),
#       lb=c(NA,0, -1.086,0), ub=c(NA,NA, 1.086,1), type="real") #, vars = c("tau")) #, vars = "vint1[n]"

#     GAL2 <- custom_family(
#       "GAL2", dpars = c("mu", "sigma","ligam", "tau"), links= c("identity","log","logit","identity"),
#       lb=c(NA,0, -15,0), ub=c(NA,NA, 15,1), type="real", vars=c("gamL", "gamU")) #, vars = c("tau")) #, vars = "vint1[n]"

#     GAL2 <- custom_family(
#      "GAL2", dpars = c("mu", "sigma","ligam", "tau","gamL","gamU"), links=c("identity","log","identity","identity", rep("identity",2)),
#      lb=c(NA,0, -15,0,NA,NA), ub=c(NA,NA, 15,1,NA,NA), type="real") #, vars = "vint1[n]"

#     GAL2 <- custom_family(
#       "GAL2", dpars = c("mu", "sigma","ligam", "tau"), links=c("identity","log","identity","identity"),
#       lb=c(NA,0, Bd[1]*.99,0), ub=c(NA,NA, Bd[2]*.99,1), type="real") #, vars = "vint1[n]"


#     stan_funs2 <- "
#   /*
#   A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + log1m_exp(fabs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2))));
#   gam = (gamU - gamL) * ligam + gamL;
#   real gam = (gamU - gamL) * ligam + gamL;
#   real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau, real gamL, real gamU){
#   real GAL2_rng(real mu, real sigma, real ligam, real tau, real gamL, real gamU){
#    */
#      /* helper function for asym_laplace_lpdf
#   * Args:
#     *   y: the response value
#   *   tau: quantile parameter in (0, 1)
#   */
#     real rho_quantile(real y, real tau) {
#       if (y < 0) {
#         return y * (tau - 1);
#       } else {
#         return y * tau;
#       }
#     }
#   /* asymmetric laplace log-PDF for a single response
#   * Args:
#     *   y: the response value
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   tau: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lpdf(real y, real mu, real sigma, real tau) {
#       return log(tau * (1 - tau)) -
#         log(sigma) -
#         rho_quantile((y - mu) / sigma, tau);
#     }
#   /* asymmetric laplace log-CDF for a single quantile
#   * Args:
#     *   y: a quantile
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   tau: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lcdf(real y, real mu, real sigma, real tau) {
#       if (y < mu) {
#         return log(tau) + (1 - tau) * (y - mu) / sigma;
#       } else {
#         return log1m((1 - tau) * exp(-tau * (y - mu) / sigma));
#       }
#     }
#   /* asymmetric laplace log-CCDF for a single quantile
#   * Args:
#     *   y: a quantile
#   *   mu: location parameter
#   *   sigma: positive scale parameter
#   *   tau: quantile parameter in (0, 1)
#   * Returns:
#     *   a scalar to be added to the log posterior
#   */
#     real asym_laplace_lccdf(real y, real mu, real sigma, real tau) {
#       if (y < mu) {
#         return log1m(tau * exp((1 - tau) * (y - mu) / sigma));
#       } else {
#         return log1m(tau) - tau * (y - mu) / sigma;
#       }
#     }

#    real GAL2_lpdf(real y, real mu, real sigma, real ligam, real tau){

#    real p_pos;
#    real p_neg;
#    real a3;
#    real a2;
#    real p;
#    real est;
#    real A;
#    real B;
#    real Res = 0;
#    real gam = ligam;
#     p = (gam < 0) + (tau - (gam < 0))/(2*Phi(-fabs(gam))*exp(.5*pow(gam, 2)));
#     p_pos = p -  (gam > 0);
#     p_neg = p -  (gam < 0);
#     est = (y - mu) / sigma;

#     if(fabs(gam) > 0){
#     a3 = p_pos * (est / fabs(gam));
#     a2 = fabs(gam) * (p_neg / p_pos);


#     if(est/gam > 0){
#       A =  0.5*pow(gam, 2)*pow(p_neg/p_pos, 2) - est*p_neg + log_diff_exp(log(Phi_approx(a2-a3)), log(Phi_approx(a2)) );
#       B =  0.5*pow(gam, 2) - p_pos*est + log(Phi_approx(-fabs(gam) + a3));
#       Res = log(2*p*(1-p)) - log(sigma) +  log_sum_exp(A, B);
#     }else{
#       Res =  log(2*p*(1-p)) - log(sigma) - p_pos * est + 0.5 * pow(gam, 2) + log(Phi_approx(-fabs(gam) ));
#     }
#     }else{
#     Res = asym_laplace_lpdf( y | mu, sigma, tau);
#     }

#     return Res;
#    }

#   real GAL2_rng(real mu, real sigma, real ligam, real tau){

#      real A;
#      real B;
#      real C;
#      real p;
#      real hi;
#      real nui;
#      real mui=0;
#      real Up = uniform_rng(.5, 1.0);

#      real gam = ligam;
#      p = (gam < 0) + (tau - (gam < 0))/(2*Phi_approx(-fabs(gam))*exp(.5*pow(gam, 2)));
#      A = (1 - 2*p)/(p - pow(p,2));
#      B = 2/(p - pow(p,2));
#      C = 1/((gam > 0) - p);

#       hi = sigma * inv_Phi(Up);
#      nui = sigma * exponential_rng(1);
#      mui += mu + A * nui + C * fabs(gam) * hi;

#      return normal_rng(mui, sqrt(sigma*B*nui));
#   }
#   "
#   }



# #--- Now define all of these here
# stanvars <- stanvar(scode = stan_funs, block = "functions")
# stanvars2 <- stanvar(scode = stan_funs2, block = "functions")

# #--- Fit the quantile regression
# #-- quantile

# #fit_fakeq.5 <- brms::brm(bf(as.formula(form), quantile = .5), data = df_sampleWt, family = asym_laplace(), prior = set_prior(horseshoe(df=1,par_ratio = 0.1), class="b"))
# #fit_fakeq.5 <- brms::brm(bf(as.formula(form), quantile = .5), data = df_sampleWt, family = asym_laplace())

# # fit_fakeqGal.5 <- brms::brm(bf(as.formula(form), tau = .5, gamL = Bd[1], gamU = Bd[2]), data = df_sampleWt, family = GAL,
# #                             stanvars = stanvars,chains = 3,iter = 10000, control = list(adapt_delta = 0.97, max_treedepth = 10))


# #fit_fakeqGal2.5 <- brms::brm(bf(as.formula(form), tau = .5, gamL = Bd[1], gamU = Bd[2]), data = df_sampleWt, family = GAL2, stanvars = stanvars2,chains = 3,iter = 1000)
# #fit_fakeqGal2.5 <- brms::brm(bf(as.formula(form), tau = .5,gamL = Bd[1], gamU = Bd[2]), data = df_sampleWt, family = GAL2, stanvars = stanvars2,chains = 2,iter = 2000, init = 0.1, control = list(adapt_delta = 0.97))

# fit_fakeqGal2.5 <- brms::brm(bf(as.formula(form), tau = .5), data = df_sampleWt, family = GAL2, stanvars = stanvars2,chains = 2,iter = 5000, control = list(adapt_delta = 0.97)) # init = 0.1,



# brms::conditional_smooths(fit_fakeqGal2.5)
# #fit_fakeq.5


# #--- Old function

# stan_funsOld <- "
#    real GAL_lpmf(real y, real mu, real sigma, real gam, real tau){

#    real p_pos;
#    real p_neg;
#    real a3;
#    real a2;
#    real p;
#    real est;
#    real A;
#    real B;
#     p = (gam < 0) + (tau - (gam < 0))/(2*Phi_approx(-fabs(gam))*exp(.5*pow(gam, 2)));
#     p_pos = p - 1.*(gam > 0);
#     p_neg = p - 1.*(gam < 0);
#     est = y/sigma;

#     a3 = p_pos*est/fabs(gam);
#     a2 = fabs(gam)*p_neg/p_pos;
#     A = -est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + log1m_exp(fabs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2)))); // Rmpfr::log1mexp(abs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2))) )
#     B = -p_pos*est + 0.5*pow(gam, 2) + log(Phi_approx(-abs(gam) + a3));

#     if(est/gam > 0){
#       #(2*p*(1-p)/sqrt(sigma))*((log(Phi_approx(q = a2-a3) - log(Phi_approx(q = a2))*exp(-est*p_neg + .5*(gam^2)*(p_neg/p_pos)^2) + log(Phi_approx(q = -abs(gam) + a3)*exp(-p_pos*est + 0.5*(gam^2)))
#       (2*p*(1-p)/sigma)*(exp(-est*p_neg + .5*pow(gam, 2)*pow(p_neg/p_pos, 2) + log(Phi_approx(a2-a3)) + Rmpfr::log1mexp(abs(log(Phi_approx(a2-a3)) - log(Phi_approx(a2))) )) + exp(-p_pos*est + 0.5*pow(gam, 2) + log(Phi_approx(-abs(gam) + a3)) )
#       log(2*p*(1-p)) -log(sigma) + log_sum_exp(A, B);
#     }else{
#       log(2*p*(1-p)) -log(sigma) - p_pos*est + 0.5*pow(gam, 2) + log(Phi_approx(-fabs(gam)))
#     }
#    }

#   "
