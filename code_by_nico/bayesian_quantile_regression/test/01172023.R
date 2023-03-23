library(brms)

# --------
# simulating data
# --------

# Simulates data such that q_p(y | x) = beta * x + qnorm(p, 0, 1) * x

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
        )
        # append them to data
        for (j in 1:n_i) {
            x <- append(x, x_i)
            y <- append(y, y_i[j])
        }
    }
}

simulated_data <- data.frame(x = x, y = y)



#------------
# fitting
#-------------


gal <- custom_family("gal",
    dpars = c("mu", "sigma", "gamma"),
    lb = c(NA, 0, NA), type = "real", vars = "vreal1[1]"
)

stan_funs <- "
    real gal_lpdf(real y, real mu, real sigma, real gamma, real p_0){
        real g;
        real p;
        real p_gamma_plus;
        real p_gamma_minus;
        real y_star;
        real aux_1;
        real aux_2;
        real aux_3;
        real aux_4;
        real aux_5;
        real aux_6;
        real aux_7;
        real aux_8;
        real indic_1;
        real indic_2;
        real indic_3;
        g = 2*Phi(-fabs(gamma))*exp(gamma^2/2);
        indic_1 = 0;
        if(gamma > 0){
            indic_1 = 1;
        }
        indic_2 = 0;
        if(gamma < 0){
            indic_2 = 1;
        }
        p = indic_2 + (p_0 - indic_2)/g;
        p_gamma_plus = p - indic_1;
        p_gamma_minus = p - indic_2;
        y_star = (y - mu)/sigma;
        indic_3 = 0;
        if(y_star/gamma > 0){
            indic_3 = 1;
        }
        aux_1 = p_gamma_minus*fabs(gamma)/p_gamma_plus;
        aux_2 = y_star*p_gamma_plus/fabs(gamma);
        aux_3 = - aux_2 + aux_1;
        aux_4 = -y_star*p_gamma_minus + gamma^2/2*(p_gamma_minus/p_gamma_plus)^2;
        aux_5 = -fabs(gamma) + aux_2*indic_3;
        aux_6 = -y_star*p_gamma_plus + gamma^2/2;
        aux_7 = (Phi(aux_3) - Phi(aux_1))*exp(aux_4)*indic_3;
        aux_8 = Phi(aux_5)*exp(aux_6);
        return 2*p*(1 - p)*(aux_7 + aux_8);
    }
"
stanvars <- stanvar(scode = stan_funs, block = "functions")

fit <- brm(y | vreal(0.25) ~ x,
    data = simulated_data,
    family = gal, stanvars = stanvars
)

summary(fit)
