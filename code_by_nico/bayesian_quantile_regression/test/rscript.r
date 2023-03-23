library(invgamma)
library(rstan)

set.seed(1234567890)

# ------
# Generating synthetic data
# -----

# hyperparameters
alpha_shape <- 2
beta_shape <- 5
beta_0 <- 1
B_0 <- 2
n_0 <- 1
d_0 <- 1
delta_0 <- 0.5
D_0 <- 1
p_0 <- 0.5
mu_data <- 0
sigma_data <- 1
J <- 4
c <- 1
n <- 1000

# choosing the gamma parameter
L <- -sqrt(1 - p_0)
U <- sqrt(p_0)
gamma_par <- rbeta(1, alpha_shape, beta_shape) * (U - L) + L

# Constructing the constants A, B and C
g <- function(x) 2 * pnorm(-abs(x)) * exp(x^2 / 2)
indic <- ifelse(gamma_par < 0, 1, 0)
p <- indic + (p_0 - indic) / g(gamma_par)
A <- (1 - 2 * p) / (p * (1 - p))
B <- 2 / (p * (1 - p))
C <- (1 - indic - p)^(-1)

# choosing other parameters
sigma_par_2 <- rinvgamma(1, n_0 / 2, d_0 / 2)
sigma_par <- sqrt(sigma_par_2)
beta_par <- rnorm(1, beta_0, B_0)
delta_par <- rnorm(1, delta_0, D_0)

# constructing the cutoffs
xi_2 <- c
xi_3 <- xi_2 + exp(delta_par)

# choosing values of the covariates
x <- rnorm(n, mu_data, sigma_data)

# observing values for the other random variables
h_rvar <- abs(rnorm(n, 0, sigma_par))
nu_rvar <- rexp(n, sigma_par)
u_rvar <- rnorm(n)

# deducing values of the latent random variable
z <- x * beta_par + A * nu_rvar + C * abs(gamma_par) * h_rvar +
    sqrt(sigma_par * B * nu_rvar) * u_rvar

# deducing values of the response
y <- cut(z, breaks = c(-Inf, 0, xi_2, xi_3, Inf), labels = FALSE)



# ---
# Fitting the model
# ---
fbqror_data <- list(
    x = x, y = y, p_0 = p_0, n = n, c = c, beta_0 = beta_0,
    B_0 = B_0, d_0 = d_0, delta_0 = delta_0, D_0 = D_0,
    shape_1 = alpha_shape, shape_2 = beta_shape, n_0 = n_0
)
inits <- function() {
    list(p = 0.5, nu_rvar = rep(1, n))
}
fit <- stan(
    file = "C:/bfda/code/fbqror_test/stanCode.stan",
    data = fbqror_data, init = inits, chains = 1, iter = 10000
)
print(fit)
plot(fit)
fit@sim
traceplot(fit)
