library(nimble)
library(splines)

# Parameters of the model
M <- 10
D <- 3
n <- 100
lK <- M + D
a <- 0.1
b <- 0.1
a_sigma <- 0.1
b_sigma <- 0.1

# Creating synthetic data to estimate
synthetic_function <- function(x) {
    (x - .5)^2
}
grid_x <- seq(0, 1, 1 / (n - 1))
epsilon <- rnorm(n, sd = 0.003)
y <- synthetic_function(grid_x) + epsilon

# Creating the design matrix, X
X <- bs(grid_x, df = lK)


# Creating the penalty matrix, K
# D_2 <- matrix(nrow = lK - 2, ncol = lK)
# for (i in 1:(lK - 2)) {
#     D_2[i, ] <- c(rep(0, i - 1), 1, -2, 1, rep(0, lK - i - 2))
# }
# K <- t(D_2) %*% D_2

# Creating the prior for gamma_coef in Nimble
# dimproperprior <- nimbleFunction(
#     run = function(x = double(lK), tau_2 = double(0), K = nimMatrix(nrow = lK, ncol = lK), log = integer(0, default = 0)) {
#         returnType(double(0))
#         logProb <- -((lK - 2) / 2) * log(tau_2) -
#             (1 / (2 * tau_2)) * t(x) %*% K %*% x
#         if (log) {
#             return(logProb)
#         } else {
#             return(exp(logProb))
#         }
#     }
# )


# Coding the model in Nimble
splinesCode <- nimbleCode({
    aux_cov[1:n, 1:n] <- sigma_2 * diag(n)[1:n, 1:n]
    zeroes_n[1:n] <- rep(0, n)[1:n]
    epsilon_vector[1:n] ~ dmnorm(mean = zeroes_n[1:n], cov = aux_cov[1:n, 1:n])
    epsilon[1:n, 1] <- nimMatrix(epsilon_vector[1:n], ncol = 1)
    y_vec[1:n, 1] <- X[1:n, 1:lK] %*% gamma_coef[1:lK, 1] + epsilon[1:n, 1]
    sigma_2 ~ dinvgamma(a_sigma, b_sigma)

    # fake prior
    # zeroes_lK[1:lK] <- rep(0, lK)[1:lK]
    # I_lK[1:lK, 1:lK] <- diag(lK)[1:lK, 1:lK]
    # # gamma_coef_vector[1:lK] ~ dmnorm(zeroes_lK[1:lK], I_lK[1:lK, 1:lK])
    # gamma_coef_vector[1:lK] ~ dimproperprior(tau_2, K)

    prec <- 1 / tau_2

    for (i in 2:lK) {
        gamma_coef_vector[i] ~ dnorm(
            mean = gamma_coef_vector[i - 1],
            tau = prec
        )
    }
    # gamma_coef_vector[1] ~ dnorm(mean = 0, tau = 1)
    gamma_coef_vector[1] ~ dnorm(mean = 0, tau = 1)


    gamma_coef[1:lK, 1] <- nimMatrix(gamma_coef_vector[1:lK], ncol = 1)
    tau_2 ~ dinvgamma(a, b)
    # tau_2 ~ dnorm(mean = 0, tau = 1)
})

splinesConsts <- list(
    lK = lK, n = n, a = a, b = b, a_sigma = a_sigma,
    b_sigma = b_sigma, X = X
    # K = K
)

splinesData <- list(y_vec = matrix(y, ncol = 1))

splinesInit <- list(
    # sigma_2 = rinvgamma(1, a_sigma, b_sigma),
    # tau_2 = rinvgamma(1, a, b),
    sigma_2 = 0.01,
    tau_2 = 0.01,
    gamma_coef_vector = rep(0, lK),
    # gamma_coef_vector[1] = 0,
    epsilon_vector = rep(0, n), epsilon = matrix(rep(0, n), ncol = 1)
)

splinesModel <- nimbleModel(
    code = splinesCode, name = "splinesModel",
    constants = splinesConsts, data = splinesData,
    inits = splinesInit
)

Csplines <- compileNimble(splinesModel)

splinesConf <- configureMCMC(splinesModel, print = TRUE)

splinesConf$addMonitors(c("sigma_2", "tau_2", "gamma_coef_vector"))

splinesMCMC <- buildMCMC(splinesConf)

CsplinesMCMC <- compileNimble(splinesMCMC, project = splinesModel)

niter <- 1000
set.seed(1)
samples <- runMCMC(CsplinesMCMC, niter = niter)

# par(mfrow = c(3, 1))
# plot(samples[1:1000, "sigma_2"], type="l")
# plot(samples[1:1000, "tau_2"], type="l")
# plot(samples[1:1000, "gamma_coef_vector[1]"], type="l")

# splines$plotGraph()


# mcmc.out <- nimbleMCMC(
#     code = splinesCode, constants = splinesConsts,
#     data = splinesData, inits = splinesInit, summary = TRUE,
#     monitors = c("gamma_coef_vector")
#     # monitors = c("sigma_2")
# )

# mcmc.out$summary




mean_gamma_coef <- vector(length = lK)

names_vectors <- paste0("gamma_coef_vector[", 1:lK, "]")

for (i in 1:lK) {
    mean_gamma_coef[i] <- mean(samples[, names_vectors[i]])
}

y_hat <- X %*% matrix(mean_gamma_coef, ncol = 1)

plot(y)
lines(y_hat)


# plot(y)
# for (i in 100:100) {
#     lines(X %*% matrix(samples[i, names_vectors], ncol = 1))
# }

# plot(X %*% matrix(samples[500, names_vectors], ncol = 1))
