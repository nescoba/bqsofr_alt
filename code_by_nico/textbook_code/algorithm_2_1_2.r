library(splines)
library(nimble)

# creating grid of x values
X <- seq(-5, 5, .1)

# creating spline basis
B <- t(bs(x = X, degree = D, knots = seq(-5, 5, 1), intercept = TRUE))

# parameters
num_data <- length(X)
num_basis <- nrow(B)

# creating synthetic data
a0 <- 0.2
a <- rnorm(num_basis, 0, 1)
Y_true <- as.vector(a0 * X + a %*% B)
Y <- Y_true + rnorm(num_data, 0, 0.2)


# creating cauchy distribution
# dcauchy <- nimbleFunction(
#     run = function(x = double(0), log = integer(0, default = 0)) {
#         returnType(double(0))
#         logProb <- -log(pi) - log(1 + x^2)
#         if (log) {
#             return(logProb)
#         } else {
#             return(exp(logProb))
#         }
#     }
# )

# creating the model
splinesCode <- nimbleCode({
    a[1:num_basis] <- tau * a_raw[1:num_basis]
    Y_hat[1:num_data] <- a0 * X[1:num_data] +
        a[1:num_basis] %*% B[1:num_basis, 1:num_data]

    I_num_basis[1:num_basis, 1:num_basis] <-
        diag(num_basis)[1:num_basis, 1:num_basis]
    zeroes_num_basis[1:num_basis] <- rep(0, num_basis)
    a_raw[1:num_basis] ~ dmnorm(
        mean = zeroes_num_basis[1:num_basis],
        cov = I_num_basis[1:num_basis, 1:num_basis]
    )
    a0 ~ dnorm(mean = 0, sd = 1)
    tau ~ T(dnorm(mean = 0, sd = 1), 0, Inf)
    sigma ~ T(dnorm(mean = 0, sd = 1), 0, Inf)

    cov_matr[1:num_data, 1:num_data] <-
        sigma * diag(num_data)[1:num_data, 1:num_data]
    Y[1:num_data] ~ dmnorm(
        mean = Y_hat[1:num_data],
        cov = cov_matr[1:num_data, 1:num_data]
    )
})

splinesConsts <- list(
    num_data = num_data, num_basis = num_basis,
    X = X, B = B
)

splinesData <- list(Y = Y)

splinesInits <- list(
    a0 = 0.2, a_raw = rep(0, num_basis),
    tau = 1, sigma = 0.2
)

splinesModel <- nimbleModel(
    code = splinesCode, name = "splinesModel",
    constants = splinesConsts, data = splinesData,
    inits = splinesInits
)

Csplines <- compileNimble(splinesModel)

splinesConf <- configureMCMC(splinesModel, print = TRUE)

splinesConf$addMonitors(c("sigma", "tau", "a0", "a"))

splinesMCMC <- buildMCMC(splinesConf)

CsplinesMCMC <- compileNimble(splinesMCMC, project = splinesModel)

niter <- 500
set.seed(1)
samples <- runMCMC(CsplinesMCMC, niter = niter)

# trace plots
# par(mfrow = c(3, 1))
# plot(samples[, "sigma"], type = "l")
# plot(samples[, "tau"], type = "l")
# plot(samples[, "gamma_coef[1]"], type = "l")


# plotting the fit
par(mfrow = c(1, 1))
mean_a <- vector(length = num_basis)

names_vectors <- paste0("a[", 1:num_basis, "]")

for (i in 1:num_basis) {
    mean_a[i] <- mean(samples[, names_vectors[i]])
}

mean_a0 <- mean(samples[, "a0"])

y_fit <- mean_a0 * X + as.vector(a %*% B)

plot(X, Y)
lines(X, y_fit)
