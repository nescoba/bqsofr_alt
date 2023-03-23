library(splines)
library(nimble)

# parameter
D <- 3

# useful function
inv.logit <- function(x) 1 / (1 + exp(-x))

# grid of x values
x <- seq(0, 1, .01)

# creating spline basis
X <- bs(x = x, degree = D, knots = seq(0, 1, .1))

# computed parameters
K <- ncol(X)
num_data <- length(x)

# creating synthetic data
true_gamma <- rnorm(n = K, sd = 1)
true_gamma_col <- matrix(true_gamma, ncol = 1)
true_trend <- X %*% true_gamma_col
true_trend_vec <- as.vector(true_trend)
y <- rbinom(n = num_data, size = 1, prob = inv.logit(true_trend_vec))

# implementing inv.logit as a nimble function
inv.logit.nimble <- nimbleFunction(
    run = function(x = double(1)) {
        ans <- 1 / (1 + exp(-x))
        return(ans)
        returnType(double(1))
    }
)

# creating the model
splinesCode <- nimbleCode({
    # smoothing prior for coefficients
    gamma_est[1] ~ dflat()
    for (i in 2:K) {
        gamma_est[i] ~ dnorm(mean = gamma_est[i - 1], var = tau_2)
    }

    # defining the trend
    trend[1:num_data, 1] <- X[1:num_data, 1:K] %*% asCol(gamma_est[1:K])

    # prior for tau_2
    tau_2 ~ dinvgamma(shape = a, scale = b)

    # computing the estimated probabilites
    p_est[1:num_data] <- inv.logit.nimble(trend[1:num_data, 1])

    # sampling model
    for (i in 1:num_data) {
        y[i] ~ dbern(p_est[i])
    }
})

splinesConsts <- list(num_data = num_data, X = X, a = 0.1, b = 0.1, K = K)

splinesData <- list(y = y)

splinesInits <- list(gamma_est = rep(0, K), tau_2 = 1)

splinesModel <- nimbleModel(
    code = splinesCode, name = "splinesModel",
    constants = splinesConsts, data = splinesData,
    inits = splinesInits
)

Csplines <- compileNimble(splinesModel)

splinesConf <- configureMCMC(splinesModel, print = TRUE)

splinesConf$addMonitors(c("tau_2", "gamma_est"))

splinesMCMC <- buildMCMC(splinesConf)

CsplinesMCMC <- compileNimble(splinesMCMC, project = splinesModel)

niter <- 1000
set.seed(1)
samples <- runMCMC(CsplinesMCMC, niter = niter)

mean_gamma_est <- vector(length = K)

names_vectors <- paste0("gamma_est[", 1:K, "]")

for (i in 1:K) {
    mean_gamma_est[i] <- mean(samples[, names_vectors[i]])
}

mean_trend <- X %*% matrix(mean_gamma_est, ncol = 1)

plot(x, true_trend_vec)
lines(x, mean_trend)
