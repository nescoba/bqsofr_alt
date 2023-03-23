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

fit <- brm(bf(y ~ x), data = simulated_data)

stancode(fit)
