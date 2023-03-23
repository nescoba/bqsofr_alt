library("rstan")

schools_dat <- list(
    J = 8,
    y = c(28, 8, -3, 7, -1, 1, 18, 12),
    sigma = c(15, 10, 16, 11, 9, 11, 10, 18)
)

fit <- stan(file = "testing_stan/schools.stan", data = schools_dat)

install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
example(stan_model, package = "rstan", run.dontrun = TRUE)
install.packages("cli")
remove.packages("cli")
