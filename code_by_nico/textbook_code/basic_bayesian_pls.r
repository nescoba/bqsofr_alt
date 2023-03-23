library(nimble)

# parameters of the model
final_time <- 100
time <- 1:final_time
true_gamma <- sin(time * pi / 600)
sigma_2 <- .01
tau_2 <- .01

# creating the data
epsilon <- rnorm(final_time, sd = sigma_2)
true_y <- true_gamma + epsilon

# Creating the D, K and associated matrices
diff_matrix <- matrix(nrow = final_time - 2, ncol = final_time)
for (i in 1:(final_time - 2)) {
    diff_matrix[i, ] <- c(
        rep(0, i - 1), 1, -2, 1,
        rep(0, final_time - 2 - i)
    )
}
kernel_matrix <- t(diff_matrix) %*% diff_matrix
matrix_1 <- sigma_2 * diag(final_time)
matrix_2 <- kernel_matrix / tau_2
zeroes <- rep(0, final_time)


# implementing the code in NIMBLE
model_specification <- nimbleCode({
    y[1:100] ~ dmnorm(mean = gamma[1:100], cov = matrix_1[1:100, 1:100])
    gamma[1:100] ~ dmnorm(
        mean = zeroes[1:100],
        prec = I[1:100, 1:100]
    )
})

model_constants <- list(
    matrix_1 = matrix_1,
    zeroes = zeroes, matrix_2 = matrix_2, I = diag(100)
)
model_data <- list(y = true_y)
model_inits <- list(gamma = rep(0, final_time))

mcmc_out <- nimbleMCMC(
    code = model_specification, data = model_data, inits = model_inits,
    monitors = c("gamma"), constants = model_constants, summary = TRUE
)

mcmc_out$summary
