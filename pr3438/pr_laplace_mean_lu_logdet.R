# pr_laplace_mean_lu_logdet.R
#
# Examples for a stan-dev/math PR (self-contained; the Stan programs are
# inline):
#
#   A. laplace_latent_solve_tol() returned K * grad log p(y | theta_hat)
#      instead of the final Newton iterate theta_hat.  Before exact
#      convergence the two differ by K * grad log p(theta | y), which is
#      large when the posterior is much narrower than the prior.
#   B. laplace_marginal_tol() with solver 3 returned NaN when the LU
#      factorization of B = I + K W has a negative pivot (also for a
#      Gaussian likelihood).
#   C. At a saddle point of log p(theta | y), the solver-3 log determinant
#      now reports that det(B) is not positive.
#
#   Rscript pr_laplace_mean_lu_logdet.R
#   CMDSTAN_PATH=<CmdStan tree> STAN_MATH=<Stan Math tree>/ \
#     Rscript pr_laplace_mean_lu_logdet.R

library(cmdstanr)
library(posterior)
if (nzchar(Sys.getenv("CMDSTAN_PATH"))) {
  set_cmdstan_path(Sys.getenv("CMDSTAN_PATH"))
}
cpp <- list()
if (nzchar(Sys.getenv("STAN_MATH"))) {
  cpp$MATH <- Sys.getenv("STAN_MATH")
}
compile <- function(code) {
  cmdstan_model(write_stan_file(code), cpp_options = cpp,
                force_recompile = length(cpp) > 0)
}
cat(sprintf("CmdStan %s at %s, Stan Math %s\n", cmdstan_version(),
            cmdstan_path(), if (is.null(cpp$MATH)) "of the CmdStan" else cpp$MATH))

# One fixed_param draw; returns the values of `vars` (NA if the draw failed)
# and the "Exception" lines of the output.
run_once <- function(mod, data, vars) {
  out <- capture.output(
    fit <- mod$sample(data = data, fixed_param = TRUE, chains = 1,
                      iter_sampling = 1, refresh = 0, show_messages = TRUE,
                      sig_figs = 12),
    type = "output")
  vals <- tryCatch(as.numeric(as_draws_matrix(fit$draws(vars))),
                   error = function(e) NA_real_)
  list(vals = vals,
       exception = unique(sub(".*Exception: (.*?) \\(in '.*", "\\1",
                              grep("Exception", out, value = TRUE),
                              perl = TRUE)))
}

# ── A. Returned mean ───────────────────────────────────────────────────────
tumor_code <- "
functions {
  real ll(vector theta, vector y, vector t, vector mu, real sigma) {
    real bsld = mu[1] * exp(theta[1]);
    real d = mu[2] * exp(theta[2]);
    real g = mu[3] * exp(theta[3]);
    real phi = inv_logit(logit(mu[4]) + theta[4]);
    vector[num_elements(y)] m
        = bsld * (phi * exp(-d * t) + (1 - phi) * exp(g * t));
    return normal_lpdf(y | m, sigma * m);
  }
  matrix k_fun(vector omega) {
    return diag_matrix(square(omega));
  }
}
data {
  int<lower=1> n;
  vector[n] y;
  vector[n] t;
  vector[4] mu;
  vector<lower=0>[4] omega;
  real<lower=0> sigma;
  vector[4] theta0;
  real<lower=0> tol;
}
generated quantities {
  vector[4] mode;
  {
    tuple(vector[4], matrix[4, 4]) lap = laplace_latent_solve_tol(
      ll, (y, t, mu, sigma), 4, k_fun, (omega,),
      (theta0, tol, 100, 3, 100, 0));
    mode = lap.1;
  }
}
"
y <- c(201.04, 78.30, 48.01, 44.39, 52.91, 50.49, 48.63, 67.29, 52.21,
       72.43, 102.43, 105.06)
t <- seq(0, 693, by = 63)
mu <- c(117.4, 9.121e-3, 1.913e-3, 0.3781)
omega <- c(0.4925, 0.9217, 1.120, 1.299)
sigma <- 0.09826
log_f <- function(th) {
  p <- plogis(qlogis(mu[4]) + th[4])
  m <- mu[1] * exp(th[1]) * (p * exp(-mu[2] * exp(th[2]) * t) +
                               (1 - p) * exp(mu[3] * exp(th[3]) * t))
  sum(dnorm(y, m, sigma * m, log = TRUE)) + sum(dnorm(th, 0, omega, log = TRUE))
}
mode_ref <- optim(rep(0, 4), function(z) -log_f(z), method = "BFGS",
                  control = list(reltol = 1e-15, maxit = 5000))$par
for (it in 1:5) {
  mode_ref <- mode_ref - solve(numDeriv::hessian(log_f, mode_ref),
                               numDeriv::grad(log_f, mode_ref))
}
post_sd <- sqrt(diag(solve(-numDeriv::hessian(log_f, mode_ref))))
cat("\nA. Returned mean.  Mode (optim + Newton):",
    sprintf("%.4f", mode_ref), "\n")
cat("   posterior sd:", sprintf("%.3f", post_sd), " prior sd:", omega, "\n")
mod_a <- compile(tumor_code)
for (tol in c(1e-6, 1.49012e-8, 1e-10)) {
  r <- run_once(mod_a, list(n = length(y), y = y, t = t, mu = mu,
                            omega = omega, sigma = sigma,
                            theta0 = c(0.3636, -0.3385, 1.9794, -1.9322),
                            tol = tol), "mode")
  cat(sprintf("   tol %-9g returned %s: %.2g posterior sd from the mode, log density %.3g lower\n",
              tol, paste(sprintf("%.4f", r$vals), collapse = " "),
              max(abs(r$vals - mode_ref) / post_sd),
              log_f(mode_ref) - log_f(r$vals)))
}

# ── B. Solver 3 with a negative LU pivot ───────────────────────────────────
normal_code <- "
functions {
  real ll(vector theta, vector y, real s) {
    return normal_lpdf(y | theta, s);
  }
  matrix k_fun(matrix K) {
    return K;
  }
}
data {
  vector[2] y;
  real<lower=0> s;
  matrix[2, 2] K;
  int<lower=1, upper=3> solver;
}
generated quantities {
  real lml = laplace_marginal_tol(ll, (y, s), 2, k_fun, (K,),
                                  (rep_vector(0.0, 2), 1e-10, 100, solver,
                                   100, 0));
  real lml_exact = multi_normal_lpdf(y | rep_vector(0.0, 2),
                                     K + diag_matrix(rep_vector(square(s), 2)));
}
"
K_b <- matrix(c(1, 1.5, 1.5, 4), 2)
s_b <- sqrt(0.1)
B <- diag(2) + K_b %*% diag(1 / s_b^2, 2)
cat(sprintf("\nB. Negative LU pivot.  B = I + K W: det %.4g, LU diagonal (partial pivoting) %s\n",
            det(B), paste(sprintf("%.4g", diag(as.matrix(Matrix::expand(Matrix::lu(B))$U))),
                          collapse = " ")))
mod_b <- compile(normal_code)
for (sv in 1:3) {
  r <- run_once(mod_b, list(y = c(0.3, -0.2), s = s_b, K = K_b, solver = sv),
                c("lml", "lml_exact"))
  cat(sprintf("   solver %d: laplace_marginal_tol %.6f, exact %.6f\n", sv,
              r$vals[1], r$vals[2]))
}

# ── C. Saddle point ────────────────────────────────────────────────────────
saddle_code <- "
functions {
  real ll(vector theta, vector y, real s, real z, real sz) {
    return cauchy_lpdf(y | theta[1], s) + normal_lpdf(z | theta[2], sz);
  }
  matrix k_fun(matrix K) {
    return K;
  }
}
data {
  vector[2] y;
  real<lower=0> s;
  real z;
  real<lower=0> sz;
  matrix[2, 2] K;
  vector[2] theta0;
  int<lower=0, upper=1> marginal;
}
generated quantities {
  real lml = not_a_number();
  vector[2] mode = rep_vector(not_a_number(), 2);
  if (marginal) {
    lml = laplace_marginal_tol(ll, (y, s, z, sz), 2, k_fun, (K,),
                               (theta0, 1e-10, 100, 3, 100, 0));
  } else {
    tuple(vector[2], matrix[2, 2]) lap = laplace_latent_solve_tol(
      ll, (y, s, z, sz), 2, k_fun, (K,), (theta0, 1e-10, 100, 3, 100, 0));
    mode = lap.1;
  }
}
"
mod_c <- compile(saddle_code)
data_c <- list(y = c(-3, 3), s = 0.3, z = 0, sz = 0.1,
               K = matrix(c(4, 1, 1, 1), 2), theta0 = c(0.01, 3))
cat("\nC. Saddle point at (0, 0); start (0.01, 3), solver 3\n")
r <- run_once(mod_c, c(data_c, list(marginal = 1L)), "lml")
cat(sprintf("   laplace_marginal_tol: %s\n",
            if (length(r$exception)) paste("exception:", r$exception)
            else sprintf("returned %.6g", r$vals)))
r <- run_once(mod_c, c(data_c, list(marginal = 0L)), "mode")
cat(sprintf("   laplace_latent_solve_tol: %s\n",
            if (length(r$exception)) paste("exception:", r$exception)
            else sprintf("returned %s", paste(sprintf("%.4g", r$vals), collapse = " "))))
