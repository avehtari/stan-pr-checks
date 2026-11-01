# stan-dev/stan#3088: offset/multiplier "trapped at low sigma".
# Step 1: partial inits (R1).  Step 2: default inits, built-in vs manual
# non-centering (R2).  Step 3: gradients of the generated code.
# Run on Triton: Rscript run_3088.R <work_dir> [step ...]
# Needs CMDSTAN set to the CmdStan develop tree.

library(cmdstanr)
library(posterior)

args <- commandArgs(trailingOnly = TRUE)
WORK <- if (length(args) >= 1) args[1] else getwd()
STEPS <- if (length(args) >= 2) args[-1] else c("1", "2", "3")
N_SEEDS <- as.integer(Sys.getenv("N_SEEDS", "20"))
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
TMP <- Sys.getenv("TMPDIR", tempdir())
OUT <- file.path(WORK, "out")
BIN <- file.path(WORK, paste0("bin_", Sys.getenv("SLURM_JOB_ID", "local")))
dir.create(OUT, showWarnings = FALSE)
dir.create(BIN, showWarnings = FALSE)
set_cmdstan_path(Sys.getenv("CMDSTAN"))
cat(sprintf("cmdstan %s, cmdstanr %s, cores %d, seeds %d\n", cmdstan_version(),
            as.character(packageVersion("cmdstanr")), N_CORES, N_SEEDS))

# ── Models ──
S2_MODELS <- strsplit(Sys.getenv("S2_MODELS",
  "j_ncp,j_cov,j_om,j_m,k_om,k_ncp,d_om,d_ncp,p_m,p_ncp"), "[,:]")[[1]]
METHODS_MODELS <- c("j_ncp", "j_cov", "j_om", "j_m")
MODEL_NAMES <- unique(c(if ("1" %in% STEPS) c("r1_m", "r1_manual", "r1_om", "r1_lower"),
                        if ("2" %in% STEPS) S2_MODELS,
                        if ("3" %in% STEPS) METHODS_MODELS))
compile_one <- function(name) {
  src <- file.path(WORK, paste0(name, ".stan"))
  file.copy(src, file.path(BIN, basename(src)), overwrite = TRUE)
  cmdstan_model(file.path(BIN, basename(src)))
}
mods <- parallel::mclapply(MODEL_NAMES, function(n) {
  tryCatch(compile_one(n), error = function(e) conditionMessage(e))
}, mc.cores = min(N_CORES, 8))
names(mods) <- MODEL_NAMES
bad <- MODEL_NAMES[!vapply(mods, inherits, logical(1), "CmdStanModel")]
if (length(bad) > 0) stop(sprintf("compile failed: %s", paste(bad, collapse = ", ")))
# reload in the master process (the mclapply objects hold child paths only)
mods <- lapply(MODEL_NAMES, function(n) {
  with_methods <- n %in% METHODS_MODELS
  cmdstan_model(file.path(BIN, paste0(n, ".stan")), compile_model_methods = with_methods,
                force_recompile = with_methods)
})
names(mods) <- MODEL_NAMES

# ── Data ──
make_j_data <- function() {
  set.seed(3)
  n_grp <- 1000
  n_per_grp <- 5
  re <- rnorm(n_grp, 0, 0.5)
  grp_id <- rep(1:n_grp, n_per_grp)
  list(N = n_grp * n_per_grp, n_grp = n_grp, y = rnorm(n_grp * n_per_grp, re[grp_id], 1),
       grp_id = grp_id)
}
make_d_data <- function() {
  set.seed(1)
  I <- 10000
  p <- plogis(rnorm(I, mean = 0.69, sd = 0.5))
  N <- rpois(I, rgamma(I, 20, 1))
  list(I = I, N = N, y = rbinom(I, N, p))
}
P_N <- as.integer(Sys.getenv("P_N", "50000"))
make_p_data <- function() {
  set.seed(10)
  list(n = P_N, y = rnorm(P_N))
}
DATA <- list(j = make_j_data(), d = make_d_data(), p = make_p_data(), k = list(),
             r1 = list())
data_for <- function(name) DATA[[sub("_.*$", "", name)]]

# ── Step 1: partial inits ──
# Calls the CmdStan binary directly: cmdstanr stops when the autodiff and
# finite-difference gradients disagree, and it hides the init messages.
diag_one <- function(name, seed, init) {
  init_file <- tempfile(fileext = ".json")
  jsonlite::write_json(init, init_file, auto_unbox = TRUE, digits = NA)
  out <- suppressWarnings(system2(mods[[name]]$exe_file(),
    c("diagnose", "test=gradient", "error=1e300", sprintf("init=%s", init_file),
      "random", sprintf("seed=%d", seed), "output", sprintf("file=%s", tempfile())),
    stdout = TRUE, stderr = TRUE))
  n_reject <- sum(grepl("Rejecting initial value", out))
  lp_line <- grep("Log probability=", out, value = TRUE)
  idx <- grep("param idx", out)
  if (length(idx) == 0 || length(lp_line) == 0) {
    return(list(model = name, seed = seed, ok = FALSE, n_reject = n_reject,
                msg = paste(tail(out, 3), collapse = " | ")))
  }
  rows <- out[(idx + 1):length(out)]
  rows <- rows[grepl("^\\s*[0-9]+\\s", rows)]
  tab <- do.call(rbind, lapply(strsplit(trimws(rows), "\\s+"), as.numeric))
  list(model = name, seed = seed, ok = TRUE, n_reject = n_reject,
       lp = as.numeric(sub(".*Log probability=", "", lp_line[1])),
       value = tab[, 2], grad = tab[, 3])
}

if ("1" %in% STEPS) {
  cat("\n== Step 1: partial inits, diagnose at the initial point ==\n")
  cases <- list(
    list(name = "r1_m", init = list(sigma = 1e-20)),
    list(name = "r1_manual", init = list(sigma = 1e-20)),
    list(name = "r1_om", init = list(sigma = 1e-20)),
    list(name = "r1_om", init = list(mu = 5)),
    list(name = "r1_lower", init = list(a = 10)),
    list(name = "r1_lower", init = list(a = -10)))
  s1 <- list()
  for (k in seq_along(cases)) {
    cs <- cases[[k]]
    r <- lapply(seq_len(N_SEEDS), function(s) diag_one(cs$name, s, cs$init))
    lbl <- sprintf("%s[%s]", cs$name, paste(names(cs$init), unlist(cs$init), sep = "=", collapse = ","))
    s1[[lbl]] <- r
    ok <- vapply(r, function(x) isTRUE(x$ok), logical(1))
    n_rej <- vapply(r, function(x) x$n_reject, numeric(1))
    cat(sprintf("%-28s init ok %2d/%d; rejected inits per seed: min %d max %d\n", lbl, sum(ok),
                length(ok), min(n_rej), max(n_rej)))
    if (any(ok)) {
      v <- t(vapply(r[ok], function(x) x$value, numeric(length(r[ok][[1]]$value))))
      lp <- vapply(r[ok], function(x) as.numeric(x$lp)[1], numeric(1))
      for (j in seq_len(ncol(v))) {
        cat(sprintf("   unconstrained[%d]: |min| %9.3g  |max| %9.3g\n", j, min(abs(v[, j])),
                    max(abs(v[, j]))))
      }
      cat(sprintf("   lp at init: min %9.3g  max %9.3g\n", min(lp), max(lp)))
    } else {
      cat(sprintf("   first error: %s\n", substr(r[[1]]$msg, 1, 200)))
    }
  }
  saveRDS(s1, file.path(OUT, "step1_diagnose.rds"))

  cat("\n== Step 1b: sampling with init sigma = 1e-20 ==\n")
  samp1 <- function(name, seed) {
    fit <- mods[[name]]$sample(seed = seed, init = list(list(sigma = 1e-20)), chains = 1,
                               iter_warmup = 1000, iter_sampling = 1000, refresh = 0,
                               show_messages = FALSE, show_exceptions = FALSE,
                               output_dir = TMP)
    d <- fit$draws(c("lp__", "sigma"), format = "draws_matrix")
    sd_ <- fit$sampler_diagnostics(format = "draws_matrix")
    out <- data.frame(model = name, seed = seed, lp_mean = mean(d[, "lp__"]),
                      sd_sigma = sd(d[, "sigma"]), med_sigma = median(d[, "sigma"]),
                      stepsize = sd_[1, "stepsize__"], accept = mean(sd_[, "accept_stat__"]),
                      n_div = sum(sd_[, "divergent__"]), stringsAsFactors = FALSE)
    unlink(fit$output_files())
    out
  }
  grid <- expand.grid(name = c("r1_m", "r1_manual"), seed = seq_len(N_SEEDS),
                      stringsAsFactors = FALSE)
  r <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
    tryCatch(samp1(grid$name[i], grid$seed[i]), error = function(e) {
      data.frame(model = grid$name[i], seed = grid$seed[i], lp_mean = NA, sd_sigma = NA,
                 med_sigma = NA, stepsize = NA, accept = NA, n_div = NA)
    })
  }, mc.cores = N_CORES)
  s1b <- do.call(rbind, r)
  saveRDS(s1b, file.path(OUT, "step1_sample.rds"))
  for (m in unique(s1b$model)) {
    x <- s1b[s1b$model == m, ]
    cat(sprintf("%-10s stuck (sd_sigma == 0) %2d/%d; lp_mean range %9.3g %9.3g; stepsize med %9.3g\n",
                m, sum(x$sd_sigma == 0, na.rm = TRUE), nrow(x), min(x$lp_mean, na.rm = TRUE),
                max(x$lp_mean, na.rm = TRUE), median(x$stepsize, na.rm = TRUE)))
  }
}

# ── Step 2: default inits, warmup behavior ──
SCALE <- c(j = "sigma", d = "sigma", p = "sigma", k = "s")
run_one <- function(name, seed) {
  t0 <- Sys.time()
  big <- sub("_.*$", "", name) %in% c("d", "p")
  fit <- mods[[name]]$sample(data = data_for(name), seed = seed, chains = 1,
                             iter_warmup = 1000, iter_sampling = if (big) 200 else 1000,
                             save_warmup = TRUE, refresh = 0, show_messages = FALSE,
                             show_exceptions = FALSE, output_dir = TMP)
  sc <- SCALE[[sub("_.*$", "", name)]]
  csv <- read_cmdstan_csv(fit$output_files(), variables = c("lp__", sc),
                          sampler_diagnostics = c("stepsize__", "accept_stat__",
                                                  "divergent__", "treedepth__"))
  unlink(fit$output_files())
  w <- as_draws_matrix(csv$warmup_draws)
  s <- as_draws_matrix(csv$post_warmup_draws)
  wd <- as_draws_matrix(csv$warmup_sampler_diagnostics)
  sdg <- as_draws_matrix(csv$post_warmup_sampler_diagnostics)
  sw <- as.numeric(w[, sc])
  ss <- as.numeric(s[, sc])
  list(row = data.frame(model = name, seed = seed, min_scale_w = min(sw),
                        it_min = which.min(sw), first_above = {
                          k <- which(seq_along(sw) > which.min(sw) & sw > 0.1 * median(ss))
                          if (length(k) > 0) k[1] else NA_integer_
                        },
                        med_scale_s = median(ss), sd_scale_s = sd(ss),
                        n_div_w = sum(wd[, "divergent__"]), n_div_s = sum(sdg[, "divergent__"]),
                        stepsize = as.numeric(sdg[1, "stepsize__"]),
                        accept_s = mean(sdg[, "accept_stat__"]),
                        treedepth_s = mean(sdg[, "treedepth__"]),
                        secs = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                        stringsAsFactors = FALSE),
       trace_w = sw)
}

if ("2" %in% STEPS) {
  cat("\n== Step 2: default inits ==\n")
  grid <- expand.grid(seed = seq_len(N_SEEDS), name = S2_MODELS, stringsAsFactors = FALSE)
  S2_DIR <- file.path(OUT, "step2")
  dir.create(S2_DIR, showWarnings = FALSE)
  r <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
    tag <- if (startsWith(grid$name[i], "p_")) sprintf("%s_n%d", grid$name[i], P_N) else grid$name[i]
    f <- file.path(S2_DIR, sprintf("%s_seed%02d.rds", tag, grid$seed[i]))
    if (file.exists(f)) return(readRDS(f))
    x <- tryCatch(run_one(grid$name[i], grid$seed[i]),
                  error = function(e) list(row = data.frame(model = grid$name[i],
                                                            seed = grid$seed[i],
                                                            error = conditionMessage(e),
                                                            stringsAsFactors = FALSE)))
    if (is.null(x$row$error)) saveRDS(x, f)
    x
  }, mc.cores = N_CORES, mc.preschedule = FALSE)
  rows <- lapply(r, function(x) x$row)
  ok <- vapply(rows, function(x) is.null(x$error), logical(1))
  if (any(!ok)) {
    cat(sprintf("%d runs failed; first: %s\n", sum(!ok), rows[!ok][[1]]$error))
  }
  s2 <- do.call(rbind, rows[ok])
  cat(sprintf("%-6s %4s %9s %9s %9s %6s %6s %9s %5s\n", "model", "runs", "minscW",
              "n<1e-12", "stuck", "divW", "divS", "stepsize", "tdS"))
  for (m in S2_MODELS) {
    x <- s2[s2$model == m, ]
    if (nrow(x) == 0) next
    cat(sprintf("%-6s %4d %9.2g %9d %9d %6d %6d %9.3g %5.2f\n", m, nrow(x), min(x$min_scale_w),
                sum(x$min_scale_w < 1e-12), sum(x$med_scale_s < 1e-6),
                sum(x$n_div_w), sum(x$n_div_s), median(x$stepsize), mean(x$treedepth_s)))
  }
}

# ── Step 3: gradients of the generated code at small sigma ──
if ("3" %in% STEPS) {
  cat("\n== Step 3: grad_log_prob of the J models against j_ncp ==\n")
  jd <- DATA$j
  fits <- lapply(METHODS_MODELS, function(n) {
    f <- mods[[n]]$sample(data = jd, seed = 1, chains = 1, iter_warmup = 20, iter_sampling = 5,
                          refresh = 0, show_messages = FALSE, output_dir = TMP)
    f$init_model_methods()
    f
  })
  names(fits) <- METHODS_MODELS
  set.seed(20261031)
  raw <- runif(jd$n_grp, -2, 2)
  rows <- list()
  for (mu in c(0, 1, 1000)) {
    for (e in seq(-2, -30, by = -2)) {
      u <- c(mu, e * log(10), 0, raw)
      g <- lapply(fits, function(f) f$grad_log_prob(u, jacobian = TRUE))
      ref <- g$j_ncp
      for (n in setdiff(METHODS_MODELS, "j_ncp")) {
        d <- abs(g[[n]] - ref)
        rows[[length(rows) + 1]] <- data.frame(model = n, mu = mu, sigma = 10^e,
                                               err_mu = d[1], err_logsigma = d[2],
                                               err_raw_max = max(d[-(1:3)]),
                                               ref_mu = ref[1], ref_logsigma = ref[2],
                                               stringsAsFactors = FALSE)
      }
    }
  }
  s3 <- do.call(rbind, rows)
  saveRDS(s3, file.path(OUT, "step3_grad.rds"))
  cat(sprintf("%-6s %5s %8s %10s %10s %10s %10s\n", "model", "mu", "sigma", "ref_mu",
              "err_mu", "err_lsig", "err_raw"))
  for (i in seq_len(nrow(s3))) {
    x <- s3[i, ]
    cat(sprintf("%-6s %5g %8.0e %10.3g %10.2g %10.2g %10.2g\n", x$model, x$mu, x$sigma,
                x$ref_mu, x$err_mu, x$err_logsigma, x$err_raw_max))
  }
}
cat("\nEXIT=0\n")
