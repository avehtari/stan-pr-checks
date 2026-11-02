# Prototype check of the stanc3 offset/multiplier std_normal rewrite
# (stanc3 branch offset-multiplier-std-normal, flag
# -foffset-multiplier-std-normal).
# Variants: "base" = the prototype stanc without the flag, "om" = with it.
# make runs stanc as $(WINE) bin/stanc ...; the driver sets WINE to a
# wrapper that runs the prototype stanc instead.
# Run with tools/cmdstan_run.sh; files: stanc_om_static, models, this script.
# Rscript run_om_std_normal.R <work_dir> [step ...]   steps: 1 2 3 4
# Step 4 also needs om_nodata.stan and om_probe_main.hpp.

library(cmdstanr)
library(posterior)

args <- commandArgs(trailingOnly = TRUE)
WORK <- if (length(args) >= 1) normalizePath(args[1]) else getwd()
STEPS <- if (length(args) >= 2) args[-1] else c("1", "2", "3", "4")
N_SEEDS <- as.integer(Sys.getenv("N_SEEDS", "20"))
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))
TMP <- Sys.getenv("TMPDIR", tempdir())
OUT <- file.path(WORK, "out")
dir.create(OUT, showWarnings = FALSE)
set_cmdstan_path(Sys.getenv("CMDSTAN"))
cat(sprintf("cmdstan %s, cmdstanr %s, cores %d, seeds %d\n", cmdstan_version(),
            as.character(packageVersion("cmdstanr")), N_CORES, N_SEEDS))

# ── stanc wrappers ──
STANC <- file.path(WORK, "stanc_om_static")
Sys.chmod(STANC, "0755")
FLAGS <- c(base = "", om = "-foffset-multiplier-std-normal")
WRAP <- vapply(names(FLAGS), function(v) {
  f <- file.path(WORK, sprintf("stanc_%s.sh", v))
  writeLines(c("#!/bin/bash", "# make calls $(WINE) bin/stanc ARGS; drop bin/stanc", "shift",
               sprintf("echo \"%s $*\" >> %s", v, file.path(WORK, "stanc_calls.log")),
               sprintf("exec %s %s \"$@\"", STANC, FLAGS[[v]])), f)
  Sys.chmod(f, "0755")
  f
}, character(1))
cat(sprintf("stanc: %s\n", system2(STANC, "--version", stdout = TRUE)[1]))

# Does the rewrite fire?  Count std_normal_lpdf calls in the generated C++.
MODELS <- c("j_om", "k_om", "d_om", "j_ncp", "om_cases")
for (n in MODELS) {
  cnt <- vapply(names(FLAGS), function(v) {
    hpp <- tempfile(fileext = ".hpp")
    a <- c(if (nzchar(FLAGS[[v]])) FLAGS[[v]], sprintf("--o=%s", hpp),
           file.path(WORK, paste0(n, ".stan")))
    system2(STANC, a)
    sum(grepl("std_normal_lpdf", readLines(hpp)))
  }, numeric(1))
  cat(sprintf("%-9s std_normal_lpdf lines: base %d, om %d\n", n, cnt[["base"]], cnt[["om"]]))
}

# ── Compile ──
# cmdstanr builds model methods with bin/stanc directly, not through make,
# so steps 1 and 2 use the log_prob method of the executables instead.
grid <- rbind(expand.grid(name = c("j_om", "k_om", "d_om", "om_cases"),
                          variant = names(FLAGS), stringsAsFactors = FALSE),
              data.frame(name = "j_ncp", variant = "base", stringsAsFactors = FALSE))
grid$key <- sprintf("%s__%s", grid$name, grid$variant)
stan_path <- function(i) file.path(WORK, paste0("bin_", grid$variant[i]), paste0(grid$name[i], ".stan"))
compile_one <- function(i) {
  dir.create(dirname(stan_path(i)), showWarnings = FALSE)
  file.copy(file.path(WORK, paste0(grid$name[i], ".stan")), stan_path(i), overwrite = TRUE)
  cmdstan_model(stan_path(i), cpp_options = list(WINE = WRAP[[grid$variant[i]]]),
                force_recompile = TRUE)
}
r <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
  tryCatch({ compile_one(i); "ok" }, error = function(e) conditionMessage(e))
}, mc.cores = min(N_CORES, nrow(grid)), mc.preschedule = FALSE)
bad <- grid$key[unlist(r) != "ok"]
if (length(bad) > 0) {
  cat(sprintf("compile failed: %s\n%s\n", paste(bad, collapse = ", "), r[[which(unlist(r) != "ok")[1]]]))
  quit(status = 1)
}
mods <- lapply(seq_len(nrow(grid)), function(i) {
  cmdstan_model(stan_path(i), cpp_options = list(WINE = WRAP[[grid$variant[i]]]))
})
names(mods) <- grid$key
calls <- readLines(file.path(WORK, "stanc_calls.log"))
cat(sprintf("wrapper calls: base %d, om %d\n", sum(startsWith(calls, "base ")),
            sum(startsWith(calls, "om "))))

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
DATA <- list(j = make_j_data(), d = make_d_data(), k = list(),
             om = list(N = 3, s_data = 0.7, s_vec_data = c(0.5, 1, 2)))
data_for <- function(name) DATA[[sub("_.*$", "", name)]]

DATA_FILE <- vapply(c("j", "om"), function(g) {
  f <- file.path(WORK, sprintf("data_%s.json", g))
  write_stan_json(DATA[[g]], f)
  f
}, character(1))
U_FILE <- file.path(WORK, "u.json")
# log density and gradient from the executable (CmdStan log_prob method)
lp_grad <- function(key, jacobian) {
  of <- tempfile(fileext = ".csv")
  out <- system2(mods[[key]]$exe_file(),
                 c("log_prob", sprintf("unconstrained_params=%s", U_FILE),
                   sprintf("jacobian=%d", as.integer(jacobian)), "data",
                   sprintf("file=%s", DATA_FILE[[sub("_.*$", "", key)]]), "output",
                   sprintf("file=%s", of), "sig_figs=18"), stdout = TRUE, stderr = TRUE)
  if (!file.exists(of)) stop(sprintf("%s log_prob failed: %s", key, paste(tail(out, 3), collapse = " | ")))
  v <- as.numeric(read.csv(of, comment.char = "#")[1, ])
  unlink(of)
  list(lp = v[1], grad = v[-1])
}
write_u <- function(u) jsonlite::write_json(list(params_r = u), U_FILE, digits = NA)

# ── Step 1: same log_prob and gradient on om_cases ──
if ("1" %in% STEPS) {
  cat("\n== Step 1: om_cases, om against base at random unconstrained points ==\n")
  # mu, sigma, sigma_vec[3], a, b, c, e, f, g, h, k [3 each], d, m[2, 2]
  n_par <- 1 + 1 + 3 + 8 * 3 + 1 + 4
  set.seed(1101)
  rows <- list()
  for (k in 1:20) {
    u <- rnorm(n_par, 0, 1.5)
    write_u(u)
    for (jac in c(TRUE, FALSE)) {
      b <- lp_grad("om_cases__base", jac)
      o <- lp_grad("om_cases__om", jac)
      lb <- b$lp
      lo <- o$lp
      gb <- b$grad
      go <- o$grad
      rows[[length(rows) + 1]] <- data.frame(point = k, jacobian = jac, lp = lb,
                                             lp_diff = lo - lb,
                                             grad_rel = max(abs(go - gb) / pmax(1, abs(gb))))
    }
  }
  s1 <- do.call(rbind, rows)
  saveRDS(s1, file.path(OUT, "step1_om_cases.rds"))
  for (jac in c(TRUE, FALSE)) {
    x <- s1[s1$jacobian == jac, ]
    cat(sprintf("jacobian=%-5s points %d: max |lp diff| %.2g (|lp| up to %.3g), max grad rel diff %.2g\n",
                jac, nrow(x), max(abs(x$lp_diff)), max(abs(x$lp)), max(x$grad_rel)))
  }
}

# ── Step 2: gradients of j_om against j_ncp at small sigma ──
if ("2" %in% STEPS) {
  cat("\n== Step 2: grad_log_prob of j_om (base, om) against j_ncp ==\n")
  keys <- c("j_ncp__base", "j_om__base", "j_om__om")
  jd <- DATA$j
  set.seed(20261031)
  raw <- runif(jd$n_grp, -2, 2)
  rows <- list()
  for (mu in c(0, 1, 1000)) {
    for (e in seq(-2, -30, by = -4)) {
      u <- c(mu, e * log(10), 0, raw)
      write_u(u)
      g <- lapply(keys, function(k) lp_grad(k, TRUE)$grad)
      names(g) <- keys
      ref <- g$j_ncp__base
      for (n in keys[-1]) {
        d <- abs(g[[n]] - ref)
        rows[[length(rows) + 1]] <- data.frame(model = n, mu = mu, sigma = 10^e,
                                               err_mu = d[1], err_logsigma = d[2],
                                               err_raw_max = max(d[-(1:3)]),
                                               ref_mu = ref[1], stringsAsFactors = FALSE)
      }
    }
  }
  s2 <- do.call(rbind, rows)
  saveRDS(s2, file.path(OUT, "step2_grad.rds"))
  cat(sprintf("%-11s %5s %8s %10s %10s %10s %10s\n", "model", "mu", "sigma", "ref_mu",
              "err_mu", "err_lsig", "err_raw"))
  for (i in seq_len(nrow(s2))) {
    x <- s2[i, ]
    cat(sprintf("%-11s %5g %8.0e %10.3g %10.2g %10.2g %10.2g\n", x$model, x$mu, x$sigma,
                x$ref_mu, x$err_mu, x$err_logsigma, x$err_raw_max))
  }
}

# ── Step 3: default inits, 20 seeds per model and variant ──
SCALE <- c(j = "sigma", d = "sigma", k = "s")
run_one <- function(key, seed) {
  t0 <- Sys.time()
  name <- sub("__.*$", "", key)
  big <- startsWith(name, "d_")
  fit <- mods[[key]]$sample(data = data_for(name), seed = seed, chains = 1,
                            iter_warmup = 1000, iter_sampling = if (big) 200 else 1000,
                            save_warmup = TRUE, refresh = 0, show_messages = FALSE,
                            show_exceptions = FALSE, output_dir = TMP)
  sc <- SCALE[[sub("_.*$", "", name)]]
  csv <- read_cmdstan_csv(fit$output_files(), variables = c("lp__", sc),
                          sampler_diagnostics = c("stepsize__", "accept_stat__",
                                                  "divergent__", "treedepth__"))
  unlink(fit$output_files())
  sw <- as.numeric(as_draws_matrix(csv$warmup_draws)[, sc])
  ss <- as.numeric(as_draws_matrix(csv$post_warmup_draws)[, sc])
  wd <- as_draws_matrix(csv$warmup_sampler_diagnostics)
  sdg <- as_draws_matrix(csv$post_warmup_sampler_diagnostics)
  data.frame(key = key, seed = seed, min_scale_w = min(sw), med_scale_s = median(ss),
             n_div_w = sum(wd[, "divergent__"]), n_div_s = sum(sdg[, "divergent__"]),
             stepsize = as.numeric(sdg[1, "stepsize__"]),
             treedepth_s = mean(sdg[, "treedepth__"]),
             secs = as.numeric(difftime(Sys.time(), t0, units = "secs")),
             stringsAsFactors = FALSE)
}
if ("3" %in% STEPS) {
  cat("\n== Step 3: default inits (stuck: median scale after warmup < 1e-6) ==\n")
  keys <- c("j_om__base", "j_om__om", "k_om__base", "k_om__om", "d_om__base", "d_om__om")
  sg <- expand.grid(seed = seq_len(N_SEEDS), key = keys, stringsAsFactors = FALSE)
  # the slow D runs first
  sg <- sg[order(!startsWith(sg$key, "d_")), ]
  S3_DIR <- file.path(OUT, "step3")
  dir.create(S3_DIR, showWarnings = FALSE)
  r <- parallel::mclapply(seq_len(nrow(sg)), function(i) {
    f <- file.path(S3_DIR, sprintf("%s_seed%02d.rds", sg$key[i], sg$seed[i]))
    if (file.exists(f)) return(readRDS(f))
    x <- tryCatch(run_one(sg$key[i], sg$seed[i]), error = function(e) {
      data.frame(key = sg$key[i], seed = sg$seed[i], error = conditionMessage(e),
                 stringsAsFactors = FALSE)
    })
    if (is.null(x$error)) saveRDS(x, f)
    x
  }, mc.cores = N_CORES, mc.preschedule = FALSE)
  ok <- vapply(r, function(x) is.null(x$error), logical(1))
  if (any(!ok)) cat(sprintf("%d runs failed; first: %s\n", sum(!ok), r[!ok][[1]]$error))
  s3 <- do.call(rbind, r[ok])
  cat(sprintf("%-11s %4s %9s %6s %6s %6s %9s %5s %7s\n", "model", "runs", "minscW", "stuck",
              "divW", "divS", "stepsize", "tdS", "secs"))
  for (k in keys) {
    x <- s3[s3$key == k, ]
    if (nrow(x) == 0) next
    cat(sprintf("%-11s %4d %9.2g %6d %6d %6d %9.3g %5.2f %7.0f\n", k, nrow(x),
                min(x$min_scale_w), sum(x$med_scale_s < 1e-6), sum(x$n_div_w), sum(x$n_div_s),
                median(x$stepsize), mean(x$treedepth_s), median(x$secs)))
  }
}

# ── Step 4: om_nodata through a C++ probe, all (propto, jacobian) pairs ──
# The probe main replaces CmdStan's main.o.  With propto, a constant lp
# shift is expected: the rewrite keeps -L for a data multiplier (g, g2, h),
# where normal_lupdf drops it.
if ("4" %in% STEPS) {
  cat("\n== Step 4: om_nodata, C++ probe, variants against base (--O0) ==\n")
  PFLAGS <- c(base = "", om = "-foffset-multiplier-std-normal",
              om_O1 = "--O1 -foffset-multiplier-std-normal", om_Oexp = "--Oexperimental",
              base_Oexp = "--Oexperimental -fno-offset-multiplier-std-normal")
  pdir <- file.path(WORK, "probe")
  dir.create(pdir, showWarnings = FALSE)
  file.copy(file.path(WORK, "om_probe_main.hpp"), pdir, overwrite = TRUE)
  main_o <- file.path(pdir, "om_probe_main.o")
  mk <- function(target, extra = character()) {
    out <- suppressWarnings(system2("make", c("-C", cmdstan_path(), extra, target),
                                    stdout = TRUE, stderr = TRUE))
    st <- attr(out, "status")
    if (!is.null(st) && st != 0) {
      stop(sprintf("make %s failed:\n%s", target, paste(tail(out, 30), collapse = "\n")))
    }
  }
  mk(main_o)
  for (v in names(PFLAGS)) {
    hpp <- tempfile(fileext = ".hpp")
    system2(STANC, c(strsplit(PFLAGS[[v]], " ")[[1]], sprintf("--o=%s", hpp),
                     file.path(WORK, "om_nodata.stan")))
    cat(sprintf("%-9s std_normal_lpdf lines: %d\n", v, sum(grepl("std_normal_lpdf", readLines(hpp)))))
  }
  built <- parallel::mclapply(names(PFLAGS), function(v) {
    d <- file.path(pdir, v)
    dir.create(d, showWarnings = FALSE)
    w <- file.path(d, "stanc_wrap.sh")
    writeLines(c("#!/bin/bash", "shift", sprintf("exec %s %s \"$@\"", STANC, PFLAGS[[v]])), w)
    Sys.chmod(w, "0755")
    file.copy(file.path(WORK, "om_nodata.stan"), d, overwrite = TRUE)
    exe <- file.path(d, "om_nodata")
    tryCatch({
      mk(exe, c(sprintf("WINE=%s", w), sprintf("CMDSTAN_MAIN_O=%s", main_o)))
      exe
    }, error = function(e) conditionMessage(e))
  }, mc.cores = min(N_CORES, length(PFLAGS)), mc.preschedule = FALSE)
  exes <- unlist(built)
  names(exes) <- names(PFLAGS)
  bad <- !file.exists(exes)
  if (any(bad)) {
    cat(sprintf("probe build failed: %s\n%s\n", paste(names(exes)[bad], collapse = ", "),
                exes[bad][1]))
    quit(status = 1)
  }
  n_par <- as.integer(system2(exes[["base"]], stdin = "/dev/null", stdout = TRUE))
  cat(sprintf("om_nodata: %d unconstrained parameters\n", n_par))
  run_probe <- function(exe, uf) {
    out <- system2(exe, stdin = uf, stdout = TRUE)
    do.call(rbind, lapply(strsplit(out, " "), as.numeric))
  }
  set.seed(1102)
  rows <- list()
  for (k in 1:20) {
    uf <- tempfile()
    writeLines(sprintf("%.17g", rnorm(n_par, 0, 1.5)), uf)
    ref <- run_probe(exes[["base"]], uf)
    for (v in setdiff(names(PFLAGS), "base")) {
      x <- run_probe(exes[[v]], uf)
      for (r in 1:4) {
        gr <- ref[r, -(1:3)]
        rows[[length(rows) + 1]] <- data.frame(
          variant = v, point = k, propto = ref[r, 1], jacobian = ref[r, 2], lp = ref[r, 3],
          lp_diff = x[r, 3] - ref[r, 3],
          grad_rel = max(abs(x[r, -(1:3)] - gr) / pmax(1, abs(gr))),
          stringsAsFactors = FALSE)
      }
    }
    unlink(uf)
  }
  s4 <- do.call(rbind, rows)
  saveRDS(s4, file.path(OUT, "step4_probe.rds"))
  cat(sprintf("%-10s %6s %8s %11s %11s %10s\n", "variant", "propto", "jacobian", "min lpdiff",
              "max lpdiff", "grad rel"))
  for (v in setdiff(names(PFLAGS), "base")) {
    for (pr in c(1, 0)) {
      for (jc in c(1, 0)) {
        x <- s4[s4$variant == v & s4$propto == pr & s4$jacobian == jc, ]
        cat(sprintf("%-10s %6d %8d %11.4g %11.4g %10.2g\n", v, pr, jc, min(x$lp_diff),
                    max(x$lp_diff), max(x$grad_rel)))
      }
    }
  }
  cat(sprintf("expected propto lp shift: %.6g\n", -(4 * log(0.7) + 3 * log(2))))
}
cat("\nEXIT=0\n")
