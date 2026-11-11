## Step 3: Model 2 on the GeneEvolve G x Age trios from 01_simulate.py: every trio file in output/data (one per
## population and u level it exports; with the default independent design each u level has its own populations) and
## the pooled file of each mechanism x u. Each fit uses FTestValidity's testF (a multi-start free fit and an f = 0 fit),
## starting from the parents' equilibrium truth. The truth of every parameter is the parents' (adult) value; f's is .15.
## One CSV per fit in output/sim_fits/ (resumable), combined into output/sim_fits.csv: one row per fit x Model 2
## parameter (estimate, SE, truth), plus f_est, f_SE and the f = 0 LRT.
## Run with: Rscript 03_fit_sim.R [--workers 3] [--out output]   (--out: the output root holding data/; output_snpsets or
##           output_overlap for 06_simulate_overlap.py)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
getArg  <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else as.numeric(a[i + 1])
}
workers <- getArg("--workers", 3)
outArg  <- match("--out", commandArgs(trailingOnly = TRUE))
outRoot <- if (is.na(outArg)) "output" else commandArgs(trailingOnly = TRUE)[outArg + 1]   # output root (default output/)

# one thread per process, and one OpenMx thread per fit (the fitting scripts ask for parallel::detectCores())
Sys.setenv(OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(MASS); library(parallel) })
assignInNamespace("detectCores", function(...) 1L, ns = "parallel")
source(file.path(baseDir, "conditions.R"))

dataDir <- file.path(baseDir, outRoot, "data")
fitDir  <- file.path(baseDir, outRoot, "sim_fits"); dir.create(fitDir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(dataDir, "^Eq_(primary|genetic|social)_(rep[0-9]{3}|pooled)_u[0-9]{3}\\.tsv$")
if (length(files) == 0) stop("No datasets in ", dataDir, " - run 01_simulate.py first")
TRUTH <- lapply(setNames(names(MECHANISMS), names(MECHANISMS)), function(m) gxageSpec(m, F_TRUE, 0)$truth)

fitOne <- function(fn){
    tag <- sub("\\.tsv$", "", fn)
    out <- file.path(fitDir, paste0(tag, ".csv"))
    if (file.exists(out)) return(invisible(out))
    m    <- regmatches(tag, regexec("^Eq_([a-z]+)_(rep([0-9]{3})|pooled)_u([0-9]{3})$", tag))[[1]]
    mech <- m[2]
    repl <- if (m[3] == "pooled") NA_integer_ else as.integer(m[4])
    path <- file.path(dataDir, fn)
    t0 <- Sys.time()
    r  <- testF("M2_Cascade", path, multipliers("M2_Cascade", mech), TRUTH[[mech]])
    est <- attr(r, "params")
    if (is.null(est)) est <- data.frame(parameter = NA_character_, estimate = NA_real_, SE = NA_real_)
    id  <- data.frame(mechanism = mech, rep = repl, pooled = is.na(repl), u = as.integer(m[5]) / 100, f = F_TRUE,
                      n = nrow(fread(path, select = 1L)))
    res <- cbind(id[rep(1, nrow(est)), ], est, true = unname(TRUTH[[mech]][est$parameter]),
                 r[rep(1, nrow(est)), setdiff(names(r), "model")],
                 time_sec = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
    fwrite(res, out)
    cat(sprintf("%-26s n=%5d f=%.3f (SE %.3f) LRT=%6.2f status=%s (%s s)\n", tag, id$n, r$f_est, r$f_SE, r$LRT,
                r$status_free, res$time_sec[1]))
    invisible(out)
}

todo <- files[!file.exists(file.path(fitDir, sub("\\.tsv$", ".csv", files)))]
cat(length(files), "datasets,", length(todo), "to fit;", workers, "workers\n")
t0 <- Sys.time()
invisible(mclapply(todo, function(fn) tryCatch(fitOne(fn), error = function(e) cat(fn, "ERROR:", conditionMessage(e), "\n")),
                   mc.cores = workers, mc.preschedule = FALSE))
cat("Done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

all <- rbindlist(lapply(list.files(fitDir, "\\.csv$", full.names = TRUE), fread), fill = TRUE)
fwrite(all, file.path(baseDir, outRoot, "sim_fits.csv"))
cat("Wrote", nrow(all), "rows from", nrow(unique(all[, .(mechanism, rep, u)])), "fits to", file.path(outRoot, "sim_fits.csv"), "\n\n")
fr <- unique(all[, .(mechanism, rep, pooled, u, f_est, status_free)])
print(fr[pooled == FALSE, .(reps = .N, mean_f = round(mean(f_est, na.rm = TRUE), 4),
                            mc_se = round(sd(f_est, na.rm = TRUE) / sqrt(sum(!is.na(f_est))), 4)), by = .(mechanism, u)]
      [fr[pooled == TRUE, .(mechanism, u, pooled_f = round(f_est, 4))], on = .(mechanism, u)])
