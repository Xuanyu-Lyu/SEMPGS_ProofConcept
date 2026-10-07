## Step 2: Monte Carlo. For every condition, independent trio datasets are drawn from the model-implied 7x7 covariance
## (mvrnorm, sampling variability included) and f is tested with all three models (conditions.R: testF).
##   f = 0 conditions (type I error): 50 datasets at n = 2,000 and 50 at n = 10,000
##   f > 0 conditions:                50 datasets at n = 10,000
## The Monte Carlo is a finite-sample check: the precise expected rejection rates come from the population
## noncentrality in 01_population_lrt.R.
## One CSV per dataset in output/mc/fits/ (resumable: finished datasets are skipped); all rows are then combined into
## output/mc_results.csv.
## Run with: Rscript 02_monte_carlo.R [--workers 3] [--scale 1]   (--scale multiplies the numbers of datasets)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
getArg  <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else as.numeric(a[i + 1])
}
workers <- getArg("--workers", 3)
scale   <- getArg("--scale", 1)

# one thread per process: without this, OpenMP/BLAS threads in every worker oversubscribe the CPU
Sys.setenv(OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(MASS); library(parallel) })
# one OpenMx thread per fit: the fitting scripts ask for parallel::detectCores() threads, which oversubscribes the CPU
# when several fits run at once (as in OpenMxScripts_Cascade/Validation/02_fit_models.R)
assignInNamespace("detectCores", function(...) 1L, ns = "parallel")
source(file.path(baseDir, "conditions.R"))

fitDir <- file.path(baseDir, "output", "mc", "fits"); dir.create(fitDir, recursive = TRUE, showWarnings = FALSE)
tmpDir <- file.path(baseDir, "output", "mc", "tmp");  dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)

design <- rbind(data.frame(f = 0, n = 2000, reps = 50), data.frame(f = 0, n = 10000, reps = 50),
                data.frame(f = F_VALUES[F_VALUES > 0], n = 10000, reps = 50))
design$reps <- round(design$reps * scale)
tasks <- rbindlist(lapply(seq_len(nrow(GRID)), function(i){
    d <- design[design$f == GRID$f[i], ]
    rbindlist(lapply(seq_len(nrow(d)), function(j) data.table(gi = i, cond = GRID$cond[i], n = d$n[j], rep = seq_len(d$reps[j]))))
}))
tasks[, out := file.path(fitDir, sprintf("%s_n%05d_rep%03d.csv", cond, n, rep))]
todo <- which(!file.exists(tasks$out))
cat(nrow(tasks), "datasets in total,", length(todo), "to do;", workers, "workers\n")

SPECS <- lapply(seq_len(nrow(GRID)), function(i) conditionSpec(GRID[i, ]))

runOne <- function(k){
    t   <- tasks[k]
    row <- GRID[t$gi, ]
    spec <- SPECS[[t$gi]]
    set.seed(100000 * t$gi + 1000 * (t$n == 10000) + t$rep)
    dat <- mvrnorm(t$n, rep(0, 7), spec$CM); colnames(dat) <- colnames(spec$CM)
    paths <- writeData(dat, file.path(tmpDir, sprintf("%s_n%05d_rep%03d", t$cond, t$n, t$rep)))
    t0 <- Sys.time()
    r <- rbindlist(lapply(MODELS, function(m)
        tryCatch(testF(m, paths[[m]], multipliers(m, row$mating), spec$truth),
                 error = function(e) data.frame(model = m, f_est = NA, f_SE = NA, m2LL_free = NA, m2LL_f0 = NA,
                                                LRT = NA, p = NA, status_free = NA, status_f0 = NA, reseeded = NA))))
    unlink(unique(paths))
    r <- cbind(data.table(cond = t$cond, mating = row$mating, latent = row$latent, f = row$f, n = t$n, rep = t$rep), r,
               time_sec = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
    fwrite(r, t$out)
    invisible(NULL)
}

t0 <- Sys.time()
invisible(mclapply(todo, function(k) tryCatch(runOne(k), error = function(e) NULL),
                   mc.cores = workers, mc.preschedule = FALSE))
cat("Done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

files <- tasks$out[file.exists(tasks$out)]
all <- rbindlist(lapply(files, fread), fill = TRUE)
fwrite(all, file.path(baseDir, "output", "mc_results.csv"))
cat("Wrote", nrow(all), "rows from", length(files), "datasets to output/mc_results.csv\n")
print(all[, .(datasets = .N, failed = sum(is.na(p)), reject_rate = round(mean(p < .05, na.rm = TRUE), 3)),
          by = .(f, n, model)][order(f, n, model)])
