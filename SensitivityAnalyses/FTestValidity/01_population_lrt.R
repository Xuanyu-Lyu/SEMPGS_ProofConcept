## Step 1: the population value of each f test. Every condition's model-implied 7x7 covariance is turned into a
## dataset whose sample covariance equals it exactly (mvrnorm, empirical = TRUE, n = 10,000), and the three models
## are tested on it. The resulting LRT is the noncentrality lambda of that test at n = 10,000; it grows in proportion
## to n, so the expected rejection rate at any n is P(chi2_1(lambda * n / 10,000) > 3.84). The free-f estimate is
## the model's population (pseudo-true) value of f.
## Run with: Rscript 01_population_lrt.R            (writes output/population_lrt.csv)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
# at most 3 threads: the fitting scripts ask OpenMx for parallel::detectCores() threads
Sys.setenv(OMP_NUM_THREADS = 3)
suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(MASS); library(parallel) })
assignInNamespace("detectCores", function(...) 3L, ns = "parallel")
source(file.path(baseDir, "conditions.R"))

N_REF <- 10000
outDir <- file.path(baseDir, "output"); dir.create(file.path(outDir, "tmp"), recursive = TRUE, showWarnings = FALSE)
res <- rbindlist(lapply(seq_len(nrow(GRID)), function(i){
    row  <- GRID[i, ]
    spec <- conditionSpec(row)
    dat  <- mvrnorm(N_REF, rep(0, 7), spec$CM, empirical = TRUE); colnames(dat) <- colnames(spec$CM)
    paths <- writeData(dat, file.path(outDir, "tmp", paste0("pop_", row$cond)))
    r <- rbindlist(lapply(MODELS, function(m) testF(m, paths[[m]], multipliers(m, row$mating), spec$truth, extraTries = 10)))
    unlink(unique(paths))
    q <- spec$q
    r <- cbind(row[rep(1, nrow(r)), c("cond", "mating", "latent", "f")], r,
               f_true = row$f, ic_true = q$ic, gt_true = q$gt, mu_true = q$mu, VF_true = q$VF, n_ref = N_REF)
    cat(sprintf("%-24s %s\n", row$cond, paste(sprintf("%s: f=%.3f lambda=%.2f", r$model, r$f_est, r$LRT), collapse = " | ")))
    r
}))
res[, power_at_nref := pchisq(qchisq(.95, 1), 1, ncp = pmax(LRT, 0), lower.tail = FALSE)]
fwrite(res, file.path(outDir, "population_lrt.csv"))
cat("Wrote", nrow(res), "rows to output/population_lrt.csv\n")
