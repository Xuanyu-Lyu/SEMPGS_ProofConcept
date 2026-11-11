## Step 2: the population (pseudo-true) value of Model 2's f under G x Age. For every mating type x f x u, the G x Age
## population covariance (conditions.R: gxageSpec) is turned into a dataset whose sample covariance equals it exactly
## (mvrnorm, empirical = TRUE, n = 10,000), and Model 2 is fit to it with FTestValidity's testF (a multi-start free fit
## and an f = 0 fit). The free-fit estimates are Model 2's pseudo-true values, so f_est - f is the bias of f without
## sampling noise. LRT is the noncentrality of the f = 0 test at n = 10,000, and misfit_chi2 is Model 2's chi-square
## against the saturated model at n = 10,000 (whether the global fit reveals the G x Age misspecification). Both grow
## in proportion to n. As a check against local optima, Model 2 is also fit with the script's own exhaustive search
## from its default starting values (m2LL_check, f_check); it should not beat testF's -2LL.
## Writes
##   output/truth/gxage_truth.csv   analytic offspring moments and parental values per cell (also used by 04's sim check)
##   output/population_fits.csv     one row per cell x Model 2 parameter (estimate, SE, truth), plus the f-test columns
## Run with: Rscript 02_population_fits.R [--workers 3]   (resumable: finished cells in output/population/ are skipped)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
getArg  <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else as.numeric(a[i + 1])
}
workers <- getArg("--workers", 3)

# one thread per process, and one OpenMx thread per fit (the fitting scripts ask for parallel::detectCores())
Sys.setenv(OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(MASS); library(parallel) })
assignInNamespace("detectCores", function(...) 1L, ns = "parallel")
source(file.path(baseDir, "conditions.R"))

N_REF <- 10000
CELLS <- expand.grid(mating = c("none", "primary", "genetic", "social"), f = c(0, F_TRUE), u = U_VALUES,
                     stringsAsFactors = FALSE)
CELLS$cell <- with(CELLS, sprintf("%s_f%03d_u%03d", mating, round(100 * f), round(100 * u)))

outDir  <- file.path(baseDir, "output")
cellDir <- file.path(outDir, "population")
tmpDir  <- file.path(outDir, "tmp")
for (d in c(cellDir, tmpDir, file.path(outDir, "truth"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---------------- analytic truth ----------------
truth <- rbindlist(lapply(seq_len(nrow(CELLS)), function(i){
    cl <- CELLS[i, ]
    s  <- gxageSpec(cl$mating, cl$f, cl$u); q <- s$q; o <- s$o
    v  <- c(r = s$r, VG0 = s$VG0, VY_o = o$VY_o, thetaT = o$thetaT_o, thetaNT = o$thetaNT_o, Yo_Yp = o$Yo_Yp_o,
            delta_o = o$delta_o, a_o = o$a_o, VE_o = o$VE_o,
            VY = q$VY, delta = q$delta, a = q$a, VE = q$VE, f = q$f, mu = q$mu, Omega = q$Omega, Gamma = q$Gamma,
            w = q$w, v = q$v, VF = q$VF, gc = q$gc, gt = q$gt, ic = q$ic, Yp_Ym = q$Yp_Ym, Yp_PGSm = q$Yp_PGSm)
    data.table(cl[rep(1, length(v)), ], quantity = names(v), value = unname(v))
}))
fwrite(truth, file.path(outDir, "truth", "gxage_truth.csv"))

## ---------------- Model 2 on the exact population covariance ----------------
runCell <- function(i){
    cl  <- CELLS[i, ]
    out <- file.path(cellDir, paste0(cl$cell, ".csv"))
    if (file.exists(out)) return(invisible(out))
    s <- gxageSpec(cl$mating, cl$f, cl$u)
    set.seed(2026)
    dat <- mvrnorm(N_REF, rep(0, 7), s$CM, empirical = TRUE); colnames(dat) <- colnames(s$CM)
    path <- file.path(tmpDir, paste0(cl$cell, ".tsv"))
    fwrite(as.data.frame(dat)[, COLS7], path, sep = "\t")
    t0 <- Sys.time()
    r  <- testF("M2_Cascade", path, multipliers("M2_Cascade", cl$mating), s$truth, extraTries = 10)
    # local-optimum check: the script's own exhaustive search from its default starting values
    chk <- NULL
    invisible(capture.output(chk <- tryCatch(suppressMessages(suppressWarnings(do.call(fitUniSEMPGS_Cascade_SameTrait_ObservedYPYM,
        c(list(data_path = path), multipliers("M2_Cascade", cl$mating), list(extraTries = 30, exhaustive = TRUE))))),
        error = function(e) NULL)))
    unlink(path)
    est <- attr(r, "params")
    if (is.null(est)) est <- data.frame(parameter = NA_character_, estimate = NA_real_, SE = NA_real_)
    # saturated model: 7 means and 28 covariances at the data's ML covariance (= CM * (n - 1)/n)
    m2LL_sat <- N_REF * (7 * log(2 * pi) + as.numeric(determinant(s$CM * (N_REF - 1) / N_REF)$modulus) + 7)
    df_model <- if (is.null(attr(r, "df"))) NA else attr(r, "df")
    res <- cbind(cl[rep(1, nrow(est)), ], est, true = unname(s$truth[est$parameter]),
                 r[rep(1, nrow(est)), setdiff(names(r), "model")],
                 misfit_chi2 = r$m2LL_free - m2LL_sat, misfit_df = df_model - (7 * N_REF - 35), n_ref = N_REF,
                 m2LL_check = if (is.null(chk)) NA else chk$Minus2LogLikelihood,
                 f_check = if (is.null(chk)) NA else chk$parameters$Estimate[chk$parameters$name == "f11"],
                 time_sec = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
    fwrite(res, out)
    cat(sprintf("%-18s f=%.3f bias=%+.4f LRT=%7.2f misfit=%7.2f | check: f=%.3f d(-2LL)=%+.4f (%s s)\n", cl$cell,
                r$f_est, r$f_est - cl$f, r$LRT, res$misfit_chi2[1], res$f_check[1], res$m2LL_check[1] - r$m2LL_free,
                res$time_sec[1]))
    invisible(out)
}

todo <- which(!file.exists(file.path(cellDir, paste0(CELLS$cell, ".csv"))))
cat(nrow(CELLS), "cells,", length(todo), "to fit;", workers, "workers\n")
t0 <- Sys.time()
invisible(mclapply(todo, function(i) tryCatch(runCell(i), error = function(e) cat(CELLS$cell[i], "ERROR:", conditionMessage(e), "\n")),
                   mc.cores = workers, mc.preschedule = FALSE))
cat("Done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

all <- rbindlist(lapply(list.files(cellDir, "\\.csv$", full.names = TRUE), fread), fill = TRUE)
fwrite(all, file.path(outDir, "population_fits.csv"))
cat("Wrote", nrow(all), "rows from", uniqueN(all$cell), "cells to output/population_fits.csv\n\nBias of f (pseudo-true f - f):\n")
fb <- unique(all[, .(mating, f, u, bias = round(f_est - f, 4))])
print(dcast(fb, mating + f ~ u, value.var = "bias"))
