## Step 2: fit all 20 Cascade scripts (5 designs x {Eq, DisEq} x {continuous, binary}) under each of the
## three mating mechanisms, and join every estimate to its iterative-math truth by OpenMx label.
##
##   --mode mvn   exact-covariance coding check (as in the original test_*Models.R): data are drawn from the
##                model-implied covariance matrix of the iterative math (mvrnorm empirical = TRUE; continuous
##                n = 32,000, binary n = 50,000 dichotomized at fixed thresholds). A correctly coded,
##                identified continuous model must recover the truth to near machine precision.
##   --mode sim   the pooled GeneEvolve trio data written by 01_simulate.py (output/data/).
##   --binary 0|1 restrict to continuous (0) or binary (1) models; default both.
##   --workers K  fits in parallel (mclapply). Resumable: finished fits in output/results/fits/ are skipped.
##   --threads T  OpenMx threads per fit (default: cores / workers)
## Run with: Rscript 02_fit_models.R --mode mvn [--binary 0] [--workers 4]

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
cascadeDir <- dirname(baseDir)
getArg <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else a[i + 1]
}
mode    <- getArg("--mode", "mvn")
binSel  <- getArg("--binary", NA)
workers <- as.integer(getArg("--workers", 1))
stopifnot(mode %in% c("mvn", "sim"))

suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(MASS); library(parallel) })
# The fitting scripts ask OpenMx for parallel::detectCores() threads; with several fits in parallel that
# oversubscribes the CPU, so cap what detectCores() returns in this session (the scripts are not changed).
threads <- as.integer(getArg("--threads", max(1, parallel::detectCores() %/% workers)))
assignInNamespace("detectCores", function(...) threads, ns = "parallel")
source(file.path(baseDir, "conditions.R"))
cascadeValidationSetup(cascadeDir)
for (regime in REGIMES) for (design in DESIGNS) for (binary in c(FALSE, TRUE)) source(scriptPath(cascadeDir, regime, design, binary))

fitDir  <- file.path(baseDir, "output", "results", "fits"); dir.create(fitDir, recursive = TRUE, showWarnings = FALSE)
mvnDir  <- file.path(baseDir, "output", "mvn_data");        dir.create(mvnDir, recursive = TRUE, showWarnings = FALSE)
simDir  <- file.path(baseDir, "output", "data")

tasks <- expand.grid(regime = REGIMES, mech = names(MECHANISMS), design = DESIGNS,
                     binary = if (is.na(binSel)) c(FALSE, TRUE) else as.logical(as.integer(binSel)),
                     stringsAsFactors = FALSE)
cat(nrow(tasks), "fits, mode =", mode, ",", workers, "worker(s)\n")

mvnData <- function(tag, spec, binary){
    tsv <- file.path(mvnDir, paste0(tag, ".tsv"))
    if (file.exists(tsv)) return(tsv)
    CM <- spec$CM
    stopifnot(max(abs(CM - t(CM))) < 1e-10, !is.null(tryCatch(chol(CM), error = function(e) NULL)))
    set.seed(2026)
    dat <- as.data.frame(mvrnorm(n = if (binary) 5e4 else 32000, mu = rep(0, ncol(CM)), Sigma = CM, empirical = TRUE))
    if (binary) for (v in names(spec$thresholds)) dat[[v]] <- as.integer(dat[[v]] > spec$thresholds[[v]])
    fwrite(dat, tsv, sep = "\t")
    tsv
}

fitOne <- function(i){
    t <- tasks[i, ]
    tag <- sprintf("%s_%s_%s_%s_%s", mode, t$regime, t$mech, t$design, if (t$binary) "binary" else "continuous")
    out <- file.path(fitDir, paste0(tag, ".csv"))
    if (file.exists(out) || file.exists(sub("\\.csv$", "_origvals.csv", out))) return(invisible(out))
    tr   <- conditionTruth(t$regime, t$mech)
    spec <- designSpec(tr, t$design, t$binary)
    base_used <- "BASE"
    if (mode == "mvn" && min(eigen(spec$CM, only.values = TRUE)$values) <= 0){
        # the model-implied matrix at BASE is not positive definite, so no data can be drawn from it: record
        # that, and run the coding check at the original OpenMxScripts test values instead
        base_used <- sprintf("BASE_ORIGINAL_TEST (BASE-implied matrix not PD, min eigenvalue %.4f)",
                             min(eigen(spec$CM, only.values = TRUE)$values))
        tr   <- conditionTruth(t$regime, t$mech, base = BASE_ORIGINAL_TEST)
        spec <- designSpec(tr, t$design, t$binary)
        tag  <- paste0(tag, "_origvals")
    }
    truth <- spec$truth
    if (mode == "mvn"){
        data_path <- mvnData(sub("^mvn_", "", tag), spec, t$binary)
        truth <- c(truth, setNames(rep(0, 7), c("meanYp1","meanYm1","meanYo1","meanTp1","meanNTp1","meanTm1","meanNTm1")))
    } else {
        data_path <- file.path(simDir, sprintf("%s_%s_%s%s.tsv", t$regime, t$mech, spec$dataset, if (t$binary) "_binary" else ""))
    }
    fitArgs <- c(list(data_path = data_path), spec$args,
                 # continuous fits: exhaustive search (the coding check found local optima without it);
                 # binary fits are single-threaded ordinal FIML (10-30 min per attempt), so they use the
                 # scripts' own default search: up to 30 retries, stopping at the first converged solution
                 list(extraTries = 30, exhaustive = !t$binary))
    t0 <- Sys.time()
    res <- tryCatch({
        s <- do.call(fitFunctionName(t$regime, t$design, t$binary), fitArgs)
        p <- s$parameters
        data.frame(parameter = p$name, estimate = p$Estimate, SE = p$Std.Error, true = unname(truth[p$name]),
                   statusCode = as.character(s$statusCode)[1], minus2LL = s$Minus2LogLikelihood, n = s$numObs,
                   error = "", stringsAsFactors = FALSE)
    }, error = function(e) data.frame(parameter = NA, estimate = NA, SE = NA, true = NA, statusCode = "ERROR",
                                      minus2LL = NA, n = NA, error = conditionMessage(e), stringsAsFactors = FALSE))
    res <- cbind(mode = mode, base = base_used, regime = t$regime, mechanism = t$mech, design = t$design,
                 trait = if (t$binary) "binary" else "continuous", res,
                 time_sec = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
    out <- file.path(fitDir, paste0(tag, ".csv"))
    fwrite(res, out)
    struct <- !is.na(res$true) & !grepl("^mean", res$parameter)
    cat(sprintf("%-78s status=%-6s max|est-true|=%.2e  %5.0fs\n", tag, res$statusCode[1],
                if (any(struct)) max(abs(res$estimate[struct] - res$true[struct])) else NA, res$time_sec[1]))
    invisible(out)
}

if (workers > 1) invisible(mclapply(seq_len(nrow(tasks)), fitOne, mc.cores = workers, mc.preschedule = FALSE)) else
    invisible(lapply(seq_len(nrow(tasks)), fitOne))

all <- rbindlist(lapply(list.files(fitDir, paste0("^", mode, "_.*\\.csv$"), full.names = TRUE), fread), fill = TRUE)
fwrite(all, file.path(baseDir, "output", "results", paste0("estimates_", mode, ".csv")))
cat("\nWrote", nrow(all), "rows from", uniqueN(all[, paste(regime, mechanism, design, trait)]), "fits\n")
print(all[, .(fits = uniqueN(paste(regime, mechanism, design))), by = .(trait, statusCode)])
