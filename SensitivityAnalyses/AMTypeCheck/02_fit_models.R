## Fit the equilibrium SameTrait_ObservedYPYM model (which assumes phenotypic AM) to every trio
## dataset produced by 01_simulate.py. Resumable: fits already in output/results/fits/ are skipped.
## Run from anywhere:  Rscript 02_fit_models.R [--workers K] [--extraTries N]
##   --workers K   fit K datasets in parallel via mclapply (each fit still asks OpenMx for all cores,
##                 so K > 1 oversubscribes; try 3-4 and compare the per-fit seconds printed below)
##   --extraTries  mxTryHard restarts per fit (default 8, as in test_EquilibriumModels.R)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
getArg  <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else as.numeric(a[i + 1])
}
workers    <- getArg("--workers", 1)
extraTries <- getArg("--extraTries", 8)

suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(parallel) })
source(file.path(baseDir, "..", "..", "OpenMxScripts", "Equilibrium", "UniSEMPGS_SameTrait_ObservedYPYM.R"))

dataDir <- file.path(baseDir, "output", "data")
fitDir  <- file.path(baseDir, "output", "results", "fits")
dir.create(fitDir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(dataDir, "^rep[0-9]+_[A-Za-z]+\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No datasets in ", dataDir, " - run 01_simulate.py first")
cat(length(files), "datasets;", workers, "worker(s); extraTries =", extraTries, "\n")

fitOne <- function(f){
    tag <- sub("\\.tsv$", "", basename(f))
    out <- file.path(fitDir, paste0("fit_", tag, ".csv"))
    if (file.exists(out)) return(invisible(out))
    m    <- regmatches(tag, regexec("^rep([0-9]+)_([A-Za-z]+)$", tag))[[1]]
    rep  <- as.integer(m[2]); cond <- m[3]
    t0   <- Sys.time()
    res  <- tryCatch({
        s <- fitUniSEMPGS_SameTrait_ObservedYPYM(f, extraTries = extraTries)
        data.frame(rep = rep, condition = cond,
                   parameter = s$parameters$name, estimate = s$parameters$Estimate, SE = s$parameters$Std.Error,
                   statusCode = as.character(s$statusCode)[1], minus2LL = s$Minus2LogLikelihood, n = s$numObs,
                   error = "", stringsAsFactors = FALSE)
    }, error = function(e)
        data.frame(rep = rep, condition = cond, parameter = NA, estimate = NA, SE = NA,
                   statusCode = "ERROR", minus2LL = NA, n = NA, error = conditionMessage(e), stringsAsFactors = FALSE))
    res$time_sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    write.csv(res, out, row.names = FALSE)
    cat(sprintf("%-22s status=%-12s %5.0fs\n", tag, res$statusCode[1], res$time_sec[1]))
    invisible(out)
}

if (workers > 1) invisible(mclapply(files, fitOne, mc.cores = workers)) else invisible(lapply(files, fitOne))

all <- rbindlist(lapply(list.files(fitDir, "^fit_.*\\.csv$", full.names = TRUE), fread), fill = TRUE)
fwrite(all, file.path(baseDir, "output", "results", "estimates.csv"))
cat("\nWrote", nrow(all), "rows from", uniqueN(all[, paste(rep, condition)]), "fits to output/results/estimates.csv\n")
print(all[, .(fits = uniqueN(rep)), by = .(condition, statusCode)])
