## Step 7: the corrected ORIGINAL DisEq script (OpenMxScripts/Disequilibrium/UniSEMPGS_DisEq_SameTrait_ObservedYPYM.R)
## on the GeneEvolve primary-AM DisEq populations. Primary AM is the only mechanism the original model covers.
##   - each of the 36 populations on its own: is the mean estimate at the truth, and do the reported SEs match
##     the spread across populations?
##   - the 36 populations pooled
## -> output/results/original_diseq_check.csv
##
## Run with: Rscript 07_original_diseq_check.R [--workers 4]

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
cascadeDir <- dirname(baseDir)
getArg <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else a[i + 1]
}
workers <- as.integer(getArg("--workers", 4))

suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(parallel) })
assignInNamespace("detectCores", function(...) 1L, ns = "parallel")
source(file.path(baseDir, "conditions.R"))
cascadeValidationSetup(cascadeDir)
source(file.path(dirname(cascadeDir), "OpenMxScripts", "Disequilibrium", "UniSEMPGS_DisEq_SameTrait_ObservedYPYM.R"))

spec <- designSpec(conditionTruth("DisEq", "primary"), "SameTrait_ObservedYPYM", FALSE)
outDir <- file.path(baseDir, "output", "original_diseq"); dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
files <- sort(list.files(file.path(baseDir, "output", "sim", "parts"), "^DisEq_primary_rep\\d+_SameTrait_Observed\\.tsv$", full.names = TRUE))
reps  <- as.integer(sub(".*_rep(\\d+)_.*", "\\1", basename(files)))
jobs  <- c(lapply(seq_along(files), function(i) list(fit = "population", rep = reps[i], path = files[i], tries = 15)),
           list(list(fit = "pooled", rep = NA_integer_, path = file.path(baseDir, "output", "data", "DisEq_primary_SameTrait_Observed.tsv"), tries = 30)))

fitJob <- function(j){
    out <- file.path(outDir, if (j$fit == "pooled") "pooled.csv" else sprintf("rep%03d.csv", j$rep))
    if (!file.exists(out)){
        res <- tryCatch({
            s <- suppressMessages(fitUniSEMPGS_DisEq_SameTrait_ObservedYPYM(j$path, extraTries = j$tries, exhaustive = TRUE))
            p <- s$parameters
            data.table(parameter = p$name, estimate = p$Estimate, SE = p$Std.Error, true = unname(spec$truth[p$name]),
                       statusCode = as.character(s$statusCode)[1], n = s$numObs)
        }, error = function(e) data.table(parameter = NA, estimate = NA, SE = NA, true = NA, statusCode = "ERROR", n = NA))
        fwrite(res, out)
    }
    cbind(fit = j$fit, rep = j$rep, fread(out))
}
cat(length(jobs), "fits\n")
res <- rbindlist(mclapply(jobs, fitJob, mc.cores = workers, mc.preschedule = FALSE), fill = TRUE)
fwrite(res, file.path(baseDir, "output", "results", "original_diseq_check.csv"))
cat("done\n")
