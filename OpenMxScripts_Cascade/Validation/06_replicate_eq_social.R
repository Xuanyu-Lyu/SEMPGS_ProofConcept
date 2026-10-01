## Step 6 (fits): the replication populations of 06_replicate_eq_social.py, fit with SameTrait_ObservedYPYM
##   - each population on its own (same settings as 04_social_homogamy.R)  -> output/replication/perpop/*.csv
##   - the 24 populations pooled                                              -> output/replication/pooled.csv
## Both are collected in output/results/replication_eq_social.csv.
##
## Run with: Rscript 06_replicate_eq_social.R [--workers 8]

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
cascadeDir <- dirname(baseDir)
getArg <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else a[i + 1]
}
workers <- as.integer(getArg("--workers", 8))

suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(parallel) })
assignInNamespace("detectCores", function(...) 1L, ns = "parallel")
source(file.path(baseDir, "conditions.R"))
cascadeValidationSetup(cascadeDir)

REGIME <- "Eq"; MECH <- "social"; DESIGN <- "SameTrait_ObservedYPYM"
source(scriptPath(cascadeDir, REGIME, DESIGN, FALSE))
repDir  <- file.path(baseDir, "output", "replication")
partsDir <- file.path(repDir, "parts")
outDir  <- file.path(repDir, "perpop"); dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
spec <- designSpec(conditionTruth(REGIME, MECH), DESIGN, FALSE)

fitData <- function(data_path, out, extraTries){
    if (file.exists(out)) return(invisible(out))
    res <- tryCatch({
        s <- suppressMessages(do.call(fitFunctionName(REGIME, DESIGN, FALSE),
                                      c(list(data_path = data_path), spec$args, list(extraTries = extraTries, exhaustive = TRUE))))
        p <- s$parameters
        data.table(parameter = p$name, estimate = p$Estimate, SE = p$Std.Error, true = unname(spec$truth[p$name]),
                   statusCode = as.character(s$statusCode)[1], n = s$numObs)
    }, error = function(e) data.table(parameter = NA, estimate = NA, SE = NA, true = NA, statusCode = "ERROR", n = NA))
    fwrite(res, out)
    invisible(out)
}

files <- sort(list.files(partsDir, sprintf("^%s_%s_rep\\d+_SameTrait_Observed\\.tsv$", REGIME, MECH), full.names = TRUE))
reps  <- as.integer(sub(".*_rep(\\d+)_.*", "\\1", basename(files)))
cat(length(files), "replication populations\n")

# pooled data of the replication populations
pooled <- file.path(repDir, sprintf("%s_%s_SameTrait_Observed_reps%d-%d.tsv", REGIME, MECH, min(reps), max(reps)))
if (!file.exists(pooled)) fwrite(rbindlist(lapply(files, fread)), pooled, sep = "\t")

jobs <- c(lapply(seq_along(files), function(i) list(path = files[i], out = file.path(outDir, sprintf("rep%03d.csv", reps[i])), tries = 15)),
          list(list(path = pooled, out = file.path(repDir, "pooled.csv"), tries = 30)))
invisible(mclapply(jobs, function(j) fitData(j$path, j$out, j$tries), mc.cores = workers, mc.preschedule = FALSE))

res <- rbind(
    rbindlist(lapply(seq_along(reps), function(i) cbind(fit = "population", rep = reps[i], fread(file.path(outDir, sprintf("rep%03d.csv", reps[i])))))),
    cbind(fit = "pooled", rep = NA_integer_, fread(file.path(repDir, "pooled.csv"))))
fwrite(res, file.path(baseDir, "output", "results", "replication_eq_social.csv"))
cat("done\n")
