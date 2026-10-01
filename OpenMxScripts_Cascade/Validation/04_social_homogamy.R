## Step 4: why are some social-homogamy fits off? Two diagnostics on the continuous models.
##
##  1. Per-population fits. Each GeneEvolve population (12 Eq, 36 DisEq per mechanism; ~4,000 trios each) is
##     fit on its own. Across populations this gives the sampling distribution of every estimate: is its mean
##     at the truth (bias), and does the spread match the SEs OpenMx reports (SE calibration)?
##     -> output/results/estimates_perpop.csv
##  2. Profile likelihood of mu. On the pooled data, mu is fixed on a grid around the truth and every other
##     parameter re-optimized; the rise in -2LL shows how sharply the data pin mu down under each mechanism.
##     -> output/results/profile_mu.csv
##
## Run with: Rscript 04_social_homogamy.R [--workers 10]   (resumable; finished fits are skipped)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
cascadeDir <- dirname(baseDir)
getArg <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else a[i + 1]
}
workers <- as.integer(getArg("--workers", 10))

suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(parallel) })
assignInNamespace("detectCores", function(...) 1L, ns = "parallel")   # one OpenMx thread per fit
source(file.path(baseDir, "conditions.R"))
cascadeValidationSetup(cascadeDir)

PERPOP_DESIGNS <- c("SameTrait_ObservedYPYM", "SameTrait_LatentYPYM", "DiffTrait_ObservedYPYM_FixedAParent")
for (regime in REGIMES) for (design in PERPOP_DESIGNS) source(scriptPath(cascadeDir, regime, design, FALSE))

partsDir <- file.path(baseDir, "output", "sim", "parts")
outDir   <- file.path(baseDir, "output", "results", "perpop"); dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

## ---------------- 1. per-population fits ----------------
info  <- fread(file.path(baseDir, "output", "sim", "rep_info.csv"))
tasks <- CJ(regime = REGIMES, mech = names(MECHANISMS), design = PERPOP_DESIGNS, rep = sort(unique(info$rep)))
tasks <- tasks[info[, .(regime, mech = mechanism, rep)], on = .(regime, mech, rep), nomatch = 0]

fitPop <- function(i){
    t <- tasks[i]
    out <- file.path(outDir, sprintf("%s_%s_%s_rep%03d.csv", t$regime, t$mech, t$design, t$rep))
    if (file.exists(out)) return(invisible(out))
    spec <- designSpec(conditionTruth(t$regime, t$mech), t$design, FALSE)
    data_path <- file.path(partsDir, sprintf("%s_%s_rep%03d_%s.tsv", t$regime, t$mech, t$rep, spec$dataset))
    res <- tryCatch({
        s <- suppressMessages(do.call(fitFunctionName(t$regime, t$design, FALSE),
                                      c(list(data_path = data_path), spec$args, list(extraTries = 15, exhaustive = TRUE))))
        p <- s$parameters
        data.table(parameter = p$name, estimate = p$Estimate, SE = p$Std.Error, true = unname(spec$truth[p$name]),
                   statusCode = as.character(s$statusCode)[1], n = s$numObs)
    }, error = function(e) data.table(parameter = NA, estimate = NA, SE = NA, true = NA, statusCode = "ERROR", n = NA))
    fwrite(cbind(regime = t$regime, mechanism = t$mech, design = t$design, rep = t$rep, res), out)
    invisible(out)
}
cat(nrow(tasks), "per-population fits,", workers, "workers\n")
invisible(mclapply(seq_len(nrow(tasks)), fitPop, mc.cores = workers, mc.preschedule = FALSE))
perpop <- rbindlist(lapply(list.files(outDir, "\\.csv$", full.names = TRUE), fread), fill = TRUE)
fwrite(perpop, file.path(baseDir, "output", "results", "estimates_perpop.csv"))
cat("per-population fits:", uniqueN(perpop[, paste(regime, mechanism, design, rep)]), "\n")

## ---------------- 2. profile likelihood of mu (pooled data, SameTrait_ObservedYPYM) ----------------
# Build the model exactly as the fitting script does, but return it instead of fitting it.
buildModel <- function(regime, design, args, data_path){
    src <- paste(readLines(scriptPath(cascadeDir, regime, design, FALSE)), collapse = "\n")
    src <- sub("fitModel1 <- mxTryHard\\((?s).*", "return(Model1)\n}\n", src, perl = TRUE)
    env <- new.env(); eval(parse(text = src), envir = env)
    do.call(get(fitFunctionName(regime, design, FALSE), envir = env), c(list(data_path = data_path), args))
}
GRID <- c(.7, .8, .9, .95, 1, 1.05, 1.1, 1.2, 1.3)
ptasks <- CJ(regime = REGIMES, mech = names(MECHANISMS), g = GRID)
profOne <- function(i){
    t <- ptasks[i]
    tr <- conditionTruth(t$regime, t$mech); spec <- designSpec(tr, "SameTrait_ObservedYPYM", FALSE)
    data_path <- file.path(baseDir, "output", "data", sprintf("%s_%s_SameTrait_Observed.tsv", t$regime, t$mech))
    f <- tryCatch({
        m <- buildModel(t$regime, "SameTrait_ObservedYPYM", spec$args, data_path)
        # start every structural parameter at the truth (the script's generic start values can give a
        # non-positive-definite starting matrix once mu is fixed away from its estimate)
        fl <- intersect(names(omxGetParameters(m)), names(spec$truth))
        m <- omxSetParameters(m, labels = fl, values = unname(spec$truth[fl]))
        m <- omxSetParameters(m, labels = "mu11", free = FALSE, values = t$g * tr$q$mu)
        suppressMessages(mxTryHard(m, extraTries = 10, exhaustive = TRUE, silent = TRUE, OKstatuscodes = c(0, 1)))
    }, error = function(e) conditionMessage(e))
    ok <- inherits(f, "MxModel") && length(f$output$Minus2LogLikelihood) == 1
    data.table(regime = t$regime, mechanism = t$mech, mu_ratio = t$g, mu = t$g * tr$q$mu,
               minus2LL = if (ok) f$output$Minus2LogLikelihood else NA_real_,
               status = if (ok) paste(f$output$status$code, collapse = "") else if (inherits(f, "MxModel")) paste("FAILED", f$output$status$code) else paste("ERROR:", substr(f, 1, 200)))
}
profFile <- file.path(baseDir, "output", "results", "profile_mu.csv")
prof <- if (file.exists(profFile)) fread(profFile)[!is.na(minus2LL)] else NULL
todo <- if (is.null(prof)) seq_len(nrow(ptasks)) else which(!paste(ptasks$regime, ptasks$mech, ptasks$g) %in% paste(prof$regime, prof$mechanism, prof$mu_ratio))
if (length(todo)){
    cat(length(todo), "profile fits\n")
    res <- mclapply(todo, profOne, mc.cores = workers, mc.preschedule = FALSE)
    bad <- !vapply(res, is.data.frame, TRUE)
    if (any(bad)) print(lapply(res[bad], function(e) substr(as.character(e), 1, 300)))
    prof <- rbindlist(c(list(prof), res[!bad]), fill = TRUE)
    setorder(prof, regime, mechanism, mu_ratio)
    fwrite(prof, profFile)
}
cat("done\n")
