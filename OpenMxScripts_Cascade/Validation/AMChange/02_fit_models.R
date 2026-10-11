## Step 2: fit the twelve AMChange scripts (6 designs x {continuous, binary}) and join every estimate to its truth by
## OpenMx label.
##
##   --mode mvn        exact-covariance coding check, as in ../02_fit_models.R: data drawn from the model-implied
##                     covariance matrix of the iterative math (mvrnorm empirical = TRUE; continuous n = 32,000, binary
##                     n = 50,000 dichotomized at fixed thresholds). A correctly coded, identified continuous model must
##                     recover the truth to near machine precision. Model 1 uses BASE with all genetic variance in the
##                     PGS (delta^2 = .5, a = 0).
##   --mode sim        the pooled GeneEvolve trio data written by 01_simulate.py (output/data/); Model 1 is skipped (the
##                     simulated populations have a latent genetic score, so it is misspecified there).
##   --mode reduction  the continuous AMChange scripts fit to the exact-covariance data of the Eq and DisEq coding check
##                     (../output/mvn_data/, from ../02_fit_models.R --mode mvn): the AMChange model nests both, so it
##                     must reach the same -2LL with mu_pp = mu_p (Eq) or mu_pp = 0 (DisEq).
##   --branch B        weaker | stopped | stronger | all (mvn, sim). Default: all for continuous, weaker for binary.
##   --binary 0|1      continuous (0) or binary (1) models only; default both (reduction: continuous only).
##   --exhaustive 0|1  mxTryHard(exhaustive = ...). Default as in ../02_fit_models.R: 1 for continuous fits (the
##                     default search can stop at a local optimum), 0 for binary fits (10-30 min per attempt).
##   --design D, --mech M  restrict to these designs / mating mechanisms (comma-separated; default: all)
##   --extraTries X    mxTryHard extraTries (default 30, the scripts' default). The binary DiffTrait_ObservedYPYM fits
##                     on 48,000 trios take 10-25 min per attempt, and status-6 attempts keep the default search going.
##   --seed S          mvn mode: the seed for the binary data draw (default 2026); other seeds get their own data and
##                     fit files (tag suffix _seedS), for replicate checks of a binary fit
##   --workers K       fits in parallel (mclapply, default 3). Resumable: finished fits in output/results/fits/ are skipped.
##   --threads T       OpenMx threads per fit (default 2)
## Run with: Rscript 02_fit_models.R --mode mvn [--binary 0] [--workers 3]

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
cascadeDir <- dirname(dirname(baseDir))
getArg <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else a[i + 1]
}
mode      <- getArg("--mode", "mvn")
branchSel <- getArg("--branch", NA)
binSel    <- getArg("--binary", NA)
exhSel    <- getArg("--exhaustive", NA)
workers   <- as.integer(getArg("--workers", 3))
designSel <- getArg("--design", NA)
mechSel   <- getArg("--mech", NA)
seed      <- as.integer(getArg("--seed", 2026))
extraTries <- as.integer(getArg("--extraTries", 30))
stopifnot(mode %in% c("mvn", "sim", "reduction"))

suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(MASS); library(parallel) })
# The fitting scripts ask OpenMx for parallel::detectCores() threads; cap what detectCores() returns in this session
threads <- as.integer(getArg("--threads", 2))
assignInNamespace("detectCores", function(...) threads, ns = "parallel")
source(file.path(baseDir, "conditions_AMChange.R"))
amchangeSetup(cascadeDir)
if (mode == "reduction") cascadeValidationSetup(cascadeDir)      # the Eq and DisEq math, for their truth
for (design in AMC_DESIGNS) for (binary in c(FALSE, TRUE)) source(amchangeScriptPath(cascadeDir, design, binary))

fitDir <- file.path(baseDir, "output", "results", "fits"); dir.create(fitDir, recursive = TRUE, showWarnings = FALSE)
mvnDir <- file.path(baseDir, "output", "mvn_data");        dir.create(mvnDir, recursive = TRUE, showWarnings = FALSE)
simDir <- file.path(baseDir, "output", "data")

if (mode == "reduction"){
    tasks <- expand.grid(branch = c("Eq", "DisEq"), mech = names(MECHANISMS), design = DESIGNS, binary = FALSE,
                         stringsAsFactors = FALSE)
} else {
    tasks <- expand.grid(branch = names(BRANCHES), mech = names(MECHANISMS),
                         design = if (mode == "sim") DESIGNS else AMC_DESIGNS,
                         binary = if (is.na(binSel)) c(FALSE, TRUE) else as.logical(as.integer(binSel)),
                         stringsAsFactors = FALSE)
    tasks <- if (is.na(branchSel)) tasks[!tasks$binary | tasks$branch == MAIN_BRANCH, ]
             else if (branchSel == "all") tasks else tasks[tasks$branch == branchSel, ]
}
if (!is.na(designSel)) tasks <- tasks[tasks$design %in% strsplit(designSel, ",")[[1]], ]
if (!is.na(mechSel))   tasks <- tasks[tasks$mech %in% strsplit(mechSel, ",")[[1]], ]
tasks <- tasks[order(tasks$binary, tasks$branch != MAIN_BRANCH), ]      # quick continuous fits first
cat(nrow(tasks), "fits, mode =", mode, ",", workers, "worker(s) x", threads, "thread(s)\n")

mvnData <- function(tag, spec, binary){
    tsv <- file.path(mvnDir, paste0(tag, ".tsv"))
    if (file.exists(tsv)) return(tsv)
    CM <- spec$CM
    stopifnot(max(abs(CM - t(CM))) < 1e-10, !is.null(tryCatch(chol(CM), error = function(e) NULL)))
    set.seed(seed)
    dat <- as.data.frame(mvrnorm(n = if (binary) 5e4 else 32000, mu = rep(0, ncol(CM)), Sigma = CM, empirical = TRUE))
    if (binary) for (v in names(spec$thresholds)) dat[[v]] <- as.integer(dat[[v]] > spec$thresholds[[v]])
    fwrite(dat, tsv, sep = "\t")
    tsv
}

## Reduction mode: the Eq / DisEq truth of ../conditions.R, relabelled for the AMChange scripts, and the original fit
reductionSpec <- function(regime, mech, design){
    tr <- conditionTruth(regime, mech)
    q <- tr$q; o <- tr$o
    diff <- grepl("DiffTrait", design)
    truth <- c(AMGenMulti = q$AM_G, AMEnvMulti = q$AM_E, f11 = q$f, mu_p11 = q$mu,
               mu_pp11 = if (regime == "Eq") q$mu else 0, VY_p11 = q$VY)
    truth <- c(truth, if (!diff) c(VE11 = q$VE, delta11 = q$delta, a11 = q$a)
                      else c(VE_p11 = q$VE, VE_o11 = o$VE_o, delta_p11 = q$delta, a_p11 = q$a,
                             delta_o11 = o$delta_o, a_o11 = o$a_o, VY_o11 = o$VY_o))
    args <- multiplierArgs(mech, grepl("Latent", design))
    if (design == "SameTrait_LatentYPYM") args$h2_RDR <- h2RDR(q$a, q$delta, q$VE, q$VY)
    if (design == "DiffTrait_LatentYPYM"){
        args$h2_RDR_parent    <- h2RDR(q$a, q$delta, q$VE, q$VY)
        args$h2_RDR_offspring <- h2RDR(o$a_o, o$delta_o, o$VE_o, q$VY)
        args$mu_p_fixed       <- q$mu
    }
    if (design == "DiffTrait_ObservedYPYM_EstimatedAParent") args$h2_RDR_offspring <- h2RDR(o$a_o, o$delta_o, o$VE_o, o$VY_o)
    if (design == "DiffTrait_ObservedYPYM_FixedAParent"){
        args$h2_RDR_parent    <- h2RDR(q$a, q$delta, q$VE, q$VY)
        args$h2_RDR_offspring <- h2RDR(o$a_o, o$delta_o, o$VE_o, o$VY_o)
    }
    orig <- file.path(dirname(baseDir), "output", "results", "fits", sprintf("mvn_%s_%s_%s_continuous.csv", regime, mech, design))
    list(truth = truth, args = args,
         data_path = file.path(dirname(baseDir), "output", "mvn_data", sprintf("%s_%s_%s_continuous.tsv", regime, mech, design)),
         m2LL_original = if (file.exists(orig)) fread(orig)$minus2LL[1] else NA)
}

fitOne <- function(i){
    t <- tasks[i, ]
    tag <- sprintf("%s_%s_%s_%s_%s%s", mode, t$branch, t$mech, t$design, if (t$binary) "binary" else "continuous",
                   if (mode == "mvn" && seed != 2026) paste0("_seed", seed) else "")
    out <- file.path(fitDir, paste0(tag, ".csv"))
    if (file.exists(out)) return(invisible(out))
    m2LL_original <- NA
    if (mode == "reduction"){
        rs <- reductionSpec(t$branch, t$mech, t$design)
        truth <- rs$truth; args <- rs$args; data_path <- rs$data_path; m2LL_original <- rs$m2LL_original
        if (!file.exists(data_path)) { cat("missing", data_path, "(run ../02_fit_models.R --mode mvn)\n"); return(invisible(NULL)) }
    } else {
        tr   <- amchangeTruth(t$mech, t$branch, base = if (t$design == "Model1") BASE_MODEL1 else BASE)
        spec <- designSpecAMChange(tr, t$design, t$binary)
        truth <- spec$truth; args <- spec$args
        if (mode == "mvn"){
            data_path <- mvnData(sub("^mvn_", "", tag), spec, t$binary)
            truth <- c(truth, setNames(rep(0, 7), c("meanYp1","meanYm1","meanYo1","meanTp1","meanNTp1","meanTm1","meanNTm1")))
        } else {
            data_path <- file.path(simDir, sprintf("%s_%s_%s%s.tsv", t$branch, t$mech, spec$dataset, if (t$binary) "_binary" else ""))
        }
    }
    exhaustive <- if (!is.na(exhSel)) as.logical(as.integer(exhSel)) else !t$binary
    t0 <- Sys.time()
    res <- tryCatch({
        s <- do.call(amchangeFitName(t$design, t$binary), c(list(data_path = data_path), args,
                                                            list(extraTries = extraTries, exhaustive = exhaustive)))
        p <- s$parameters
        data.frame(parameter = p$name, estimate = p$Estimate, SE = p$Std.Error, true = unname(truth[p$name]),
                   statusCode = as.character(s$statusCode)[1], minus2LL = s$Minus2LogLikelihood, n = s$numObs,
                   m2LL_original = m2LL_original, error = "", stringsAsFactors = FALSE)
    }, error = function(e) data.frame(parameter = NA, estimate = NA, SE = NA, true = NA, statusCode = "ERROR",
                                      minus2LL = NA, n = NA, m2LL_original = m2LL_original, error = conditionMessage(e),
                                      stringsAsFactors = FALSE))
    res <- cbind(mode = mode, seed = if (mode == "mvn") seed else NA, extraTries = extraTries, branch = t$branch, mechanism = t$mech, design = t$design,
                 trait = if (t$binary) "binary" else "continuous", exhaustive = exhaustive, res,
                 time_sec = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
    fwrite(res, out)
    struct <- !is.na(res$true) & !grepl("^mean", res$parameter)
    cat(sprintf("%-80s status=%-6s max|est-true|=%.2e%s  %5.0fs\n", tag, res$statusCode[1],
                if (any(struct)) max(abs(res$estimate[struct] - res$true[struct])) else NA,
                if (mode == "reduction") sprintf("  -2LL - original = %.2e", res$minus2LL[1] - m2LL_original) else "",
                res$time_sec[1]))
    invisible(out)
}

if (workers > 1) invisible(mclapply(seq_len(nrow(tasks)), fitOne, mc.cores = workers, mc.preschedule = FALSE)) else
    invisible(lapply(seq_len(nrow(tasks)), fitOne))

all <- rbindlist(lapply(list.files(fitDir, paste0("^", mode, "_.*\\.csv$"), full.names = TRUE), fread), fill = TRUE)
fwrite(all, file.path(baseDir, "output", "results", paste0("estimates_", mode, ".csv")))
cat("\nWrote", nrow(all), "rows from", uniqueN(all[, paste(branch, mechanism, design, trait)]), "fits\n")
print(all[, .(fits = uniqueN(paste(branch, mechanism, design))), by = .(trait, statusCode)])
