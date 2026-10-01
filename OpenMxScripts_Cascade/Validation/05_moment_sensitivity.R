## Step 5: which data moments move the SameTrait_ObservedYPYM estimates away from the truth?
##
## The 7x7 covariance matrix of (Yp, Ym, Yo, Tp, NTp, Tm, NTm) has only 11 distinct moments under the model
## (VY, cov(Yp,Ym), cov(Yo,Yp), var(Yo), Omega, cov(Yp,PGSm), thetaT, thetaNT, var(T), gc, gt). For each
## regime x mechanism I fit the model to covariance matrices built from the model-implied matrix at the truth:
##   truth        the implied matrix itself (must return the truth)
##   <moment>     the implied matrix with the cells of that one moment replaced by the pooled GeneEvolve values
##   all          the pooled GeneEvolve covariance matrix (should reproduce the pooled raw-data fit)
## so the shift each single moment causes can be read off directly.   -> output/results/moment_sensitivity.csv
##
## Run with: Rscript 05_moment_sensitivity.R [--workers 6]

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
cascadeDir <- dirname(baseDir)
getArg <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else a[i + 1]
}
workers <- as.integer(getArg("--workers", 6))

suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(parallel) })
assignInNamespace("detectCores", function(...) 1L, ns = "parallel")
source(file.path(baseDir, "conditions.R"))
cascadeValidationSetup(cascadeDir)

DESIGN <- "SameTrait_ObservedYPYM"
# cells (upper triangle) of each distinct moment in the 7x7 matrix
MOMENTS <- list(VY = list(c(1,1), c(2,2)), Yp_Ym = list(c(1,2)), Yo_Yp = list(c(1,3), c(2,3)), VY_o = list(c(3,3)),
                Omega = list(c(1,4), c(1,5), c(2,6), c(2,7)), Yp_PGSm = list(c(1,6), c(1,7), c(2,4), c(2,5)),
                thetaT = list(c(3,4), c(3,6)), thetaNT = list(c(3,5), c(3,7)),
                var_T = list(c(4,4), c(5,5), c(6,6), c(7,7)), gc = list(c(4,5), c(6,7)),
                gt = list(c(4,6), c(4,7), c(5,6), c(5,7)))

buildModel <- function(regime, design, args, data_path){
    src <- paste(readLines(scriptPath(cascadeDir, regime, design, FALSE)), collapse = "\n")
    src <- sub("fitModel1 <- mxTryHard\\((?s).*", "return(Model1)\n}\n", src, perl = TRUE)
    env <- new.env(); eval(parse(text = src), envir = env)
    do.call(get(fitFunctionName(regime, design, FALSE), envir = env), c(list(data_path = data_path), args))
}

tasks <- CJ(regime = REGIMES, mech = names(MECHANISMS), which = c("truth", names(MOMENTS), "all"), sorted = FALSE)

fitOne <- function(i){
    t <- tasks[i]
    tr <- conditionTruth(t$regime, t$mech); spec <- designSpec(tr, DESIGN, FALSE)
    data_path <- file.path(baseDir, "output", "data", sprintf("%s_%s_SameTrait_Observed.tsv", t$regime, t$mech))
    d <- fread(data_path); cols <- colnames(spec$CM)
    S <- cov(as.matrix(d[, cols, with = FALSE])); n <- nrow(d)
    C <- spec$CM
    if (t$which == "all") C <- S
    else if (t$which != "truth") for (ij in MOMENTS[[t$which]]){ C[ij[1], ij[2]] <- S[ij[1], ij[2]]; C[ij[2], ij[1]] <- S[ij[2], ij[1]] }
    m <- buildModel(t$regime, DESIGN, spec$args, data_path)
    m <- mxModel(m, mxData(observed = C, type = "cov", numObs = n, means = setNames(rep(0, length(cols)), cols)))
    fl <- intersect(names(omxGetParameters(m)), names(spec$truth))
    m <- omxSetParameters(m, labels = fl, values = unname(spec$truth[fl]))
    f <- tryCatch(suppressMessages(mxTryHard(m, extraTries = 10, exhaustive = TRUE, silent = TRUE, OKstatuscodes = c(0, 1))),
                  error = function(e) NULL)
    if (is.null(f)) return(data.table(regime = t$regime, mechanism = t$mech, which = t$which, parameter = NA, status = "ERROR"))
    p <- omxGetParameters(f)
    data.table(regime = t$regime, mechanism = t$mech, which = t$which, parameter = names(p), estimate = unname(p),
               true = unname(spec$truth[names(p)]), status = as.character(f$output$status$code), n = n)
}
cat(nrow(tasks), "fits\n")
res <- rbindlist(mclapply(seq_len(nrow(tasks)), fitOne, mc.cores = workers, mc.preschedule = FALSE), fill = TRUE)
fwrite(res, file.path(baseDir, "output", "results", "moment_sensitivity.csv"))
cat("done\n")
