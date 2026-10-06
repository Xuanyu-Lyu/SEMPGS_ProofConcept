## Does per-person covariate control work in the binary scripts? Prototype:
## UniSEMPGS_Binary_SameTrait_ObservedYPYM_PerPerson.R (this folder) vs the current
## OpenMxScripts/EquilibriumBinary/UniSEMPGS_Binary_SameTrait_ObservedYPYM.R.
##
## Data: trios from the binary equilibrium SameTrait model (the values of test_BinaryEquilibriumModels.R), with
## person-specific covariates that act on that person only:
##   age_p, age_m, age_o  each person's age (N(0,1); spouses' ages correlated .5, the child's independent) shifts that
##                        person's liability by BETA_AGE, so the true threshold effect is g = -BETA_AGE
##   PC1_p, PC1_m         each parent's ancestry PC (spouses correlated .3) shifts the means of that parent's two
##                        haplotypic PGS by B_PC
## Fits (every fit starts at the true structural values with covariate effects at 0, unless noted):
##   per_person       prototype, covars = list(Yp1 = "age_p", Ym1 = "age_m", Yo1 = "age_o", Tp1/NTp1 = "PC1_p",
##                    Tm1/NTm1 = "PC1_m")                                                   (correct, 7 effects)
##   trio_all         current script, covars = c("age_p", "age_m", "age_o", "PC1_p", "PC1_m")
##                    (every column on every variable: correct but 35 effects, 28 of them truly 0)
##   trio_child_age   current script, covars = "age_o" (one column for the whole trio)       (misspecified)
##   none             current script, no covariates                                          (misspecified)
##   per_person_default  prototype as a user would run it: the script's own starting values and mxTryHard
## Plus two checks without optimizing: at the true values, per_person and trio_all must give the same -2LL (their
## extra effects are 0 there), and the prototype given a character vector must give the current script's -2LL.
## At most 3 R sessions. Run with: Rscript test_per_person_covariates.R [--n 10000] [--workers 3]

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
getArg  <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else as.numeric(a[i + 1])
}
N       <- getArg("--n", 10000)
workers <- getArg("--workers", 3)
only    <- { a <- commandArgs(trailingOnly = TRUE); i <- match("--only", a); if (is.na(i)) NULL else a[i + 1] }

Sys.setenv(OMP_NUM_THREADS = 1)
suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(MASS); library(parallel) })
assignInNamespace("detectCores", function(...) 1L, ns = "parallel")   # one OpenMx thread per fit
repo <- dirname(dirname(baseDir))
source(file.path(repo, "OpenMxScripts", "EquilibriumBinary", "00-UniIterativeMath_Binary.R"))
source(file.path(repo, "OpenMxScripts", "EquilibriumBinary", "UniSEMPGS_Binary_SameTrait_ObservedYPYM.R"))
source(file.path(baseDir, "UniSEMPGS_Binary_SameTrait_ObservedYPYM_PerPerson.R"))
outDir <- file.path(baseDir, "output"); dir.create(outDir, showWarnings = FALSE)

BETA_AGE <- .4; B_PC <- .15; THRESH <- qnorm(.75)

## ---------------- truth and data ----------------
k <- .5
eq <- uniIterativeMath_Binary(delta = .2, a = sqrt(.6), f = .15, VE = .36, am = .4, k = k, j = .5)
checkEquilibriumConsistency(eq)
cols7 <- c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1")
CM <- with(eq, rbind(
    c(VY,      Yp_Ym,   Yo_Yp,   Omega,   Omega,   Yp_PGSm, Yp_PGSm),
    c(Yp_Ym,   VY,      Yo_Yp,   Yp_PGSm, Yp_PGSm, Omega,   Omega),
    c(Yo_Yp,   Yo_Yp,   VY,      thetaT,  thetaNT, thetaT,  thetaNT),
    c(Omega,   Yp_PGSm, thetaT,  k + gc,  gc,      gt,      gt),
    c(Omega,   Yp_PGSm, thetaNT, gc,      k + gc,  gt,      gt),
    c(Yp_PGSm, Omega,   thetaT,  gt,      gt,      k + gc,  gc),
    c(Yp_PGSm, Omega,   thetaNT, gt,      gt,      gc,      k + gc)))
dimnames(CM) <- list(cols7, cols7)

set.seed(2026)
L   <- as.data.frame(mvrnorm(N, rep(0, 7), CM)); colnames(L) <- cols7
age <- mvrnorm(N, c(0, 0), matrix(c(1, .5, .5, 1), 2)); age_o <- rnorm(N)
pc  <- mvrnorm(N, c(0, 0), matrix(c(1, .3, .3, 1), 2))
dat <- data.table(
    Yp1 = as.integer(L$Yp1 + BETA_AGE * age[, 1] > THRESH),
    Ym1 = as.integer(L$Ym1 + BETA_AGE * age[, 2] > THRESH),
    Yo1 = as.integer(L$Yo1 + BETA_AGE * age_o    > THRESH),
    Tp1 = L$Tp1 + B_PC * pc[, 1], NTp1 = L$NTp1 + B_PC * pc[, 1],
    Tm1 = L$Tm1 + B_PC * pc[, 2], NTm1 = L$NTm1 + B_PC * pc[, 2],
    age_p = age[, 1], age_m = age[, 2], age_o = age_o, PC1_p = pc[, 1], PC1_m = pc[, 2])
tsv <- file.path(outDir, sprintf("trios_n%d.tsv", N)); fwrite(dat, tsv, sep = "\t")
cat("n =", N, "trios; prevalence Yp/Ym/Yo:", round(colMeans(dat[, .(Yp1, Ym1, Yo1)]), 3), "\n")

truth <- with(eq, c(VE11 = VE, delta11 = delta, a11 = a, Omega11 = Omega, Gamma11 = Gamma, mu11 = mu, f11 = f,
                    gt11 = gt, ht11 = ht, gc11 = gc, hc11 = hc, ic11 = ic, w11 = w, v11 = v,
                    thresh_Yp11 = THRESH, thresh_Ym11 = THRESH, thresh_Yo11 = THRESH,
                    g_Yp1_age_p = -BETA_AGE, g_Ym1_age_m = -BETA_AGE, g_Yo1_age_o = -BETA_AGE,
                    b_Tp1_PC1_p = B_PC, b_NTp1_PC1_p = B_PC, b_Tm1_PC1_m = B_PC, b_NTm1_PC1_m = B_PC))

## ---------------- fits ----------------
## The scripts' functions with their final mxTryHard/summary lines replaced, so they return the unfitted model and the
## harness controls the starting values (as in SensitivityAnalyses/FTestValidity/conditions.R).
unfitted <- function(path){
    src <- paste(readLines(path), collapse = "\n")
    src <- sub("fitModel1 <- mxTryHard\\(Model1,[^\n]*", "fitModel1 <- Model1", src)
    src <- sub("return(summary(fitModel1, verbose = TRUE))", "return(fitModel1)", src, fixed = TRUE)
    e <- new.env(); eval(parse(text = src), envir = e); get(ls(e)[1], envir = e)
}
current   <- unfitted(file.path(repo, "OpenMxScripts", "EquilibriumBinary", "UniSEMPGS_Binary_SameTrait_ObservedYPYM.R"))
prototype <- unfitted(file.path(baseDir, "UniSEMPGS_Binary_SameTrait_ObservedYPYM_PerPerson.R"))
allCols <- c("age_p", "age_m", "age_o", "PC1_p", "PC1_m")
perPerson <- list(Yp1 = "age_p", Ym1 = "age_m", Yo1 = "age_o", Tp1 = "PC1_p", NTp1 = "PC1_p", Tm1 = "PC1_m", NTm1 = "PC1_m")
atTruth <- function(m){                 # structural parameters at the truth, covariate effects at 0
    p <- omxGetParameters(m); p[] <- 0
    lab <- intersect(setdiff(names(truth), grep("^[bg]_", names(truth), value = TRUE)), names(p)); p[lab] <- truth[lab]
    omxSetParameters(m, labels = names(p), values = p)
}
noSE <- function(m) mxOption(mxOption(m, "Calculate Hessian", "No"), "Standard Errors", "No")

## checks at fixed parameter values (no optimizer)
m2LLat <- function(m) mxRun(m, useOptimizer = FALSE, silent = TRUE)$output$Minus2LogLikelihood
withTrue <- function(m){ p <- omxGetParameters(m); lab <- intersect(names(truth), names(p)); omxSetParameters(atTruth(m), labels = lab, values = truth[lab]) }
chk <- c(per_person   = m2LLat(withTrue(prototype(tsv, covars = perPerson))),
         trio_all     = m2LLat(withTrue(current(tsv, covars = allCols))),
         proto_vector = m2LLat(withTrue(prototype(tsv, covars = allCols))))
cat(sprintf("\n-2LL at the true values: per_person %.4f | trio_all %.4f | prototype with a vector %.4f\n", chk[1], chk[2], chk[3]))

FITS <- list(
    per_person         = function() withSE(mxRun(atTruth(prototype(tsv, covars = perPerson)), silent = TRUE)),
    per_person_default = function() suppressMessages(mxTryHard(noSE(prototype(tsv, covars = perPerson)), extraTries = 5,
                             OKstatuscodes = c(0, 1), silent = TRUE, jitterDistrib = "rnorm", loc = .2, scale = .05)),
    trio_all           = function() mxRun(noSE(atTruth(current(tsv, covars = allCols))), silent = TRUE),
    trio_child_age     = function() mxRun(noSE(atTruth(current(tsv, covars = "age_o"))), silent = TRUE),
    none               = function() mxRun(noSE(atTruth(current(tsv, covars = NULL))), silent = TRUE))
withSE <- function(fit) fit       # per_person keeps the scripts' Hessian/SE options
if (!is.null(only)) FITS <- FITS[strsplit(only, ",")[[1]]]

runFit <- function(name){
    t0 <- Sys.time()
    fit <- tryCatch(suppressWarnings(FITS[[name]]()), error = function(e) { message(name, ": ", conditionMessage(e)); NULL })
    if (is.null(fit)) return(NULL)
    s <- summary(fit)
    p <- s$parameters
    data.table(fit = name, parameter = p$name, estimate = p$Estimate,
               SE = if ("Std.Error" %in% names(p)) p$Std.Error else NA_real_, true = unname(truth[p$name]),
               m2LL = fit$output$Minus2LogLikelihood, n_free = nrow(p), status = fit$output$status$code,
               minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
}
res <- rbindlist(mclapply(names(FITS), runFit, mc.cores = min(workers, 3), mc.preschedule = FALSE), fill = TRUE)
fwrite(res, file.path(outDir, sprintf("estimates_n%d%s.csv", N, if (is.null(only)) "" else paste0("_", gsub(",", "-", only)))))

## ---------------- summary ----------------
cat("\n==== fits ====\n")
print(unique(res[, .(fit, m2LL = round(m2LL, 3), n_free, status, minutes)]))
key <- c("VE11", "delta11", "a11", "f11", "mu11", "Omega11", "w11", "thresh_Yp11", "thresh_Ym11", "thresh_Yo11")
cat("\n==== structural parameters: estimate (truth in the first column) ====\n")
w <- dcast(res[parameter %in% key], parameter ~ fit, value.var = "estimate")
w <- merge(data.table(parameter = key, true = truth[key]), w, by = "parameter", sort = FALSE)
print(w, digits = 3)
cat("\n==== covariate effects ====\n")
cv <- res[grepl("^[bg]_", parameter)]
cv[is.na(true), true := 0]
print(cv[, .(n = .N, max_abs_error = round(max(abs(estimate - true)), 3), max_abs_z = round(max(abs((estimate - true) / SE)), 2)),
         by = .(fit, kind = ifelse(true == 0, "truly 0 (cross-person)", "own covariate"))])
if (all(c("trio_all", "per_person") %in% res$fit)){
    d <- unique(res[fit %in% c("per_person", "trio_all"), .(fit, m2LL, n_free)])
    lrt <- d[fit == "per_person", m2LL] - d[fit == "trio_all", m2LL]; df <- d[fit == "trio_all", n_free] - d[fit == "per_person", n_free]
    cat(sprintf("\nper_person vs trio_all: LRT = %.2f on %d df (the cross-person effects are truly 0), p = %.3f\n",
                lrt, df, pchisq(lrt, df, lower.tail = FALSE)))
}
if (all(c("per_person", "per_person_default") %in% res$fit))
    cat(sprintf("per_person from the script's default starts: -2LL %.3f vs %.3f from the true values\n",
                unique(res[fit == "per_person_default", m2LL]), unique(res[fit == "per_person", m2LL])))
