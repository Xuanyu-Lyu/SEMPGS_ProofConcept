## Tests for the four Model 0 scripts (Balbona, Kim & Keller 2021: vertical transmission, no assortative mating, PGS
## explains all heritability, only the offspring phenotype and the parental haplotypic PGS observed):
##   Equilibrium/UniSEMPGS_Model0.R                  Disequilibrium/UniSEMPGS_DisEq_Model0.R
##   EquilibriumBinary/UniSEMPGS_Binary_Model0.R     DisequilibriumBinary/UniSEMPGS_DisEq_Binary_Model0.R
##
## A. Parameter recovery, as in test_Model1.R: the scalar iterative math is run with am = 0 and a = 0, one
##    exact-covariance MVN dataset is drawn per script (mvrnorm empirical = TRUE; Yo1 dichotomized for the binary
##    scripts), and the estimates are compared to the generating values.
## B. The continuous fit is compared with the reference perform_SEM_model0() (ReferenceCode/VT_SEM_functions.R) on the
##    same data.
## C. The paper's claim for Model 0 (eqs 6-9): without assortative mating, f, VF and w are estimated without bias
##    whatever share of the heritability the PGS captures. Data from the SameTrait_LatentYPYM generating process
##    with a > 0 and am = 0.
## One R session, at most 3 OpenMx threads. Run with: Rscript test_Model0.R  (from anywhere)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()

Sys.setenv(OMP_NUM_THREADS = 3)
suppressPackageStartupMessages({ library(OpenMx); library(MASS); library(stringr); library(parallel) })
assignInNamespace("detectCores", function(...) 3L, ns = "parallel")   # the scripts request detectCores() threads

source(file.path(baseDir, "Equilibrium", "00-UniIterativeMath.R"))
source(file.path(baseDir, "Disequilibrium", "00-UniIterativeMath_DisEq.R"))
source(file.path(baseDir, "EquilibriumBinary", "00-UniIterativeMath_Binary.R"))           # adds uniIterativeMath_Binary
source(file.path(baseDir, "DisequilibriumBinary", "00-UniIterativeMath_DisEq_Binary.R"))  # adds uniIterativeMath_DisEq_Binary
source(file.path(baseDir, "Equilibrium", "UniSEMPGS_Model0.R"))
source(file.path(baseDir, "Disequilibrium", "UniSEMPGS_DisEq_Model0.R"))
source(file.path(baseDir, "EquilibriumBinary", "UniSEMPGS_Binary_Model0.R"))
source(file.path(baseDir, "DisequilibriumBinary", "UniSEMPGS_DisEq_Binary_Model0.R"))
source(file.path(dirname(baseDir), "ReferenceCode", "VT_SEM_functions.R"))

dataDir <- file.path(baseDir, "test_data");    dir.create(dataDir, showWarnings = FALSE)
resDir  <- file.path(baseDir, "test_results"); dir.create(resDir,  showWarnings = FALSE)

n_cont <- 32000; tol_cont <- .02
n_bin  <- 50000; tol_bin  <- .05
thresh_Yo <- qnorm(0.75)                                         # ~25% prevalence on the liability scale

## ---------------- true generating parameters ----------------
# all base-population genetic variance (.64) in the PGS, no assortative mating
delta <- .8; a <- 0; VE <- .36
f <- .15; am <- 0; k <- .5; j <- .5

## ---------------- helpers ----------------
cols5 <- c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")
build5x5 <- function(VYo, thetaT, thetaNT, k, gc, gt){
    M <- rbind(
        c(VYo,     thetaT, thetaNT, thetaT, thetaNT),
        c(thetaT,  k+gc,   gc,      gt,     gt),
        c(thetaNT, gc,     k+gc,    gt,     gt),
        c(thetaT,  gt,     gt,      k+gc,   gc),
        c(thetaNT, gt,     gt,      gc,     k+gc))
    dimnames(M) <- list(cols5, cols5)
    M
}
simulate <- function(name, CMatrix, n, thresh = NULL){
    stopifnot(max(abs(CMatrix - t(CMatrix))) < 1e-10, !is.null(tryCatch(chol(CMatrix), error = function(e) NULL)))
    set.seed(2026)
    dat <- as.data.frame(MASS::mvrnorm(n = n, mu = rep(0, 5), Sigma = CMatrix, empirical = TRUE))
    colnames(dat) <- colnames(CMatrix)
    if (!is.null(thresh)) dat$Yo1 <- as.integer(dat$Yo1 > thresh)
    tsv <- file.path(dataDir, paste0(name, ".tsv"))
    write.table(dat, tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    tsv
}
compareToTruth <- function(fitSum, truth, name, tol_flag){
    pars <- fitSum$parameters
    cmp <- data.frame(parameter = pars$name, true = unname(truth[pars$name]), estimate = pars$Estimate, SE = pars$Std.Error)
    cmp$diff <- cmp$estimate - cmp$true
    structural <- !grepl("^mean", cmp$parameter)
    flagged <- structural & (is.na(cmp$true) | abs(cmp$diff) > tol_flag)
    cmp$flag <- ifelse(flagged, "FLAG", "")
    cat("\n========================================\n")
    cat("Model:", name, "| status code:", paste(as.character(fitSum$statusCode), collapse = " "), "\n")
    print(cmp, digits = 4, row.names = FALSE)
    status <- if (any(flagged)) "FLAG" else "PASS"
    cat(status, "-", name, ":", sum(flagged), "structural parameter(s) without truth or off by more than", tol_flag, "\n")
    write.csv(cmp, file.path(resDir, paste0(name, "_recovery.csv")), row.names = FALSE)
    invisible(list(table = cmp, status = status))
}
results <- list()

## ================= A. Parameter recovery =================
eq <- uniIterativeMath(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j)
stopifnot(eq$converged); checkEquilibriumConsistency(eq)
dq <- uniIterativeMath_DisEq(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j)
stopifnot(dq$converged); checkDisEqConsistency(dq)
stopifnot(abs(dq$VY_off - dq$VY) < 1e-10, abs(eq$VY - dq$VY) < 1e-8)   # no AM: the two regimes coincide
eb <- uniIterativeMath_Binary(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j)
db <- uniIterativeMath_DisEq_Binary(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j)
cat("VY =", eq$VY, "| Omega =", eq$Omega, "| w =", eq$w, "| VF =", eq$VF, "\n")

means5 <- setNames(rep(0, 5), c("meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"))
truthOf <- function(q, binary) {
    t <- c(VY11 = q$VY, VE11 = q$VE, delta11 = q$delta, Omega11 = q$Omega, f11 = q$f, w11 = q$w)
    if (binary) c(t[names(t) != "VY11"], thresh_Yo11 = thresh_Yo) else c(t, means5)
}
tsvEq <- simulate("Model0_Eq", build5x5(eq$VY, eq$thetaT, eq$thetaNT, k, 0, 0), n_cont)
results$Eq_Model0 <- compareToTruth(fitUniSEMPGS_Model0(tsvEq, extraTries = 8), truthOf(eq, FALSE), "Model0_Eq", tol_cont)
tsvDisEq <- simulate("Model0_DisEq", build5x5(dq$VY_off, dq$thetaT, dq$thetaNT, k, 0, 0), n_cont)
results$DisEq_Model0 <- compareToTruth(fitUniSEMPGS_DisEq_Model0(tsvDisEq, extraTries = 8), truthOf(dq, FALSE), "Model0_DisEq", tol_cont)
tsv <- simulate("Model0_Eq_Binary", build5x5(eb$VY, eb$thetaT, eb$thetaNT, k, 0, 0), n_bin, thresh_Yo)
results$Eq_Binary_Model0 <- compareToTruth(fitUniSEMPGS_Binary_Model0(tsv, extraTries = 5, exhaustive = TRUE),
                                           truthOf(eb, TRUE), "Model0_Eq_Binary", tol_bin)
tsv <- simulate("Model0_DisEq_Binary", build5x5(db$VY_off, db$thetaT, db$thetaNT, k, 0, 0), n_bin, thresh_Yo)
results$DisEq_Binary_Model0 <- compareToTruth(fitUniSEMPGS_DisEq_Binary_Model0(tsv, extraTries = 5, exhaustive = TRUE),
                                              truthOf(db, TRUE), "Model0_DisEq_Binary", tol_bin)

## ================= B. Reference implementation =================
cat("\n================ B. perform_SEM_model0() on the same continuous data ================\n")
mine <- fitUniSEMPGS_Model0(tsvEq, extraTries = 8)
ref  <- perform_SEM_model0(read.delim(tsvEq), "NTm1", "Tm1", "NTp1", "Tp1", "Yo1")$result_summary
est  <- setNames(mine$parameters$Estimate, mine$parameters$name)
refCheck <- data.frame(quantity = c("-2LL", "delta", "f", "VE", "VF", "w"),
                       reference = c(ref$ll, ref$deltaest, ref$f, ref$VE, ref$VF, ref$west),
                       new_script = c(mine$Minus2LogLikelihood, est[["delta11"]], est[["f11"]], est[["VE11"]],
                                      2 * est[["f11"]]^2 * est[["VY11"]], est[["w11"]]))
print(refCheck, digits = 7, row.names = FALSE)
results$Reference <- list(status = if (abs(refCheck$reference[1] - refCheck$new_script[1]) < .01 &&
                                       max(abs(refCheck$reference[-1] - refCheck$new_script[-1])) < 1e-3) "PASS" else "FLAG")

## ================= C. PGS explains part of the heritability, no assortative mating =================
cat("\n================ C. a > 0, no assortative mating: f, VF and w should be unbiased ================\n")
e0  <- uniIterativeMath(delta = .2, a = sqrt(.60), f = f, VE = .36, am = 0, k = k, j = j)
tsv <- simulate("Model0_Eq_PartialPGS", build5x5(e0$VY, e0$thetaT, e0$thetaNT, k, 0, 0), n_cont)
s   <- fitUniSEMPGS_Model0(tsv, extraTries = 8)
est <- setNames(s$parameters$Estimate, s$parameters$name)
claim <- data.frame(quantity = c("delta", "f", "VF", "w"), true = c(.2, f, e0$VF, e0$w),
                    model0 = c(est[["delta11"]], est[["f11"]], 2 * est[["f11"]]^2 * est[["VY11"]], est[["w11"]]))
print(claim, digits = 4, row.names = FALSE)
results$PartialPGS_noAM <- list(status = if (all(abs(claim$model0 - claim$true) < tol_cont)) "AS PAPER" else "NOT AS PAPER")

## ---------------- summary ----------------
cat("\n================ OVERALL ================\n")
for (nm in names(results)) cat(sprintf("   %-22s %s\n", nm, results[[nm]]$status))
