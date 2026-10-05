## Tests for the four Model 1 scripts (Balbona, Kim & Keller 2021: VT + AM, PGS explains all heritability, only the
## offspring phenotype and the parental haplotypic PGS observed):
##   Equilibrium/UniSEMPGS_Model1.R                  Disequilibrium/UniSEMPGS_DisEq_Model1.R
##   EquilibriumBinary/UniSEMPGS_Binary_Model1.R     DisequilibriumBinary/UniSEMPGS_DisEq_Binary_Model1.R
##
## A. Parameter recovery (coding + identification), as in the folders' test_*.R: the scalar iterative math is run
##    with a = 0 (the PGS carries all genetic variance), one exact-covariance MVN dataset is drawn per script
##    (mvrnorm empirical = TRUE; the binary scripts get Yo1 dichotomized), and the estimates are compared to the
##    generating values.
## B. The paper's claims about Model 1 when the PGS explains only part of the heritability (a > 0; continuous
##    scripts; data from the SameTrait_LatentYPYM generating process). Without AM, f, VF and w are still unbiased
##    (paper eq 6); with AM the unmodelled PGS-LGS covariance inflates thetaNT and so VF and w (paper, Model 2).
##    delta is unbiased in both cases (thetaT - thetaNT = delta*k per haplotype).
## Run with: Rscript test_Model1.R  (from anywhere)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()

library(MASS)

source(file.path(baseDir, "Equilibrium", "00-UniIterativeMath.R"))
source(file.path(baseDir, "Disequilibrium", "00-UniIterativeMath_DisEq.R"))
source(file.path(baseDir, "EquilibriumBinary", "00-UniIterativeMath_Binary.R"))           # adds uniIterativeMath_Binary
source(file.path(baseDir, "DisequilibriumBinary", "00-UniIterativeMath_DisEq_Binary.R"))  # adds uniIterativeMath_DisEq_Binary
source(file.path(baseDir, "Equilibrium", "UniSEMPGS_Model1.R"))
source(file.path(baseDir, "Disequilibrium", "UniSEMPGS_DisEq_Model1.R"))
source(file.path(baseDir, "EquilibriumBinary", "UniSEMPGS_Binary_Model1.R"))
source(file.path(baseDir, "DisequilibriumBinary", "UniSEMPGS_DisEq_Binary_Model1.R"))

dataDir <- file.path(baseDir, "test_data");    dir.create(dataDir, showWarnings = FALSE)
resDir  <- file.path(baseDir, "test_results"); dir.create(resDir,  showWarnings = FALSE)

n_cont    <- 32000;  extraTries_cont <- 8;   tol_cont <- .02
n_bin     <- 100000; extraTries_bin  <- 15;  tol_bin  <- .05   # ordinal ML is noisier (see ConceptProof.md)
thresh_Yo <- qnorm(0.75)                                         # ~25% prevalence on the parents' liability scale

## ---------------- true generating parameters ----------------
# Same total base-population genetic variance as the other suites (.04 + .60 = .64), but all of it in the PGS
delta <- .8; a <- 0; VE <- .36
f <- .15; am <- .4; k <- .5; j <- .5

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

## Draw one exact-covariance dataset (dichotomizing Yo1 if a threshold is given), write it, fit it.
simulateAndFit <- function(model_name, CMatrix, fit_call, n, thresh = NULL){
    stopifnot(max(abs(CMatrix - t(CMatrix))) < 1e-10)
    if (is.null(tryCatch(chol(CMatrix), error = function(e) NULL)))
        stop("Expected covariance matrix for ", model_name, " is not positive definite")
    set.seed(2026)
    dat <- as.data.frame(MASS::mvrnorm(n = n, mu = rep(0, ncol(CMatrix)), Sigma = CMatrix, empirical = TRUE))
    colnames(dat) <- colnames(CMatrix)
    if (!is.null(thresh)) dat$Yo1 <- as.integer(dat$Yo1 > thresh)
    tsv <- file.path(dataDir, paste0(model_name, ".tsv"))
    write.table(dat, tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    tryCatch(fit_call(tsv), error = function(e){ cat("\nERROR fitting", model_name, ":", conditionMessage(e), "\n"); NULL })
}

compareToTruth <- function(fitSum, truth, model_name, tol_flag){
    pars <- tryCatch(fitSum$parameters, error = function(e) NULL)
    if (!is.data.frame(pars) || nrow(pars) == 0){
        cat("\n========================================\n")
        cat("Model:", model_name, "- FIT FAILED (no parameter estimates returned)\n")
        return(invisible(list(table = NULL, status = "ERROR")))
    }
    cmp <- data.frame(parameter = pars$name,
                      true      = unname(truth[pars$name]),
                      estimate  = pars$Estimate,
                      SE        = pars$Std.Error)
    cmp$diff <- cmp$estimate - cmp$true
    structural <- !grepl("^mean", cmp$parameter)
    flagged <- structural & (is.na(cmp$true) | abs(cmp$diff) > tol_flag)
    cmp$flag <- ifelse(flagged, "FLAG", "")
    cat("\n========================================\n")
    cat("Model:", model_name, "| status code:", paste(as.character(fitSum$statusCode), collapse = " "), "\n")
    print(cmp, digits = 4, row.names = FALSE)
    status <- if (any(flagged)) "FLAG" else "PASS"
    cat(status, "-", model_name, ":", sum(flagged), "structural parameter(s) without truth or off by more than", tol_flag, "\n")
    write.csv(cmp, file.path(resDir, paste0(model_name, "_recovery.csv")), row.names = FALSE)
    invisible(list(table = cmp, status = status))
}

est <- function(fitSum, label) fitSum$parameters$Estimate[fitSum$parameters$name == label]

results <- list()

## ================= A. Parameter recovery (a = 0) =================
eq <- uniIterativeMath(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j)
stopifnot(eq$converged); checkEquilibriumConsistency(eq)
dq <- uniIterativeMath_DisEq(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j)
stopifnot(dq$converged); checkDisEqConsistency(dq)
eb <- uniIterativeMath_Binary(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j)
checkEquilibriumConsistency(eb)
db <- uniIterativeMath_DisEq_Binary(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j)
checkDisEqConsistency(db)
cat("Eq: VY =", eq$VY, "after", eq$generations, "generations | DisEq: parents' VY =", dq$VY, ", offspring VY =", dq$VY_off, "\n")

means5 <- setNames(rep(0, 5), c("meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"))
truthEq    <- function(q) c(VY11 = q$VY, VE11 = q$VE, delta11 = q$delta, Omega11 = q$Omega, mu11 = q$mu,
                            gt11 = q$gt, gc11 = q$gc, f11 = q$f, w11 = q$w)
truthDisEq <- function(q) c(VY11 = q$VY, VE11 = q$VE, delta11 = q$delta, Omega11 = q$Omega, mu11 = q$mu, f11 = q$f)

## 1. Equilibrium, continuous
results$Eq_Model1 <- compareToTruth(
    simulateAndFit("Model1_Eq", build5x5(eq$VY, eq$thetaT, eq$thetaNT, k, eq$gc, eq$gt),
                   function(tsv) fitUniSEMPGS_Model1(tsv, extraTries = extraTries_cont), n_cont),
    c(truthEq(eq), means5), "Model1_Eq", tol_cont)

## 2. Disequilibrium, continuous (Yo1 has the offspring's own variance VY_off)
results$DisEq_Model1 <- compareToTruth(
    simulateAndFit("Model1_DisEq", build5x5(dq$VY_off, dq$thetaT, dq$thetaNT, k, dq$gc, dq$gt),
                   function(tsv) fitUniSEMPGS_DisEq_Model1(tsv, extraTries = extraTries_cont), n_cont),
    c(truthDisEq(dq), means5), "Model1_DisEq", tol_cont)

## 3. Equilibrium, binary (liability scale, VY == 1)
results$Eq_Binary_Model1 <- compareToTruth(
    simulateAndFit("Model1_Eq_Binary", build5x5(eb$VY, eb$thetaT, eb$thetaNT, k, eb$gc, eb$gt),
                   function(tsv) fitUniSEMPGS_Binary_Model1(tsv, extraTries = extraTries_bin, exhaustive = TRUE),
                   n_bin, thresh = thresh_Yo),
    c(truthEq(eb)[names(truthEq(eb)) != "VY11"], thresh_Yo11 = thresh_Yo), "Model1_Eq_Binary", tol_bin)

## 4. Disequilibrium, binary (parents' liability VY == 1; the offspring's liability variance is VY_off)
results$DisEq_Binary_Model1 <- compareToTruth(
    simulateAndFit("Model1_DisEq_Binary", build5x5(db$VY_off, db$thetaT, db$thetaNT, k, db$gc, db$gt),
                   function(tsv) fitUniSEMPGS_DisEq_Binary_Model1(tsv, extraTries = extraTries_bin, exhaustive = TRUE),
                   n_bin, thresh = thresh_Yo),
    c(truthDisEq(db)[names(truthDisEq(db)) != "VY11"], thresh_Yo11 = thresh_Yo), "Model1_DisEq_Binary", tol_bin)

## ================= B. The PGS explains part of the heritability (paper's claims) =================
# The other suites' values: delta^2 = .04 of the .64 base genetic variance is in the PGS
delta_B <- .2; a_B <- sqrt(.60); VE_B <- .36
claims <- list()

## Model 1 quantities implied by the estimates. Eq: VF = 2f^2*VY*(1 + mu*VY), w free. DisEq: the offspring's
## VF_o and w_o (the parents' VF and w carry no AM term), with VY the parents' variance.
impliedEq    <- function(s) with(as.list(setNames(s$parameters$Estimate, s$parameters$name)),
                    c(delta = delta11, f = f11, VF = 2*f11^2*VY11*(1 + mu11*VY11), w = w11))
impliedDisEq <- function(s) with(as.list(setNames(s$parameters$Estimate, s$parameters$name)),
                    c(delta = delta11, f = f11, VF = 2*f11^2*VY11*(1 + mu11*VY11), w = 2*f11*Omega11*(1 + mu11*VY11)))

checkClaim <- function(name, CM, fit_call, implied, truth, expect){
    s <- simulateAndFit(name, CM, fit_call, n_cont)
    cat("\n========================================\n")
    cat("Claim check:", name, "\n")
    if (is.null(s)) return(invisible(list(status = "ERROR")))
    tab <- data.frame(quantity = names(truth), true = unname(truth), model1 = unname(implied(s)[names(truth)]))
    tab$rel_bias <- tab$model1 / tab$true - 1
    print(tab, digits = 4, row.names = FALSE)
    ok <- if (expect == "unbiased") all(abs(tab$model1 - tab$true) < tol_cont)
          else abs(tab$model1[tab$quantity == "delta"] - tab$true[tab$quantity == "delta"]) < tol_cont &&
               all(tab$rel_bias[tab$quantity %in% c("VF", "w")] > .05)
    status <- if (ok) "AS PAPER" else "NOT AS PAPER"
    cat(status, "- expected:", if (expect == "unbiased") "delta, f, VF, w unbiased" else "delta unbiased, VF and w inflated", "\n")
    write.csv(tab, file.path(resDir, paste0(name, "_claim.csv")), row.names = FALSE)
    invisible(list(table = tab, status = status))
}

## 5. Equilibrium, no AM: f, VF and w unbiased whatever the PGS r2
e0 <- uniIterativeMath(delta = delta_B, a = a_B, f = f, VE = VE_B, am = 0, k = k, j = j)
claims$Eq_noAM <- checkClaim("Model1_Eq_PartialPGS_noAM",
    build5x5(e0$VY, e0$thetaT, e0$thetaNT, k, e0$gc, e0$gt),
    function(tsv) fitUniSEMPGS_Model1(tsv, extraTries = extraTries_cont), impliedEq,
    c(delta = delta_B, f = f, VF = e0$VF, w = e0$w), "unbiased")

## 6. Equilibrium, AM: VF and w inflated
e1 <- uniIterativeMath(delta = delta_B, a = a_B, f = f, VE = VE_B, am = am, k = k, j = j)
claims$Eq_AM <- checkClaim("Model1_Eq_PartialPGS_AM",
    build5x5(e1$VY, e1$thetaT, e1$thetaNT, k, e1$gc, e1$gt),
    function(tsv) fitUniSEMPGS_Model1(tsv, extraTries = extraTries_cont), impliedEq,
    c(delta = delta_B, f = f, VF = e1$VF, w = e1$w), "inflated")

## 7. Disequilibrium, AM: VF_o and w_o inflated
d1 <- uniIterativeMath_DisEq(delta = delta_B, a = a_B, f = f, VE = VE_B, am = am, k = k, j = j)
claims$DisEq_AM <- checkClaim("Model1_DisEq_PartialPGS_AM",
    build5x5(d1$VY_off, d1$thetaT, d1$thetaNT, k, d1$gc, d1$gt),
    function(tsv) fitUniSEMPGS_DisEq_Model1(tsv, extraTries = extraTries_cont), impliedDisEq,
    c(delta = delta_B, f = f, VF = d1$VF_o, w = d1$w_o), "inflated")

## ---------------- summary ----------------
cat("\n================ OVERALL ================\n")
cat("A. Recovery with a = 0\n")
for (nm in names(results)) cat(sprintf("   %-28s %s\n", nm, results[[nm]]$status))
cat("B. PGS explains part of the heritability (a > 0)\n")
for (nm in names(claims)) cat(sprintf("   %-28s %s\n", nm, claims[[nm]]$status))
