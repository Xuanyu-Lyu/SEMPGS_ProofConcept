## Parameter-recovery tests for the 5 binary/liability-threshold disequilibrium univariate SEM-PGS
## fitting scripts. Builds the model-implied liability covariance matrix from scalar iterative math
## (VT iterated to its random-mating equilibrium, then ONE generation of assortative mating) and
## RESCALED so the liability variance is exactly 1 (uniIterativeMath_DisEq_Binary), simulates one
## exact-covariance MVN liability dataset per model (mvrnorm empirical=TRUE), dichotomizes the
## phenotype columns at a chosen threshold, fits each OpenMx liability-threshold script, and
## compares the estimates to the generating parameters.
##
## Binary/ordinal ML carries much less information per observation than continuous ML, so unlike
## the continuous test suite (which gets EXACT recovery from mvrnorm(empirical=TRUE)), recovery
## here is only approximate and improves with n. We use a large n and a looser flag tolerance to
## reflect this -- see ConceptProof.md.
##
## Run with: Rscript test_BinaryDisequilibriumModels.R  (from anywhere)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()

library(MASS)

source(file.path(baseDir, "00-UniIterativeMath_DisEq_Binary.R"))
source(file.path(baseDir, "UniSEMPGS_DisEq_Binary_SameTrait_ObservedYPYM.R"))
source(file.path(baseDir, "UniSEMPGS_DisEq_Binary_SameTrait_LatentYPYM.R"))
source(file.path(baseDir, "UniSEMPGS_DisEq_Binary_DiffTrait_LatentYPYM.R"))
source(file.path(baseDir, "UniSEMPGS_DisEq_Binary_DiffTrait_ObservedYPYM_EstimatedAParent.R"))
source(file.path(baseDir, "UniSEMPGS_DisEq_Binary_DiffTrait_ObservedYPYM_FixedAParent.R"))

dataDir <- file.path(baseDir, "test_data");    dir.create(dataDir, showWarnings = FALSE)
resDir  <- file.path(baseDir, "test_results"); dir.create(resDir,  showWarnings = FALSE)

n_samples  <- 100000
extraTries <- 15   # ordinal ML surfaces have more local optima than exact-covariance continuous ML;
                    # combined with exhaustive=TRUE below, this reliably finds the global optimum
                    # (verified: extraTries=6 without exhaustive intermittently stuck at a spurious
                    # local optimum for the DiffTrait_ObservedYPYM_EstimatedAParent models).
tol_flag   <- .05   # looser than the continuous suite's .02: ordinal ML is inherently noisier

## ---------------- true generating parameters (liability scale, VY == 1 by construction) ----------------
delta_p <- .2; a_p <- sqrt(.60); VE_p <- .36    # pre-standardization inputs (see 00-UniIterativeMath_DisEq_Binary.R)
f <- .15; am <- .4; k <- .5; j <- .5
delta_o <- .2; a_o <- .3                        # offspring path coefficients (chosen directly on the liability scale)

dq <- uniIterativeMath_DisEq_Binary(delta = delta_p, a = a_p, f = f, VE = VE_p, am = am, k = k, j = j)
stopifnot(dq$converged, abs(dq$VY - 1) < 1e-8)
checkDisEqConsistency(dq)
cat("Standardized VT equilibrium reached after", dq$generations, "generations; VY =", dq$VY,
    "(rescale factor s =", dq$scale, ") | one-generation AM applied (mu =", dq$mu, ")\n")

## ---------------- covariance-matrix builders ----------------
cols7 <- c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1")
cols5 <- c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")

build7x7 <- function(VYp, VYo, Yp_Ym, Yo_Yp, Omega, Yp_PGSm, thetaT, thetaNT, k, gc, gt){
    M <- rbind(
        c(VYp,     Yp_Ym,   Yo_Yp,   Omega,   Omega,   Yp_PGSm, Yp_PGSm),
        c(Yp_Ym,   VYp,     Yo_Yp,   Yp_PGSm, Yp_PGSm, Omega,   Omega),
        c(Yo_Yp,   Yo_Yp,   VYo,     thetaT,  thetaNT, thetaT,  thetaNT),
        c(Omega,   Yp_PGSm, thetaT,  k+gc,    gc,      gt,      gt),
        c(Omega,   Yp_PGSm, thetaNT, gc,      k+gc,    gt,      gt),
        c(Yp_PGSm, Omega,   thetaT,  gt,      gt,      k+gc,    gc),
        c(Yp_PGSm, Omega,   thetaNT, gt,      gt,      gc,      k+gc))
    dimnames(M) <- list(cols7, cols7)
    M
}
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

## ---------------- test machinery ----------------
compareToTruth <- function(fitSum, truth, model_name, tol_flag = .05){
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
    flagged <- structural & !is.na(cmp$true) & abs(cmp$diff) > tol_flag
    cmp$flag <- ifelse(flagged, "FLAG", "")
    cat("\n========================================\n")
    cat("Model:", model_name, "| status code:", paste(as.character(fitSum$statusCode), collapse = " "), "\n")
    print(cmp, digits = 4, row.names = FALSE)
    status <- if (any(flagged)) "FLAG" else "PASS"
    cat(status, "-", model_name, ":", sum(flagged), "structural parameter(s) off by more than", tol_flag, "\n")
    write.csv(cmp, file.path(resDir, paste0(model_name, "_recovery.csv")), row.names = FALSE)
    invisible(list(table = cmp, status = status))
}

## Dichotomize the liability columns of an mvrnorm(empirical=TRUE) draw at fixed thresholds.
dichotomize <- function(dat, thresholds){
    for (v in names(thresholds)) dat[[v]] <- as.integer(dat[[v]] > thresholds[[v]])
    dat
}

runModelTest <- function(model_name, CMatrix, thresholds, fit_call, truth){
    stopifnot(max(abs(CMatrix - t(CMatrix))) < 1e-10)
    ch <- tryCatch(chol(CMatrix), error = function(e) NULL)
    if (is.null(ch)) stop("Expected covariance matrix for ", model_name, " is not positive definite")
    set.seed(2026)
    dat <- MASS::mvrnorm(n = n_samples, mu = rep(0, ncol(CMatrix)), Sigma = CMatrix, empirical = TRUE)
    colnames(dat) <- colnames(CMatrix)
    dat <- as.data.frame(dat)
    dat <- dichotomize(dat, thresholds)
    tsv <- file.path(dataDir, paste0(model_name, ".tsv"))
    write.table(dat, tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    fitSum <- tryCatch(fit_call(tsv),
                       error = function(e){ cat("\nERROR fitting", model_name, ":", conditionMessage(e), "\n"); NULL })
    if (is.null(fitSum)) return(invisible(list(table = NULL, status = "ERROR")))
    compareToTruth(fitSum, truth, model_name, tol_flag = tol_flag)
}

## Thresholds (chosen on the standardized N(0,1) liability scale):
##  - SameTrait models: Yp1/Ym1/Yo1 share ONE threshold (same trait -> same prevalence).
##  - DiffTrait models: Yp1/Ym1 (parent trait) share one threshold; Yo1 (a different trait) gets
##    a distinct one, demonstrating that the shared-vs-separate threshold structure is coded right.
thresh_same    <- qnorm(0.75)   # ~25% prevalence
thresh_parent  <- qnorm(0.75)
thresh_offspr  <- qnorm(0.60)   # ~40% prevalence (deliberately different from the parent trait)

truth_thresh_same    <- c(thresh_Y11 = thresh_same)
truth_thresh_parent  <- c(thresh_Yp11 = thresh_parent, thresh_Yo11 = thresh_offspr)

results <- list()

## ---------------- 1. SameTrait Observed ----------------
CM_st <- build7x7(dq$VY, dq$VY, dq$Yp_Ym, dq$Yo_Yp, dq$Omega, dq$Yp_PGSm,
                  dq$thetaT, dq$thetaNT, k, dq$gc, dq$gt)
truth_st <- c(VE11 = dq$VE, delta11 = dq$delta, a11 = dq$a,
              Omega11 = dq$Omega, Gamma11 = dq$Gamma, mu11 = dq$mu, f11 = f, truth_thresh_same)
results$SameTrait_Observed <- runModelTest("DisEq_Binary_SameTrait_ObservedYPYM", CM_st,
    list(Yp1 = thresh_same, Ym1 = thresh_same, Yo1 = thresh_same),
    function(tsv) fitUniSEMPGS_DisEq_Binary_SameTrait_ObservedYPYM(tsv, extraTries = extraTries, exhaustive = TRUE),
    truth_st)

## ---------------- 2. SameTrait Latent (RDR for a) ----------------
h2_RDR_st <- (2*dq$a^2*j + 2*dq$delta^2*k) * (2*dq$a^2*j + 2*dq$delta^2*k + dq$VE) / dq$VY
CM_stl <- build5x5(dq$VY, dq$thetaT, dq$thetaNT, k, dq$gc, dq$gt)
truth_stl <- c(truth_st[setdiff(names(truth_st), names(truth_thresh_same))], c(thresh_Yo11 = thresh_same))
results$SameTrait_Latent <- runModelTest("DisEq_Binary_SameTrait_LatentYPYM", CM_stl,
    list(Yo1 = thresh_same),
    function(tsv) fitUniSEMPGS_DisEq_Binary_SameTrait_LatentParents(tsv, h2_RDR = h2_RDR_st, extraTries = extraTries, exhaustive = TRUE),
    truth_stl)

## ---------------- offspring-trait quantities for the DiffTrait conditions ----------------
## Solve VE_o so the offspring's implied liability variance is ALSO exactly 1 (see
## 00-UniIterativeMath_DisEq_Binary.R for the closed-form derivation).
off <- solve_VYo_DisEq_Binary(delta_o = delta_o, a_o = a_o, k = k, j = j, parent = dq, target = 1)
h2_p <- (2*dq$a^2*j + 2*dq$delta^2*k) * (2*dq$a^2*j + 2*dq$delta^2*k + dq$VE) / dq$VY
h2_o <- (2*a_o^2*j + 2*delta_o^2*k) * (2*a_o^2*j + 2*delta_o^2*k + off$VE_o) / off$VY_o

## ---------------- 3. DiffTrait Latent (VY_o constrained equal to VY_p == 1) ----------------
CM_dtl <- build5x5(off$VY_o, off$thetaT_o, off$thetaNT_o, k, dq$gc, dq$gt)
truth_dtl <- c(VE_parent = dq$VE, VE_offspring = off$VE_o,
               delta_p11 = dq$delta, a_p11 = dq$a, delta_o11 = delta_o, a_o11 = a_o,
               Omega_p11 = dq$Omega, Gamma_p11 = dq$Gamma, mu11 = dq$mu, f11 = f,
               c(thresh_Yo11 = thresh_offspr))
results$DiffTrait_Latent <- runModelTest("DisEq_Binary_DiffTrait_LatentYPYM", CM_dtl,
    list(Yo1 = thresh_offspr),
    function(tsv) fitUniSEMPGS_DisEq_Binary_DiffTrait_LatentYPYM(tsv, h2_RDR_parent = h2_p,
                    h2_RDR_offspring = h2_o, mu_fixed = dq$mu, extraTries = extraTries, exhaustive = TRUE),
    truth_dtl)

## ---------------- 4. DiffTrait Observed, estimated parental a ----------------
CM_dto <- build7x7(dq$VY, off$VY_o, dq$Yp_Ym, off$Yo_Yp_o, dq$Omega, dq$Yp_PGSm,
                   off$thetaT_o, off$thetaNT_o, k, dq$gc, dq$gt)
truth_dto_base <- c(VE_parent = dq$VE, VE_offspring = off$VE_o,
                    VE_p11 = dq$VE, VE_o11 = off$VE_o,          # FixedAParent label variants
                    delta_p11 = dq$delta, a_p11 = dq$a, delta_o11 = delta_o, a_o11 = a_o,
                    Omega_p11 = dq$Omega, Gamma_p11 = dq$Gamma, mu11 = dq$mu, f11 = f,
                    truth_thresh_parent)
results$DiffTrait_Obs_EstA <- runModelTest("DisEq_Binary_DiffTrait_ObservedYPYM_EstimatedAParent", CM_dto,
    list(Yp1 = thresh_parent, Ym1 = thresh_parent, Yo1 = thresh_offspr),
    function(tsv) fitUniSEMPGS_DisEq_Binary_DiffTrait_ObservedYPYM_EstimatedAParent(tsv,
                    h2_RDR_offspring = h2_o, extraTries = extraTries, exhaustive = TRUE),
    truth_dto_base)

## ---------------- 5. DiffTrait Observed, fixed (RDR-anchored) parental a ----------------
results$DiffTrait_Obs_FixedA <- runModelTest("DisEq_Binary_DiffTrait_ObservedYPYM_FixedAParent", CM_dto,
    list(Yp1 = thresh_parent, Ym1 = thresh_parent, Yo1 = thresh_offspr),
    function(tsv) fitUniSEMPGS_DisEq_Binary_DiffTrait_ObservedYPYM_FixedAParent(tsv,
                    h2_RDR_parent = h2_p, h2_RDR_offspring = h2_o, extraTries = extraTries, exhaustive = TRUE),
    truth_dto_base)

## ---------------- summary ----------------
cat("\n================ OVERALL ================\n")
for (nm in names(results)) cat(sprintf("%-25s %s\n", nm, results[[nm]]$status))
cat("Note: recovery here is approximate (finite-sample binary/ordinal ML), unlike the continuous\n",
    "suite's exact-covariance recovery. A FLAG can indicate an identification problem or just\n",
    "sampling noise at n =", n_samples, "-- check the per-parameter Std.Error before concluding a bug.\n")
