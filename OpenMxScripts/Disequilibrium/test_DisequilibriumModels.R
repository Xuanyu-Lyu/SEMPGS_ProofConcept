## Parameter-recovery tests for the 5 disequilibrium univariate SEM-PGS fitting scripts.
## Builds the model-implied covariance matrix from scalar iterative math (VT iterated to its
## random-mating equilibrium, then ONE generation of assortative mating), simulates one
## exact-covariance MVN dataset per model (mvrnorm empirical=TRUE), fits each OpenMx script,
## and compares the estimates to the generating parameters.
## Run with: Rscript test_DisequilibriumModels.R  (from anywhere)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()

library(MASS)

source(file.path(baseDir, "00-UniIterativeMath_DisEq.R"))
source(file.path(baseDir, "UniSEMPGS_DisEq_SameTrait_ObservedYPYM.R"))
source(file.path(baseDir, "UniSEMPGS_DisEq_SameTrait_LatentYPYM.R"))
source(file.path(baseDir, "UniSEMPGS_DisEq_DiffTrait_LatentYPYM.R"))
source(file.path(baseDir, "UniSEMPGS_DisEq_DiffTrait_ObservedYPYM_EstimatedAParent.R"))
source(file.path(baseDir, "UniSEMPGS_DisEq_DiffTrait_ObservedYPYM_FixedAParent.R"))

dataDir <- file.path(baseDir, "test_data");    dir.create(dataDir, showWarnings = FALSE)
resDir  <- file.path(baseDir, "test_results"); dir.create(resDir,  showWarnings = FALSE)

n_samples  <- 32000
extraTries <- 8

## ---------------- true generating parameters ----------------
delta_p <- .2; a_p <- sqrt(.60); VE_p <- .36
f <- .15; am <- .4; k <- .5; j <- .5
delta_o <- .3; a_o <- .65; VE_o_obs <- .5

dq <- uniIterativeMath_DisEq(delta = delta_p, a = a_p, f = f, VE = VE_p, am = am, k = k, j = j)
stopifnot(dq$converged)
checkDisEqConsistency(dq)
cat("VT equilibrium reached after", dq$generations, "generations; VY =", dq$VY,
    "| one-generation AM applied (mu =", dq$mu, ")\n")

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
compareToTruth <- function(fitSum, truth, model_name, tol_flag = .02){
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

runModelTest <- function(model_name, CMatrix, fit_call, truth){
    stopifnot(max(abs(CMatrix - t(CMatrix))) < 1e-10)
    ch <- tryCatch(chol(CMatrix), error = function(e) NULL)
    if (is.null(ch)) stop("Expected covariance matrix for ", model_name, " is not positive definite")
    set.seed(2026)
    dat <- MASS::mvrnorm(n = n_samples, mu = rep(0, ncol(CMatrix)), Sigma = CMatrix, empirical = TRUE)
    colnames(dat) <- colnames(CMatrix)
    tsv <- file.path(dataDir, paste0(model_name, ".tsv"))
    write.table(dat, tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    fitSum <- tryCatch(fit_call(tsv),
                       error = function(e){ cat("\nERROR fitting", model_name, ":", conditionMessage(e), "\n"); NULL })
    if (is.null(fitSum)) return(invisible(list(table = NULL, status = "ERROR")))
    compareToTruth(fitSum, truth, model_name)
}

means7 <- setNames(rep(0, 7), c("meanYp1","meanYm1","meanYo1","meanTp1","meanNTp1","meanTm1","meanNTm1"))
means5 <- setNames(rep(0, 5), c("meanYo1","meanTp1","meanNTp1","meanTm1","meanNTm1"))
results <- list()

## ---------------- 1. SameTrait Observed ----------------
CM_st <- build7x7(dq$VY, dq$VY, dq$Yp_Ym, dq$Yo_Yp, dq$Omega, dq$Yp_PGSm,
                  dq$thetaT, dq$thetaNT, k, dq$gc, dq$gt)
truth_st <- c(VY11 = dq$VY, VE11 = VE_p, delta11 = delta_p, a11 = a_p,
              Omega11 = dq$Omega, Gamma11 = dq$Gamma, mu11 = dq$mu, f11 = f, means7)
results$SameTrait_Observed <- runModelTest("DisEq_SameTrait_ObservedYPYM", CM_st,
    function(tsv) fitUniSEMPGS_DisEq_SameTrait_ObservedYPYM(tsv, extraTries = extraTries),
    truth_st)

## ---------------- 2. SameTrait Latent (RDR for a) ----------------
h2_RDR_st <- (2*a_p^2*j + 2*delta_p^2*k) * (2*a_p^2*j + 2*delta_p^2*k + VE_p) / dq$VY
CM_stl <- build5x5(dq$VY, dq$thetaT, dq$thetaNT, k, dq$gc, dq$gt)
truth_stl <- c(truth_st[setdiff(names(truth_st), names(means7))], means5)
results$SameTrait_Latent <- runModelTest("DisEq_SameTrait_LatentYPYM", CM_stl,
    function(tsv) fitUniSEMPGS_DisEq_SameTrait_LatentParents(tsv, h2_RDR = h2_RDR_st, extraTries = extraTries),
    truth_stl)

## ---------------- offspring-trait quantities for the DiffTrait conditions ----------------
# disequilibrium diff-trait formulas (same algebra as the DisEq fitting scripts)
VYo_terms <- 2*delta_o^2*k + 2*a_o^2*j + 2*delta_o*dq$w + 2*a_o*dq$v + 2*f^2*dq$VY
thetaNT_o <- a_o*dq$itlo + a_o*dq$itol + 2*delta_o*dq$gt + dq$w
thetaT_o  <- 2*delta_o*k + thetaNT_o
Yo_Yp_o   <- delta_o*dq$Omega + a_o*dq$Gamma + f*dq$VY +
             (delta_o*dq$Omega + a_o*dq$Gamma)*dq$mu*dq$VY
h2_p <- (2*a_p^2*j + 2*delta_p^2*k) * (2*a_p^2*j + 2*delta_p^2*k + VE_p) / dq$VY

## ---------------- 3. DiffTrait Latent (VY_o constrained equal to VY_p) ----------------
VE_o_lat <- dq$VY - VYo_terms
stopifnot(VE_o_lat > 0)
h2_o_lat <- (2*a_o^2*j + 2*delta_o^2*k) * (2*a_o^2*j + 2*delta_o^2*k + VE_o_lat) / dq$VY
CM_dtl <- build5x5(dq$VY, thetaT_o, thetaNT_o, k, dq$gc, dq$gt)
truth_dtl <- c(VY11 = dq$VY, VE_parent = VE_p, VE_offspring = VE_o_lat,
               delta_p11 = delta_p, a_p11 = a_p, delta_o11 = delta_o, a_o11 = a_o,
               Omega_p11 = dq$Omega, Gamma_p11 = dq$Gamma, mu11 = dq$mu, f11 = f, means5)
results$DiffTrait_Latent <- runModelTest("DisEq_DiffTrait_LatentYPYM", CM_dtl,
    function(tsv) fitUniSEMPGS_DisEq_DiffTrait_LatentYPYM(tsv, h2_RDR_parent = h2_p,
                    h2_RDR_offspring = h2_o_lat, mu_fixed = dq$mu, extraTries = extraTries),
    truth_dtl)

## ---------------- 4. DiffTrait Observed, estimated parental a ----------------
VY_o <- VYo_terms + VE_o_obs
h2_o_obs <- (2*a_o^2*j + 2*delta_o^2*k) * (2*a_o^2*j + 2*delta_o^2*k + VE_o_obs) / VY_o
CM_dto <- build7x7(dq$VY, VY_o, dq$Yp_Ym, Yo_Yp_o, dq$Omega, dq$Yp_PGSm,
                   thetaT_o, thetaNT_o, k, dq$gc, dq$gt)
truth_dto_base <- c(VE_parent = VE_p, VE_offspring = VE_o_obs,
                    VE_p11 = VE_p, VE_o11 = VE_o_obs,          # FixedAParent label variants
                    delta_p11 = delta_p, a_p11 = a_p, delta_o11 = delta_o, a_o11 = a_o,
                    Omega_p11 = dq$Omega, Gamma_p11 = dq$Gamma, mu11 = dq$mu, f11 = f, means7)
results$DiffTrait_Obs_EstA <- runModelTest("DisEq_DiffTrait_ObservedYPYM_EstimatedAParent", CM_dto,
    function(tsv) fitUniSEMPGS_DisEq_DiffTrait_ObservedYPYM_EstimatedAParent(tsv,
                    h2_RDR_offspring = h2_o_obs, extraTries = extraTries),
    c(VY_parent = dq$VY, VY_offspring = VY_o, truth_dto_base))

## ---------------- 5. DiffTrait Observed, fixed (RDR-anchored) parental a ----------------
results$DiffTrait_Obs_FixedA <- runModelTest("DisEq_DiffTrait_ObservedYPYM_FixedAParent", CM_dto,
    function(tsv) fitUniSEMPGS_DisEq_DiffTrait_ObservedYPYM_FixedAParent(tsv,
                    h2_RDR_parent = h2_p, h2_RDR_offspring = h2_o_obs, extraTries = extraTries),
    c(VY_p11 = dq$VY, VY_o11 = VY_o, truth_dto_base))

## ---------------- summary ----------------
cat("\n================ OVERALL ================\n")
for (nm in names(results)) cat(sprintf("%-25s %s\n", nm, results[[nm]]$status))
cat("Note: a FLAG can indicate an identification problem of that model variant rather than a coding bug.\n")
