## Model 0 of Balbona, Kim & Keller (2021), as a univariate SEM-PGS model.
## Vertical transmission but no assortative mating; the PGS is assumed to capture all of the trait's heritability (no
## latent genetic score); only the offspring phenotype and the parents' PGS are observed.
## Equilibrium: vertical transmission has gone on for many generations under random mating.
## Binary trait: fit as a liability-threshold model with the liability variance fixed at 1.
##
## Estimates (returned by summary(fit, verbose = TRUE)):
##   delta       effect of the haplotypic PGS on the phenotype
##   f           vertical transmission from each parent's phenotype to the offspring's environment F
##   VE          residual variance
##   thresholds  liability thresholds (thresh_Yo)
##   plus Omega, w: covariance terms held to their model values by constraints.
##
## Input:
##   data_path   text file with a header row (read by data.table::fread) and the columns
##                 Yo1         offspring phenotype (0/1)
##                 Tp1, NTp1   father's transmitted and non-transmitted haplotypic PGS
##                 Tm1, NTm1   mother's transmitted and non-transmitted haplotypic PGS
##               Other columns (e.g. covariates) may be present.
##               Each haplotypic PGS is assumed to have variance k = .5 (and no covariance with the others).
##   covars      optional covariates, each applied only to the person it belongs to: a named list from data variables
##               to that person's covariate columns, e.g. list(Yo1 = c("age_o", "sex_o"), Tp1 = "PC1_p", NTp1 = "PC1_p",
##               Tm1 = "PC1_m", NTm1 = "PC1_m"). A covariate shifts the threshold or PGS mean of only the variables it
##               is listed under; variables not listed get no covariates.
##   feaTol, optTol, extraTries, exhaustive, jitterMean, jitterVar  optimizer settings (NPSOL, mxTryHard)

fitUniSEMPGS_Binary_Model0 <- function(data_path, covars = NULL, feaTol = 1e-6, optTol = 1e-8, jitterMean = .2, jitterVar = .05, extraTries = 30, exhaustive = F){
    library(OpenMx)
    library(data.table)
    library(stringr)

    mxOption(NULL,"Calculate Hessian","Yes")
    mxOption(NULL,"Standard Errors","Yes")
    mxOption(NULL,"Default optimizer","NPSOL")
    mxOption(NULL,"Feasibility tolerance",as.character(feaTol))
    mxOption(NULL,"Optimality tolerance",as.character(optTol))
    mxOption(NULL,"Number of Threads", value = parallel::detectCores())

    Example_Data  <- fread(data_path, header = T)
    obs_vars    <- c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")   # extra columns (e.g. covariates) are allowed
    binary_vars <- "Yo1"
    Example_Data[["Yo1"]] <- mxFactor(Example_Data[["Yo1"]], levels = c(0, 1))
    # covars: a named list from data variables to that person's own covariate columns (see the header); each covariate
    # adjusts only the variables it is listed under
    if (length(covars) > 0 && (!is.list(covars) || is.null(names(covars)) || any(names(covars) == "")))
        stop("covars must be a named list from data variables to that person's covariate columns, e.g. ",
             "list(Yo1 = c(\"age_o\", \"sex_o\"), Tp1 = \"PC1_p\", NTp1 = \"PC1_p\")")
    badVar <- setdiff(names(covars), obs_vars)
    if (length(badVar) > 0) stop("covars names must be data variables (", paste(obs_vars, collapse = ", "), "): ",
                                 paste(badVar, collapse = ", "))
    covCols <- unique(unlist(covars, use.names = FALSE))
    if (length(covCols) > 0) {
        missing_cov <- setdiff(covCols, colnames(Example_Data))
        if (length(missing_cov) > 0) stop("Covariate column(s) not found in data: ", paste(missing_cov, collapse = ", "))
        Example_Data <- Example_Data[complete.cases(Example_Data[, covCols, with = FALSE]), ]
    }

    # Phenotypic and Residual Variance
    # Liability variance is FIXED at 1 (not free): identifies the scale of the binary/threshold model.
    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1, label="VY11", name="VY")
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="VE11", name="VE", lbound = .001)

    # Scalar Algebra for Variances (paper eqs 4-5)
    VY_Algebra <- mxAlgebra(2 * delta * Omega + w * delta + VF_Algebra + VE, name="VY_Algebra")
    VF_Algebra <- mxAlgebra(2 * f^2 * VY, name="VF_Algebra")
    VY_Constraint <- mxConstraint(VY == VY_Algebra, name='VY_Constraint')

    # Genetic effects
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="delta11", name="delta")
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="Omega11", name="Omega")

    # Omega = cov(Y*, [N]T*) = delta*k + w/2
    Omega_Algebra <- mxAlgebra(delta * k + 0.5 * w , name="Omega_Algebra")
    Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')

    # Vertical transmission effects (paper eq 3)
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="f11", name="f")
    w     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.05, label="w11", name="w")

    w_Algebra     <- mxAlgebra(2 * f * Omega, name="w_Algebra")
    w_constraint  <- mxConstraint(w == w_Algebra, name='w_constraint')

    # Between-people covariances (Excluding Parent Phenotypes)
    # thetaNT / thetaT: covariance of Yo with ONE parental haplotype, i.e. with one data column (NTp1, Tp1, ...).
    # In the paper theta is the SUM over the father's and the mother's haplotype (eqs 1-2): twice these.
    thetaNT <- mxAlgebra(.5 * w, name="thetaNT")
    thetaT  <- mxAlgebra(delta * k + thetaNT, name="thetaT")

    # Expected covariances matrix (5x5: Yo1, Tp1, NTp1, Tm1, NTm1); without assortative mating the four haplotypic
    # PGS are uncorrelated
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_Algebra, thetaT,  thetaNT, thetaT,  thetaNT), # Yo1
        cbind(thetaT,     k,       0,       0,       0),       # Tp1
        cbind(thetaNT,    0,       k,       0,       0),       # NTp1
        cbind(thetaT,     0,       0,       k,       0),       # Tm1
        cbind(thetaNT,    0,       0,       0,       k)),      # NTm1
        dimnames = list(c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1"),
                        c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")), name="expCov")

    # Yo1's mean is fixed at 0 (identification requires fixing either the mean or the threshold);
    # the PGS means stay free.
    # Covariates (if any) enter as definition variables, as in r2_omx_partial(): a PGS mean becomes mean + sum(b * covar)
    # and a threshold becomes t0 + sum(g * covar), over the covariates covars lists for that variable; every other
    # covariate-variable effect is fixed at 0.
    Means0 <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = c(F, T, T, T, T), values = 0,
        label = c("meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"),
        dimnames = list(NULL, obs_vars), name = "Means0")
    Th0 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = 0,
        labels = "thresh_Yo11", name = "Th0")
    if (length(covCols) > 0) {
        nCov   <- length(covCols)
        uses   <- matrix(sapply(obs_vars, function(vr) covCols %in% covars[[vr]]), nrow = nCov,
                         dimnames = list(covCols, obs_vars))                 # does covariate i adjust variable v?
        bFree  <- uses & matrix(!(obs_vars %in% binary_vars), nCov, length(obs_vars), byrow = TRUE)
        gFree  <- uses[, binary_vars, drop = FALSE]
        defCov <- mxMatrix(type = "Full", nrow = 1, ncol = nCov, free = FALSE, labels = paste0("data.", covCols), name = "defCov")
        bCov   <- mxMatrix(type = "Full", nrow = nCov, ncol = length(obs_vars), free = bFree, values = 0,
            labels = ifelse(bFree, outer(covCols, obs_vars, function(cv, vr) paste0("b_", vr, "_", cv)), NA), name = "bCov")
        gCov   <- mxMatrix(type = "Full", nrow = nCov, ncol = length(binary_vars), free = gFree, values = 0,
            labels = ifelse(gFree, outer(covCols, binary_vars, function(cv, vr) paste0("g_", vr, "_", cv)), NA), name = "gCov")
        Means   <- mxAlgebra(Means0 + defCov %*% bCov, dimnames = list(NULL, obs_vars), name = "expMeans")
        Th      <- mxAlgebra(Th0 + defCov %*% gCov, name = "Th")
        CovObjs <- list(Means0, Th0, defCov, bCov, gCov)
    } else {
        Means   <- mxAlgebra(Means0, dimnames = list(NULL, obs_vars), name = "expMeans")
        Th      <- mxAlgebra(Th0, name = "Th")
        CovObjs <- list(Means0, Th0)
    }

    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans", dimnames=obs_vars,
                                            thresholds="Th", threshnames=binary_vars)
    Example_Data_Mx   <- mxData(observed=Example_Data, type="raw")
    FitFunctionML     <- mxFitFunctionML(jointConditionOn = "continuous")

    Params <- list(
                VY, VE, delta, k, Omega, f, w,
                VY_Algebra, VF_Algebra, Omega_Algebra, w_Algebra,
                VY_Constraint, Omega_Constraint, w_constraint,
                thetaNT, thetaT,
                CovMatrix, CovObjs, Means, Th, ModelExpectations, FitFunctionML)

    Model1 <- mxModel("UniSEM_Binary_Model0", Params, Example_Data_Mx)
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
