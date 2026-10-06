## Model 0 of Balbona, Kim & Keller (2021), as a univariate SEM-PGS model.
## Vertical transmission but no assortative mating; the PGS is assumed to capture all of the trait's heritability (no
## latent genetic score); only the offspring phenotype and the parents' PGS are observed.
## Equilibrium: vertical transmission has gone on for many generations under random mating.
##
## Estimates (returned by summary(fit, verbose = TRUE)):
##   delta       effect of the haplotypic PGS on the phenotype
##   f           vertical transmission from each parent's phenotype to the offspring's environment F
##   VY, VE      phenotypic and residual variance
##   plus Omega, w: covariance terms held to their model values by constraints.
##
## Input:
##   data_path   text file with a header row (read by data.table::fread) and the columns
##                 Yo1         offspring phenotype
##                 Tp1, NTp1   father's transmitted and non-transmitted haplotypic PGS
##                 Tm1, NTm1   mother's transmitted and non-transmitted haplotypic PGS
##               The file should hold only these columns.
##               Each haplotypic PGS is assumed to have variance k = .5 (and no covariance with the others).
##   feaTol, optTol, extraTries, exhaustive, jitterMean, jitterVar  optimizer settings (NPSOL, mxTryHard)

fitUniSEMPGS_Model0 <- function(data_path, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
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
    # Ensure data only contains: Yo1, Tp1, NTp1, Tm1, NTm1

    # Starting values from the sample moments. Model 0 is just-identified, so these put the first attempt near the
    # solution: thetaT - thetaNT = delta*k, Omega = thetaT, w = 2*thetaNT.
    S <- cov(Example_Data[, c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")], use = "pairwise.complete.obs")
    thetaT_s  <- mean(S["Yo1", c("Tp1", "Tm1")])
    thetaNT_s <- mean(S["Yo1", c("NTp1", "NTm1")])
    delta_s   <- (thetaT_s - thetaNT_s) / .5               # k = .5

    # Phenotypic and Residual Variance
    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=S["Yo1", "Yo1"], label="VY11", name="VY", lbound = .001)
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=S["Yo1", "Yo1"] / 2, label="VE11", name="VE", lbound = .001)

    # Scalar Algebra for Variances (paper eqs 4-5)
    VY_Algebra <- mxAlgebra(2 * delta * Omega + w * delta + VF_Algebra + VE, name="VY_Algebra")
    VF_Algebra <- mxAlgebra(2 * f^2 * VY, name="VF_Algebra")
    VY_Constraint <- mxConstraint(VY == VY_Algebra, name='VY_Constraint')

    # Genetic effects
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=delta_s, label="delta11", name="delta")
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    # k is fixed at .5, which assumes the haplotypic PGS has variance 1/2 (the full PGS is standardised). Without
    # assortative mating this is the same in every generation.
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=thetaT_s, label="Omega11", name="Omega")

    # Omega = cov(Y*, [N]T*) = delta*k + w/2 (paper: Omega = (delta + w)/2 with k = 1/2)
    Omega_Algebra <- mxAlgebra(delta * k + 0.5 * w , name="Omega_Algebra")
    Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')

    # Vertical transmission effects (paper eq 3)
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="f11", name="f")
    w     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=2 * thetaNT_s, label="w11", name="w")

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

    Means <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = TRUE, values = 0,
        label = c("meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"),
        dimnames = list(NULL, c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")),
        name = "expMeans")

    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans",
                                            dimnames=c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1"))
    Example_Data_Mx   <- mxData(observed=Example_Data, type="raw")
    FitFunctionML     <- mxFitFunctionML()

    Params <- list(
                VY, VE, delta, k, Omega, f, w,
                VY_Algebra, VF_Algebra, Omega_Algebra, w_Algebra,
                VY_Constraint, Omega_Constraint, w_constraint,
                thetaNT, thetaT,
                CovMatrix, Means, ModelExpectations, FitFunctionML)

    Model1 <- mxModel("UniSEM_Model0", Params, Example_Data_Mx)
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
