## Binary/liability-threshold version of UniSEMPGS_SameTrait_ObservedYPYM.R: parent and offspring
## express the SAME binary trait (liability-threshold model) and the parental phenotypes are
## observed. Yp1/Ym1/Yo1 are dichotomous (0/1); Tp1/NTp1/Tm1/NTm1 (the transmitted / non-transmitted
## PGS) stay continuous. Identification of the liability scale follows the standard
## liability-threshold convention: total liability variance VY is fixed to 1 (not estimated), and
## Yp1/Ym1/Yo1 share a single free threshold (same trait => same population prevalence assumed).

fitUniSEMPGS_Binary_SameTrait_ObservedYPYM <- function(data_path, feaTol = 1e-6, optTol = 1e-8, jitterMean = .2, jitterVar = .05, extraTries = 30, exhaustive = F){
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
    binary_vars <- c("Yp1", "Ym1", "Yo1")
    for (v in binary_vars) Example_Data[[v]] <- mxFactor(Example_Data[[v]], levels = c(0, 1))

    # Phenotypic and Residual Variance
    # Liability variance is FIXED at 1 (not free): identifies the scale of the binary/threshold model.
    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1, label="VY11", name="VY")
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="VE11", name="VE", lbound = .001)

    # Scalar Algebra for Variances
    VY_Algebra <- mxAlgebra(2 * delta * Omega + 2 * a * Gamma + w * delta + v * a + VF_Algebra + VE, name="VY_Algebra")
    VF_Algebra <- mxAlgebra(2 * f^2 * VY + 2 * f^2 * VY^2 * mu, name="VF_Algebra")

    VY_Constraint    <- mxConstraint(VY == VY_Algebra, name='VY_Constraint')

    # Genetic effects
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="delta11", name="delta")
    a     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="a11", name="a", lbound = .001)
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    j     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="Omega11", name="Omega")
    Gamma <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="Gamma11", name="Gamma")

    Omega_Algebra <- mxAlgebra(2 * delta * gc + 2 * a * ic + delta * k + 0.5 * w , name="Omega_Algebra")
    Gamma_Algebra <- mxAlgebra(2 * a * hc + 2 * delta * ic + a * j + 0.5 * v, name="Gamma_Algebra")

    Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
    Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')

    adelta_Constraint_Algebra <- mxAlgebra(delta, name = "adelta_Constraint_Algebra")
    adelta_Constraint <- mxConstraint(a == delta, name = "adelta_Constraint")

    j_Algebra    <- mxAlgebra(k, name = "j_Algebra")
    j_constraint <- mxConstraint(j == j_Algebra, name = "j_constraint")

    # Assortative mating effects
    mu    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1,  label="mu11", name="mu")
    gt    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.02, label="gt11", name="gt")
    ht    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.02, label="ht11", name="ht")
    gc    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.02, label="gc11", name="gc")
    hc    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.02, label="hc11", name="hc")

    gt_Algebra <- mxAlgebra(Omega * mu * Omega, name="gt_Algebra")
    ht_Algebra <- mxAlgebra(Gamma * mu * Gamma, name="ht_Algebra")
    gc_Algebra <- mxAlgebra(gt, name="gc_Algebra") # Simplified for scalar
    hc_Algebra <- mxAlgebra(ht, name="hc_Algebra")
    gchc_constraint_Algebra <- mxAlgebra(hc * ( (2 * delta^2 * k) / (2 * a^2 * j) ), name = "gchc_constraint_Algebra")

    itlo  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.02, label="itlo11", name="itlo")
    itol  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.02, label="itol11", name="itol")
    ic    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.02, label="ic11",   name="ic")

    itlo_Algebra <- mxAlgebra(Gamma * mu * Omega, name="itlo_Algebra")
    itol_Algebra <- mxAlgebra(Omega * mu * Gamma, name="itol_Algebra")
    ic_Algebra   <- mxAlgebra(.5 * (itlo + itol), name="ic_Algebra")

    gt_constraint   <- mxConstraint(gt == gt_Algebra, name='gt_constraint')
    ht_constraint   <- mxConstraint(ht == ht_Algebra, name='ht_constraint')
    gc_constraint   <- mxConstraint(gc == gc_Algebra, name='gc_constraint')
    hc_constraint   <- mxConstraint(hc == hc_Algebra, name='hc_constraint')
    gchc_constraint <- mxConstraint(gc == gchc_constraint_Algebra, name='gchc_constraint')
    itlo_constraint <- mxConstraint(itlo == itlo_Algebra, name='itlo_constraint')
    itol_constraint <- mxConstraint(itol == itol_Algebra, name='itol_constraint')
    ic_constraint   <- mxConstraint(ic == ic_Algebra, name='ic_constraint')

    # Vertical transmission effects
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="f11", name="f")
    w     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.05, label="w11", name="w")
    v     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.05, label="v11", name="v")

    w_Algebra     <- mxAlgebra(2 * f * Omega + 2 * f * VY * mu * Omega, name="w_Algebra")
    v_Algebra     <- mxAlgebra(2 * f * Gamma + 2 * f * VY * mu * Gamma, name="v_Algebra")
    wv_constraint_algebra <- mxAlgebra((w * sqrt(2 * delta^2 * k) / sqrt(2 * a^2 * j)), name='wv_constraint_algebra')

    v_constraint  <- mxConstraint(v == v_Algebra, name='v_constraint')
    w_constraint  <- mxConstraint(w == w_Algebra, name='w_constraint')
    wv_constraint <- mxConstraint(v == wv_constraint_algebra, name='wv_constraint')

    # Between-people covariances
    thetaNT <- mxAlgebra(2 * delta * gc + 2 * a * ic + .5 * w, name="thetaNT")
    thetaT  <- mxAlgebra(delta * k + thetaNT, name="thetaT")
    Yp_PGSm <- mxAlgebra(VY * mu * Omega, name="Yp_PGSm")
    Ym_PGSp <- mxAlgebra(VY * mu * Omega, name="Ym_PGSp")
    Yp_Ym   <- mxAlgebra(VY * mu * VY,    name="Yp_Ym")
    Ym_Yp   <- mxAlgebra(VY * mu * VY,    name="Ym_Yp")
    Yo_Yp   <- mxAlgebra(delta * Omega + a * Gamma + (delta * Omega + a * Gamma) * mu * VY + f * VY + f * VY^2 * mu, name = "Yo_Yp")
    Yo_Ym   <- mxAlgebra(delta * Omega + a * Gamma + (delta * Omega + a * Gamma) * mu * VY + f * VY + f * VY^2 * mu, name = "Yo_Ym")

    # Expected covariances matrix
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_Algebra, Yp_Ym,      Yo_Yp,     Omega,   Omega,   Yp_PGSm, Yp_PGSm),
        cbind(Ym_Yp,      VY_Algebra, Yo_Ym,     Ym_PGSp, Ym_PGSp, Omega,   Omega),
        cbind(Yo_Yp,      Yo_Ym,      VY_Algebra,thetaT,  thetaNT, thetaT,  thetaNT),
        cbind(Omega,      Ym_PGSp,    thetaT,    k+gc,    gc,      gt,      gt),
        cbind(Omega,      Ym_PGSp,    thetaNT,   gc,      k+gc,    gt,      gt),
        cbind(Yp_PGSm,    Omega,      thetaT,    gt,      gt,      k+gc,    gc),
        cbind(Yp_PGSm,    Omega,      thetaNT,   gt,      gt,      gc,      k+gc)),
        dimnames=list(colnames(Example_Data),colnames(Example_Data)), name="expCov")

    # Means: fixed at 0 for the binary phenotypes (identification requires fixing either the mean
    # or the threshold; we fix the mean and let the threshold float), free for the continuous PGS.
    Means <- mxMatrix(type = "Full", nrow = 1, ncol = 7, free = c(F, F, F, T, T, T, T), values = 0,
        label = c("meanYp1", "meanYm1", "meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"),
        dimnames = list(NULL, c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1")),
        name = "expMeans")

    # Threshold matrix: Yp1, Ym1 and Yo1 share ONE free threshold (same trait => same prevalence).
    Th <- mxMatrix(type = "Full", nrow = 1, ncol = 3, free = TRUE, values = 0,
        labels = rep("thresh_Y11", 3), name = "Th")
    threshVars <- c("Yp1", "Ym1", "Yo1")

    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans", dimnames=c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1"),
                                            thresholds="Th", threshnames=threshVars)
    Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )
    FitFunctionML   <- mxFitFunctionML()

    Params <- list(
                VY, VE, delta, a, k, j, Omega, Gamma, mu, gt, ht, gc, hc, itlo, itol, ic, f, w, v,
                VY_Algebra, VF_Algebra, Omega_Algebra, Gamma_Algebra, adelta_Constraint_Algebra, j_Algebra, gt_Algebra, ht_Algebra, gc_Algebra, hc_Algebra, gchc_constraint_Algebra, itlo_Algebra, itol_Algebra, ic_Algebra, w_Algebra, v_Algebra, wv_constraint_algebra,
                VY_Constraint, Gamma_Constraint, j_constraint, ht_constraint, hc_constraint, itlo_constraint, itol_constraint, ic_constraint, v_constraint, w_constraint,
                thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Ym_Yp, Yo_Yp, Yo_Ym,
                CovMatrix, Means, Th, ModelExpectations, FitFunctionML)

    Model1 <- mxModel("UniSEM_Binary_SameTrait_ObservedYPYM", Params, Example_Data_Mx)
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1), intervals=T, silent=T, exhaustive = exhaustive, jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
