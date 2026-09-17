## Binary/liability-threshold version of UniSEMPGS_DisEq_SameTrait_ObservedYPYM.R: parent and
## offspring express the SAME binary trait (liability-threshold model), the parental phenotypes are
## observed, and the trait is NOT in equilibrium (gc = hc = ic = 0, one generation of assortative
## mating). Yp1/Ym1/Yo1 are dichotomous (0/1) and share a single free threshold (same trait => same
## prevalence); Tp1/NTp1/Tm1/NTm1 stay continuous. Liability variance VY is fixed to 1.

fitUniSEMPGS_DisEq_Binary_SameTrait_ObservedYPYM <- function(data_path, feaTol = 1e-6, optTol = 1e-8, jitterMean = .2, jitterVar = .05, extraTries = 30, exhaustive = F){
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
    for (v in c("Yp1","Ym1","Yo1")) Example_Data[[v]] <- mxFactor(Example_Data[[v]], levels = c(0, 1))

    # --- Variances ---
    # Liability variance is FIXED at 1 (not free): identifies the scale of the binary/threshold model.
    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1, label="VY11", name="VY")
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="VE11", name="VE", lbound = .001)

    # Equation for VY from screenshot: 2*Omega*delta + 2*Gamma*a + delta*w + a*v + VF + VE
    VY_Algebra <- mxAlgebra(2 * Omega * delta + 2 * Gamma * a + delta * w + a * v + VF_Algebra + VE, name="VY_Algebra")
    # Equation for VF: f*VY*f (assuming father/mother equality)
    VF_Algebra <- mxAlgebra(2 * f^2 * VY, name="VF_Algebra")

    VY_Constraint <- mxConstraint(VY == VY_Algebra, name='VY_Constraint')

    # --- Genetic effects & Path Coefficients ---
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="delta11", name="delta")
    a     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="a11", name="a", lbound = .001)
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    j     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")

    # Omega and Gamma Equations from screenshot
    # Omega = delta*k + 0.5*w
    # Gamma = a*j + 0.5*v
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="Omega11", name="Omega")
    Gamma <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="Gamma11", name="Gamma")

    Omega_Algebra <- mxAlgebra(delta * k + 0.5 * w , name="Omega_Algebra")
    Gamma_Algebra <- mxAlgebra(a * j + 0.5 * v, name="Gamma_Algebra")

    Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
    Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')

    # --- Assortative Mating (Non-Equilibrium) ---
    mu    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="mu11", name="mu")

    # From screenshot: gc = hc = ic = 0 (No AM in previous generation)
    gc    <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=0, label="gc11", name="gc")
    hc    <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=0, label="hc11", name="hc")

    # gt, ht, ic are functions of mu and current Omega/Gamma. (itlo and itol -- the
    # transmitted-vs-non-transmitted haplotype covariances -- collapse to a single quantity 'ic'
    # here because Gamma*mu*Omega == Omega*mu*Gamma for scalars; the itlo/itol split only matters
    # once Omega and Gamma become asymmetric matrices in the bivariate model.)
    gt_Algebra <- mxAlgebra(Omega * mu * Omega, name="gt_Algebra")
    ht_Algebra <- mxAlgebra(Gamma * mu * Gamma, name="ht_Algebra")
    ic_Algebra <- mxAlgebra(Omega * mu * Gamma, name="ic_Algebra")

    gt <- mxAlgebra(gt_Algebra, name="gt")
    ht <- mxAlgebra(ht_Algebra, name="ht")
    ic <- mxAlgebra(ic_Algebra, name="ic")

    # --- Vertical Transmission ---
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="f11", name="f")

    # Equation: w = f*Omega + f*Omega (since parent traits/effects are equal)
    w_Algebra     <- mxAlgebra(2 * f * Omega, name="w_Algebra")
    v_Algebra     <- mxAlgebra(2 * f * Gamma, name="v_Algebra")

    w <- mxAlgebra(w_Algebra, name="w")
    v <- mxAlgebra(v_Algebra, name="v")

    # --- Offspring Covariances ---
    # thetaNT = 2*a*ic + 2*delta*gt + w   (ic == itlo == itol in the univariate model)
    thetaNT <- mxAlgebra(2*a*ic + 2*delta*gt + w, name="thetaNT")
    # thetaT = 2*delta*k + thetaNT
    thetaT  <- mxAlgebra(2*delta*k + thetaNT, name="thetaT")

    # Cross-parental expectations
    Yp_PGSm <- mxAlgebra(VY * mu * Omega, name="Yp_PGSm")
    Ym_PGSp <- mxAlgebra(VY * mu * Omega, name="Ym_PGSp")
    Yp_Ym   <- mxAlgebra(VY * mu * VY,    name="Yp_Ym")

    # Offspring-Parent: delta*Omega + a*Gamma + f*VY (Simplified for one generation)
    Yo_Yp   <- mxAlgebra(delta*Omega + a*Gamma + f*VY + (delta*Omega + a*Gamma)*mu*VY, name = "Yo_Yp")

    # --- Expected Covariance Matrix ---
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY,      Yp_Ym,   Yo_Yp,   Omega,   Omega,   Yp_PGSm, Yp_PGSm),
        cbind(Yp_Ym,   VY,      Yo_Yp,   Ym_PGSp, Ym_PGSp, Omega,   Omega),
        cbind(Yo_Yp,   Yo_Yp,   VY,      thetaT,  thetaNT, thetaT,  thetaNT),
        cbind(Omega,   Ym_PGSp, thetaT,  k+gc,    gc,      gt,      gt),
        cbind(Omega,   Ym_PGSp, thetaNT, gc,      k+gc,    gt,      gt),
        cbind(Yp_PGSm, Omega,   thetaT,  gt,      gt,      k+gc,    gc),
        cbind(Yp_PGSm, Omega,   thetaNT, gt,      gt,      gc,      k+gc)),
        dimnames=list(colnames(Example_Data),colnames(Example_Data)), name="expCov")

    # Yp1/Ym1/Yo1 means fixed at 0 (identification requires fixing either the mean or the
    # threshold); the PGS means stay free.
    Means <- mxMatrix(type = "Full", nrow = 1, ncol = 7, free = c(F, F, F, T, T, T, T), values = 0,
        label = c("meanYp1", "meanYm1", "meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"),
        dimnames = list(NULL, c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1")),
        name = "expMeans")

    # Yp1, Ym1 and Yo1 share ONE free threshold (same trait => same prevalence).
    Th <- mxMatrix(type = "Full", nrow = 1, ncol = 3, free = TRUE, values = 0,
        labels = rep("thresh_Y11", 3), name = "Th")

    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans", dimnames=colnames(Example_Data),
                                            thresholds="Th", threshnames=c("Yp1","Ym1","Yo1"))
    Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )
    FitFunctionML   <- mxFitFunctionML()

    Params <- list(
                VY, VE, delta, a, k, j, Omega, Gamma, mu, gc, hc, f,
                VY_Algebra, VF_Algebra, Omega_Algebra, Gamma_Algebra,
                gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
                gt, ht, ic, w, v,
                VY_Constraint, Omega_Constraint, Gamma_Constraint,
                thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Yo_Yp,
                CovMatrix, Means, Th, ModelExpectations, FitFunctionML)

    Model1 <- mxModel("UniSEM_Binary_NonEquilibrium", Params, Example_Data_Mx)
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1), intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
