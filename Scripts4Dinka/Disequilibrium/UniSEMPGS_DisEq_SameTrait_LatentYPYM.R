## This script is a function that fits a version univariate SEM-PGS where parent and offspring have the same trait but the parental phenotypes are latent.
## This script is for trait that is not in equilibrium
## To identify this model, the optimal solution is to use RDR to find a for both parents and offspring.

fitUniSEMPGS_DisEq_SameTrait_LatentParents <- function(data_path, h2_RDR, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
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

    # --- Variances ---
    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=1.5, label="VY11", name="VY", lbound = .001)
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.5, label="VE11", name="VE", lbound = .001)

    # Equation for VY: 2*Omega*delta + 2*Gamma*a + delta*w + a*v + VF + VE
    VY_Algebra <- mxAlgebra(2 * Omega * delta + 2 * Gamma * a + delta * w + a * v + VF_Algebra + VE, name="VY_Algebra")
    # Equation for VF: f*VY*f (assuming father/mother equality)
    VF_Algebra <- mxAlgebra(2 * f^2 * VY, name="VF_Algebra")

    VY_Constraint <- mxConstraint(VY == VY_Algebra, name='VY_Constraint')

    # --- Genetic effects & Path Coefficients ---
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.25, label="delta11", name="delta")
    a     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.7, label="a11", name="a", lbound = .001)
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    j     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")

    # Omega = delta*k + 0.5*w
    # Gamma = a*j + 0.5*v
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="Omega11", name="Omega")
    Gamma <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="Gamma11", name="Gamma")

    Omega_Algebra <- mxAlgebra(delta * k + 0.5 * w , name="Omega_Algebra")
    Gamma_Algebra <- mxAlgebra(a * j + 0.5 * v, name="Gamma_Algebra")

    Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
    Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')

    # RDR Constraint to identify 'a'
    # Current notation: sigma2 is VY, e^2 is VE
    rdr_left  <- mxAlgebra((2*a^2*j + 2*delta^2*k) * (2*a^2*j + 2*delta^2*k + VE), name="rdr_left")
    h2mat     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR, name="h2mat")
    rdr_right <- mxAlgebra(h2mat * VY, name="rdr_right")
    rdrCon    <- mxConstraint(rdr_left == rdr_right, name="rdrCon")

    # --- Assortative Mating (Non-Equilibrium) ---
    mu    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="mu11", name="mu")

    # gc = hc = ic = 0 (No AM in previous generation)
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
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="f11", name="f")

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

    # Expected covariances matrix (5x5: Yo1, Tp1, NTp1, Tm1, NTm1)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY,      thetaT,  thetaNT, thetaT,  thetaNT), # Yo1
        cbind(thetaT,  k+gc,    gc,      gt,      gt),      # Tp1
        cbind(thetaNT, gc,      k+gc,    gt,      gt),      # NTp1
        cbind(thetaT,  gt,      gt,      k+gc,    gc),      # Tm1
        cbind(thetaNT, gt,      gt,      gc,      k+gc)),   # NTm1
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
                VY, VE, delta, a, k, j, Omega, Gamma, mu, gc, hc, f,
                VY_Algebra, VF_Algebra, Omega_Algebra, Gamma_Algebra,
                gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
                gt, ht, ic, w, v,
                VY_Constraint, Omega_Constraint, Gamma_Constraint,
                h2mat, rdr_left, rdr_right, rdrCon, # Identification for 'a'
                thetaNT, thetaT,
                CovMatrix, Means, ModelExpectations, FitFunctionML)

    Model1 <- mxModel("UniSEM_NonEquilibrium_LatentParents", Params, Example_Data_Mx)
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
