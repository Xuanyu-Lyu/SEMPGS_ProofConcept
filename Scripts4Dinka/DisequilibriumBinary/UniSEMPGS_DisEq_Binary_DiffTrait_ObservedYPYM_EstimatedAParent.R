## Binary/liability-threshold version of UniSEMPGS_DisEq_DiffTrait_ObservedYPYM_EstimatedAParent.R:
## parent and offspring express DIFFERENT binary traits, the parental phenotypes ARE observed, and
## the trait is NOT in equilibrium (gc = hc = ic = 0, one generation of assortative mating). Yp1/Ym1
## (parent trait) share one free threshold; Yo1 (offspring trait) gets its own. RDR identifies the
## offspring's 'a'; the parental 'a' is estimated from the observed parental data. Each trait's
## liability variance (VY_p, VY_o) is independently fixed to 1.

fitUniSEMPGS_DisEq_Binary_DiffTrait_ObservedYPYM_EstimatedAParent <- function(data_path, h2_RDR_offspring, feaTol = 1e-6, optTol = 1e-8, jitterMean = .2, jitterVar = .05, extraTries = 30, exhaustive = F){
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

    # 1. Phenotypic and Residual Variances (Separated)
    # Each trait's liability variance is FIXED at 1 (not free), independently for parent/offspring.
    VY_p  <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1,   label="VY_parent", name="VY_p")
    VY_o  <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1,   label="VY_offspring", name="VY_o")
    VE_p  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4,  label="VE_parent", name="VE_p", lbound = .001)
    VE_o  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4,  label="VE_offspring", name="VE_o", lbound = .001)

    # 2. Genetic effects (Generation Specific)
    delta_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="delta_p11", name="delta_p")
    a_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="a_p11", name="a_p", lbound = .001)

    delta_o <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="delta_o11", name="delta_o")
    a_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="a_o11", name="a_o", lbound = .001)

    k <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    j <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")

    # 3. Covariances and Assortment
    Omega_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="Omega_p11", name="Omega_p")
    Gamma_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="Gamma_p11", name="Gamma_p")

    mu <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="mu11", name="mu")

    # gc = hc = 0 (No AM in previous generation)
    gc <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=0, label="gc11", name="gc")
    hc <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=0, label="hc11", name="hc")

    # 4. Vertical Transmission
    f <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="f11", name="f")

    # w = f*Omega_p + f*Omega_p (since parent traits/effects are equal); v = f*Gamma_p + f*Gamma_p
    w_Algebra <- mxAlgebra(2 * f * Omega_p, name="w_Algebra")
    v_Algebra <- mxAlgebra(2 * f * Gamma_p, name="v_Algebra")

    w <- mxAlgebra(w_Algebra, name="w")
    v <- mxAlgebra(v_Algebra, name="v")

    # 5. Scalar Algebra for Parent Trait (Non-Equilibrium)
    VF_p_Algebra <- mxAlgebra(2 * f^2 * VY_p, name="VF_p_Algebra")
    VY_p_Algebra <- mxAlgebra(2 * delta_p * Omega_p + 2 * a_p * Gamma_p + w * delta_p + v * a_p + VF_p_Algebra + VE_p, name="VY_p_Algebra")

    Omega_p_Algebra <- mxAlgebra(delta_p * k + 0.5 * w, name="Omega_p_Algebra")
    Gamma_p_Algebra <- mxAlgebra(a_p * j + 0.5 * v, name="Gamma_p_Algebra")

    # 6. Scalar Algebra for Offspring Trait (gc = hc = ic = 0)
    VY_o_Algebra <- mxAlgebra(2 * delta_o^2 * k + 2 * a_o^2 * j + 2 * delta_o * w + 2 * a_o * v + 2 * f^2 * VY_p + VE_o, name="VY_o_Algebra")

    # Identification of a_o via RDR
    rdr_left_o  <- mxAlgebra((2*a_o^2*j + 2*delta_o^2*k) * (2*a_o^2*j + 2*delta_o^2*k + VE_o), name="rdr_left_o")
    h2mat_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_offspring, name="h2mat_o")
    rdr_right_o <- mxAlgebra(h2mat_o * VY_o, name="rdr_right_o")

    # 7. Assortment and Transmission Algebra (functions of mu and current Omega_p/Gamma_p). (itlo
    # and itol -- the transmitted-vs-non-transmitted haplotype covariances -- collapse to a single
    # quantity 'ic' here because Gamma_p*mu*Omega_p == Omega_p*mu*Gamma_p for scalars; the
    # itlo/itol split only matters once Omega_p and Gamma_p become asymmetric bivariate matrices.)
    gt_Algebra <- mxAlgebra(Omega_p * mu * Omega_p, name="gt_Algebra")
    ht_Algebra <- mxAlgebra(Gamma_p * mu * Gamma_p, name="ht_Algebra")
    ic_Algebra <- mxAlgebra(Omega_p * mu * Gamma_p, name="ic_Algebra")

    gt <- mxAlgebra(gt_Algebra, name="gt")
    ht <- mxAlgebra(ht_Algebra, name="ht")
    ic <- mxAlgebra(ic_Algebra, name="ic")

    # thetaNT = 2*a_o*ic + 2*delta_o*gt + w   (ic == itlo == itol in the univariate model)
    thetaNT <- mxAlgebra(2*a_o*ic + 2*delta_o*gt + w, name="thetaNT")
    # thetaT = 2*delta_o*k + thetaNT
    thetaT  <- mxAlgebra(2*delta_o*k + thetaNT, name="thetaT")

    # Intergenerational Phenotypic Covariances (Between Parent Trait and Child Trait)
    Yp_Ym   <- mxAlgebra(VY_p * mu * VY_p, name="Yp_Ym")
    Yo_Yp   <- mxAlgebra(delta_o * Omega_p + a_o * Gamma_p + f * VY_p + (delta_o * Omega_p + a_o * Gamma_p) * mu * VY_p, name = "Yo_Yp")

    # Cross-person PGS-Phenotype Covariances
    Yp_PGSm <- mxAlgebra(VY_p * mu * Omega_p, name="Yp_PGSm")

    # 8. Expected Covariance Matrix (7x7)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_p_Algebra, Yp_Ym,      Yo_Yp,     Omega_p, Omega_p, Yp_PGSm, Yp_PGSm), # Yp
        cbind(Yp_Ym,        VY_p_Algebra, Yo_Yp,   Yp_PGSm, Yp_PGSm, Omega_p, Omega_p), # Ym
        cbind(Yo_Yp,        Yo_Yp,      VY_o_Algebra, thetaT, thetaNT, thetaT,  thetaNT), # Yo
        cbind(Omega_p,      Yp_PGSm,    thetaT,    k+gc,    gc,      gt,      gt), # Tp
        cbind(Omega_p,      Yp_PGSm,    thetaNT,   gc,      k+gc,    gt,      gt), # NTp
        cbind(Yp_PGSm,      Omega_p,    thetaT,    gt,      gt,      k+gc,    gc), # Tm
        cbind(Yp_PGSm,      Omega_p,    thetaNT,   gt,      gt,      gc,      k+gc)), # NTm
        dimnames=list(colnames(Example_Data), colnames(Example_Data)), name="expCov")

    # 9. Means and Expectations
    # Yp1/Ym1/Yo1 means fixed at 0 (identification requires fixing either the mean or the
    # threshold); the PGS means stay free.
    Means <- mxMatrix(type = "Full", nrow = 1, ncol = 7, free = c(F, F, F, T, T, T, T), values = 0,
        label = c("meanYp1", "meanYm1", "meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"),
        dimnames = list(NULL, colnames(Example_Data)), name = "expMeans")

    # Yp1/Ym1 share one free threshold (same parent trait); Yo1 (a different trait) gets its own.
    Th <- mxMatrix(type = "Full", nrow = 1, ncol = 3, free = TRUE, values = 0,
        labels = c("thresh_Yp11", "thresh_Yp11", "thresh_Yo11"), name = "Th")

    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans", dimnames=colnames(Example_Data),
                                            thresholds="Th", threshnames=c("Yp1","Ym1","Yo1"))

    # 10. Constraints
    Constraints <- list(
        mxConstraint(VY_p == VY_p_Algebra, name="VY_p_eq"),
        mxConstraint(VY_o == VY_o_Algebra, name="VY_o_eq"),
        mxConstraint(Omega_p == Omega_p_Algebra, name="Om_p_eq"),
        mxConstraint(Gamma_p == Gamma_p_Algebra, name="Ga_p_eq"),
        mxConstraint(rdr_left_o == rdr_right_o, name="rdr_o_con")
    )

    Params <- list(
        VY_p, VY_o, VE_p, VE_o, delta_p, a_p, delta_o, a_o, k, j, Omega_p, Gamma_p, mu, gc, hc, f,
        VY_p_Algebra, VF_p_Algebra, VY_o_Algebra, Omega_p_Algebra, Gamma_p_Algebra,
        gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
        gt, ht, ic, w, v,
        h2mat_o, rdr_left_o, rdr_right_o, thetaNT, thetaT, Yp_Ym, Yo_Yp, Yp_PGSm,
        CovMatrix, Means, Th, ModelExpectations, mxFitFunctionML(), Constraints
    )

    Model1 <- mxModel("UniSEM_Binary_NonEquilibrium_DiffTrait_Observed_RDR", Params, mxData(observed=Example_Data, type="raw"))
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
