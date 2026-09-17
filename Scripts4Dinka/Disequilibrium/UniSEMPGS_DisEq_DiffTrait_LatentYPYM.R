## This script is a function that fits a version univariate SEM-PGS where parent and offspring have different traits and the parental phenotypes are latent.
## This script is for trait that is not in equilibrium
## To identify this model, the optimal solution is to use RDR to find two different 'a' parameters for both parents and offspring and to equate the phenotypic variance of the latent parental traits to the phenotypic variance of the observed offspring trait.
##
## IDENTIFICATION NOTE: with both traits fully latent on the parent side, the 5-variable data
## (Yo1, Tp1, NTp1, Tm1, NTm1) supply only 4 independent moments (VY, thetaT, thetaNT, gt), one
## fewer than this model's free dimensions even after both RDR constraints. Concretely, the
## gt(=gc) moment only pins down the PRODUCT Omega_p^2*mu (via gt_Algebra = Omega_p*mu*Omega_p),
## not Omega_p and mu individually -- fixing an unrelated parameter (e.g. 'f') does NOT resolve
## this, since the Omega_p/mu split stays free either way (confirmed empirically in the
## equilibrium version: fixing f converges cleanly to a confidently WRONG point, not merely a
## noisy one). Passing an externally-known `mu_fixed` (e.g. the known/assumed assortative-mating
## co-path for this trait) breaks exactly that degeneracy -- Omega_p becomes uniquely determined by
## gt == Omega_p^2*mu_fixed -- and makes the model point-identified; leaving it NULL reproduces the
## original (under-identified) behavior.

fitUniSEMPGS_DisEq_DiffTrait_LatentYPYM <- function(data_path, h2_RDR_parent, h2_RDR_offspring, mu_fixed = NULL, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
    library(OpenMx)
    library(data.table)
    library(stringr)

    if (is.null(mu_fixed)) warning("DiffTrait_LatentYPYM is under-identified (see ConceptProof.md) unless mu_fixed is supplied; 'mu' will be freely estimated and the fit may converge to a confidently wrong point (not just large SEs).")

    mxOption(NULL,"Calculate Hessian","Yes")
    mxOption(NULL,"Standard Errors","Yes")
    mxOption(NULL,"Default optimizer","NPSOL")
    mxOption(NULL,"Feasibility tolerance",as.character(feaTol))
    mxOption(NULL,"Optimality tolerance",as.character(optTol))
    mxOption(NULL,"Number of Threads", value = parallel::detectCores())

    Example_Data  <- fread(data_path, header = T)

    # 1. Phenotypic and Residual Variances
    # We define one VY matrix and equate both generations to it for identification
    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=1.5, label="VY11", name="VY", lbound = .001)
    VE_p  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.5,  label="VE_parent", name="VE_p", lbound = .001)
    VE_o  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.5,  label="VE_offspring", name="VE_o", lbound = .001)

    # 2. Genetic effects (Generation Specific)
    delta_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.25, label="delta_p11", name="delta_p", lbound = .001)
    a_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.7, label="a_p11", name="a_p", lbound = .001)

    delta_o <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="delta_o11", name="delta_o")
    a_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.65, label="a_o11", name="a_o", lbound = .001)

    k <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    j <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")

    # 3. Covariances and Assortment
    Omega_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="Omega_p11", name="Omega_p", lbound = .001)
    Gamma_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="Gamma_p11", name="Gamma_p")

    # 'mu' is free unless an external value is supplied via mu_fixed (see IDENTIFICATION NOTE above).
    mu <- mxMatrix(type="Full", nrow=1, ncol=1, free=is.null(mu_fixed),
                   values=if (is.null(mu_fixed)) .15 else mu_fixed, label="mu11", name="mu")

    # gc = hc = 0 (No AM in previous generation)
    gc <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=0, label="gc11", name="gc")
    hc <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=0, label="hc11", name="hc")

    # 4. Vertical Transmission
    f <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="f11", name="f", lbound = .001)

    # w = f*Omega_p + f*Omega_p (since parent traits/effects are equal); v = f*Gamma_p + f*Gamma_p
    w_Algebra <- mxAlgebra(2 * f * Omega_p, name="w_Algebra")
    v_Algebra <- mxAlgebra(2 * f * Gamma_p, name="v_Algebra")

    w <- mxAlgebra(w_Algebra, name="w")
    v <- mxAlgebra(v_Algebra, name="v")

    # 5. Scalar Algebra
    # Parental Non-Equilibrium logic
    VF_p_Algebra <- mxAlgebra(2 * f^2 * VY, name="VF_p_Algebra")
    VY_p_Algebra <- mxAlgebra(2 * delta_p * Omega_p + 2 * a_p * Gamma_p + w * delta_p + v * a_p + VF_p_Algebra + VE_p, name="VY_p_Algebra")

    # Offspring Variance logic (gc = hc = ic = 0)
    VY_o_Algebra <- mxAlgebra(2 * delta_o^2 * k + 2 * a_o^2 * j + 2 * delta_o * w + 2 * a_o * v + 2 * f^2 * VY + VE_o, name="VY_o_Algebra")

    # Within-person PGS-Phenotype Covariance (Parent Trait)
    Omega_p_Algebra <- mxAlgebra(delta_p * k + 0.5 * w, name="Omega_p_Algebra")
    Gamma_p_Algebra <- mxAlgebra(a_p * j + 0.5 * v, name="Gamma_p_Algebra")

    # Identification via RDR heritability
    rdr_left_p  <- mxAlgebra((2*a_p^2*j + 2*delta_p^2*k) * (2*a_p^2*j + 2*delta_p^2*k + VE_p), name="rdr_left_p")
    h2mat_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_parent, name="h2mat_p")
    rdr_right_p <- mxAlgebra(h2mat_p * VY, name="rdr_right_p")

    rdr_left_o  <- mxAlgebra((2*a_o^2*j + 2*delta_o^2*k) * (2*a_o^2*j + 2*delta_o^2*k + VE_o), name="rdr_left_o")
    h2mat_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_offspring, name="h2mat_o")
    rdr_right_o <- mxAlgebra(h2mat_o * VY, name="rdr_right_o")

    # 6. Assortative Mating PGS Covariances (functions of mu and current Omega_p/Gamma_p). (itlo
    # and itol -- the transmitted-vs-non-transmitted haplotype covariances -- collapse to a single
    # quantity 'ic' here because Gamma_p*mu*Omega_p == Omega_p*mu*Gamma_p for scalars; the
    # itlo/itol split only matters once Omega_p and Gamma_p become asymmetric bivariate matrices.)
    gt_Algebra <- mxAlgebra(Omega_p * mu * Omega_p, name="gt_Algebra")
    ht_Algebra <- mxAlgebra(Gamma_p * mu * Gamma_p, name="ht_Algebra")
    ic_Algebra <- mxAlgebra(Omega_p * mu * Gamma_p, name="ic_Algebra")

    gt <- mxAlgebra(gt_Algebra, name="gt")
    ht <- mxAlgebra(ht_Algebra, name="ht")
    ic <- mxAlgebra(ic_Algebra, name="ic")

    # 7. Offspring-PGS Covariances
    # thetaNT = 2*a_o*ic + 2*delta_o*gt + w   (ic == itlo == itol in the univariate model)
    thetaNT <- mxAlgebra(2*a_o*ic + 2*delta_o*gt + w, name="thetaNT")
    # thetaT = 2*delta_o*k + thetaNT
    thetaT  <- mxAlgebra(2*delta_o*k + thetaNT, name="thetaT")

    # 8. Expected Covariance Matrix (5x5: Yo1, Tp1, NTp1, Tm1, NTm1)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_o_Algebra, thetaT,  thetaNT, thetaT,  thetaNT),
        cbind(thetaT,       k+gc,    gc,      gt,      gt),
        cbind(thetaNT,      gc,      k+gc,    gt,      gt),
        cbind(thetaT,       gt,      gt,      k+gc,    gc),
        cbind(thetaNT,      gt,      gt,      gc,      k+gc)),
        dimnames = list(c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1"),
                        c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")), name="expCov")

    Means <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = TRUE, values = 0,
        label = c("meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"),
        dimnames = list(NULL, c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")), name = "expMeans")

    # 9. Final Model and Constraints
    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans",
                                            dimnames=c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1"))

    Constraints <- list(
        mxConstraint(VY == VY_p_Algebra, name="VY_p_eq"),
        mxConstraint(VY == VY_o_Algebra, name="VY_o_eq"),
        mxConstraint(Omega_p == Omega_p_Algebra, name="Om_p_eq"),
        mxConstraint(Gamma_p == Gamma_p_Algebra, name="Ga_p_eq"),
        mxConstraint(rdr_left_p == rdr_right_p, name="rdr_p_con"),
        mxConstraint(rdr_left_o == rdr_right_o, name="rdr_o_con")
    )

    Params <- list(
        VY, VE_p, VE_o, delta_p, a_p, delta_o, a_o, k, j, Omega_p, Gamma_p, mu, gc, hc, f,
        VY_p_Algebra, VF_p_Algebra, VY_o_Algebra, Omega_p_Algebra, Gamma_p_Algebra,
        gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
        gt, ht, ic, w, v,
        h2mat_p, h2mat_o, rdr_left_p, rdr_right_p, rdr_left_o, rdr_right_o,
        thetaNT, thetaT, CovMatrix, Means, ModelExpectations, mxFitFunctionML(), Constraints
    )

    Model1 <- mxModel("UniSEM_NonEquilibrium_DiffTrait_Latent", Params, mxData(observed=Example_Data, type="raw"))
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
