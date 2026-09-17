## This script is a function that fits a version univariate SEM-PGS where parent and offspring have different traits and the parental phenotypes are latent.
## To identify this model, the optimal solution is to use RDR to find two different 'a' parameters for both parents and offspring and to equate the phenotypic variance of the latent parental traits to the phenotypic variance of the observed offspring trait.
##
## IDENTIFICATION NOTE: with both traits fully latent on the parent side, the 5-variable data
## (Yo1, Tp1, NTp1, Tm1, NTm1) supply only 4 independent moments (VY, thetaT, thetaNT, gt), one
## fewer than this model's free dimensions even after both RDR constraints. Concretely, gc_eq only
## pins down the PRODUCT Omega_p^2*mu (via gc == Omega_p*mu*Omega_p), not Omega_p and mu
## individually -- fixing an unrelated parameter (e.g. 'f') does NOT resolve this, since the
## Omega_p/mu split stays free either way (confirmed empirically: fixing f converges cleanly to a
## confidently WRONG point, not merely a noisy one). Passing an externally-known `mu_fixed` (e.g.
## the known/assumed assortative-mating co-path for this trait) breaks exactly that degeneracy --
## Omega_p becomes uniquely determined by gc == Omega_p^2*mu_fixed -- and makes the model
## point-identified (verified to recover every parameter to ~1e-4); leaving it NULL reproduces the
## original (under-identified) behavior.

fitUniSEMPGS_DiffTrait_LatentYPYM <- function(data_path, h2_RDR_parent, h2_RDR_offspring, mu_fixed = NULL, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
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
    delta_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="delta_p11", name="delta_p", lbound = .001) 
    a_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.6, label="a_p11", name="a_p", lbound = .001)     
    
    delta_o <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="delta_o11", name="delta_o") 
    a_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.6, label="a_o11", name="a_o", lbound = .001)     
    
    k <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")     
    j <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")     

    # 3. Covariances and Assortment
    Omega_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="Omega_p11", name="Omega_p", lbound = .001) 
    Gamma_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="Gamma_p11", name="Gamma_p") 
    
    # 'mu' is free unless an external value is supplied via mu_fixed (see IDENTIFICATION NOTE above).
    mu <- mxMatrix(type="Full", nrow=1, ncol=1, free=is.null(mu_fixed),
                   values=if (is.null(mu_fixed)) .1 else mu_fixed, label="mu11", name="mu")
    ic <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.02, label="ic11", name="ic")
    gc <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.05, label="gc11", name="gc")
    hc <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.01, label="hc11", name="hc")  

    # 4. Vertical Transmission
    f <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="f11", name="f", lbound = .001)
    w <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="w11", name="w")
    v <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="v11", name="v") 

    # 5. Scalar Algebra
    # Parental Equilibrium logic
    VF_p_Algebra <- mxAlgebra(2 * f^2 * VY + 2 * f^2 * VY^2 * mu, name="VF_p_Algebra")
    VY_p_Algebra <- mxAlgebra(2 * delta_p * Omega_p + 2 * a_p * Gamma_p + w * delta_p + v * a_p + VF_p_Algebra + VE_p, name="VY_p_Algebra")
    
    # Offspring Variance logic
    VY_o_Algebra <- mxAlgebra(2 * delta_o^2 * k + 4 * delta_o^2 * gc + 2 * a_o^2 * j + 4 * a_o^2 * hc + 8 * delta_o * ic * a_o + 2 * delta_o * w + 2 * a_o * v + 2 * f^2 * VY * (1 + VY * mu) + VE_o, name="VY_o_Algebra")

    # Within-person PGS-Phenotype Covariance (Parent Trait)
    Omega_p_Algebra <- mxAlgebra(2 * delta_p * gc + 2 * a_p * ic + delta_p * k + 0.5 * w, name="Omega_p_Algebra") 
    Gamma_p_Algebra <- mxAlgebra(2 * a_p * hc + 2 * delta_p * ic + a_p * j + 0.5 * v, name="Gamma_p_Algebra") 

    # Identification via RDR heritability
    rdr_left_p  <- mxAlgebra((2*a_p^2*j + 2*delta_p^2*k) * (2*a_p^2*j + 2*delta_p^2*k + VE_p), name="rdr_left_p")
    h2mat_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_parent, name="h2mat_p")
    rdr_right_p <- mxAlgebra(h2mat_p * VY, name="rdr_right_p")

    rdr_left_o  <- mxAlgebra((2*a_o^2*j + 2*delta_o^2*k) * (2*a_o^2*j + 2*delta_o^2*k + VE_o), name="rdr_left_o")
    h2mat_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_offspring, name="h2mat_o")
    rdr_right_o <- mxAlgebra(h2mat_o * VY, name="rdr_right_o")

    # 6. Assortative Mating PGS Covariances (Equilibrium)
    gt_Algebra <- mxAlgebra(Omega_p * mu * Omega_p, name="gt_Algebra")
    ht_Algebra <- mxAlgebra(Gamma_p * mu * Gamma_p, name="ht_Algebra")

    # Constraints tying ic, w and v to the equilibrium recursion (as in the same-trait script)
    ic_Algebra <- mxAlgebra(.5 * (Gamma_p * mu * Omega_p + Omega_p * mu * Gamma_p), name="ic_Algebra")
    w_Algebra  <- mxAlgebra(2 * f * Omega_p + 2 * f * VY * mu * Omega_p, name="w_Algebra")
    v_Algebra  <- mxAlgebra(2 * f * Gamma_p + 2 * f * VY * mu * Gamma_p, name="v_Algebra")

    # 7. Offspring-PGS Covariances
    thetaNT <- mxAlgebra(2 * delta_o * gc + 2 * a_o * ic + .5 * w, name="thetaNT")
    thetaT  <- mxAlgebra(delta_o * k + thetaNT, name="thetaT")

    # 8. Expected Covariance Matrix (5x5: Yo1, Tp1, NTp1, Tm1, NTm1)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_o_Algebra, thetaT,  thetaNT, thetaT,  thetaNT), 
        cbind(thetaT,       k+gc,    gc,      gt_Algebra, gt_Algebra), 
        cbind(thetaNT,      gc,      k+gc,    gt_Algebra, gt_Algebra), 
        cbind(thetaT,       gt_Algebra, gt_Algebra, k+gc,    gc),      
        cbind(thetaNT,      gt_Algebra, gt_Algebra, gc,      k+gc)),   
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
        mxConstraint(gc == gt_Algebra, name="gc_eq"),
        mxConstraint(hc == ht_Algebra, name="hc_eq"),
        mxConstraint(ic == ic_Algebra, name="ic_eq"),
        mxConstraint(w == w_Algebra, name="w_eq"),
        mxConstraint(v == v_Algebra, name="v_eq"),
        mxConstraint(rdr_left_p == rdr_right_p, name="rdr_p_con"),
        mxConstraint(rdr_left_o == rdr_right_o, name="rdr_o_con")
    )

    Params <- list(
        VY, VE_p, VE_o, delta_p, a_p, delta_o, a_o, k, j, Omega_p, Gamma_p, mu, ic, gc, hc, f, w, v,
        VY_p_Algebra, VF_p_Algebra, VY_o_Algebra, Omega_p_Algebra, Gamma_p_Algebra,
        gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
        h2mat_p, h2mat_o, rdr_left_p, rdr_right_p, rdr_left_o, rdr_right_o,
        thetaNT, thetaT, CovMatrix, Means, ModelExpectations, mxFitFunctionML(), Constraints
    )

    Model1 <- mxModel("UniSEM_DiffTrait_Latent_Final", Params, mxData(observed=Example_Data, type="raw"))
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1), 
                          intervals=T, silent=T, exhaustive = exhaustive, 
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)
    
    return(summary(fitModel1, verbose = TRUE))
}