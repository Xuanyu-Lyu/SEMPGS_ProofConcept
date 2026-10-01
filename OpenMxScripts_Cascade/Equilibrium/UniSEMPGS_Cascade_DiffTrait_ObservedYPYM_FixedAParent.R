## SEM-PGS CASCADE version of OpenMxScripts/Equilibrium/UniSEMPGS_DiffTrait_ ObservedYPYM_FixedAParent.R: the same model, except that spouses
## assort on the latent mating phenotype gamma~ = delta~*(T+NT) + a~*(LT+LNT) + 1~*F + 1~*E instead of on Y
## (SEM_PGS_Cascade_model.pdf, Part III). Following the ETFD Cascade MVN scripts, each path into
## gamma~ is the path into Y times a multiplier: AM_G on the genetic paths (delta~ = AM_G*delta, a~ = AM_G*a)
## and AM_E on the non-genetic paths (the PDF's 1~ on F and E). AM_G = AM_E = 1: primary phenotypic AM;
## AM_E = 0: genetic homogamy; AM_G = 0: social homogamy. gamma~ has no scale of its own, so at least one
## multiplier must be fixed (default AM_E = 1, the analogue of the MVN script's fixed AM_U); mu is the
## copath between spouses' gamma~. By default AM_G is estimated; for genetic homogamy fix AM_G = 1 and free AM_E instead.
##
## This script is a function that fits a version univariate SEM-PGS where parent and offspring have different traits and the parental phenotypes are observed 
## To identify this model, the optimal solution is to use RDR to find  'a' for offspring and use both RDR for parental 'a' and the rest of info to estimate 'a' for parents.

fitUniSEMPGS_Cascade_DiffTrait_ObservedYPYM_FixedAParent <- function(data_path, h2_RDR_parent, h2_RDR_offspring, AM_G_value = .5, AM_G_free = TRUE, AM_E_value = 1, AM_E_free = FALSE, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
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

    # 1. Phenotypic and Residual Variances (Independently estimated)
    VY_p  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=1.5, label="VY_p11", name="VY_p", lbound = .001) 
    VY_o  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=1.5, label="VY_o11", name="VY_o", lbound = .001) 
    VE_p  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.5,  label="VE_p11", name="VE_p", lbound = .001) 
    VE_o  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.5,  label="VE_o11", name="VE_o", lbound = .001) 

    # 2. Genetic effects (Generation Specific)
    delta_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="delta_p11", name="delta_p") 
    a_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.6, label="a_p11", name="a_p", lbound = .001)     
    
    delta_o <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="delta_o11", name="delta_o") 
    a_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.6, label="a_o11", name="a_o", lbound = .001)     
    
    k <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")     
    j <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")     
    # k and j are fixed at .5, which assumes the PGS (and, by convention, the latent LGS) is standardised
    # in the base population. Empirical PGS are usually scaled in the analysed sample, i.e. at equilibrium,
    # where the haplotypic PGS variance is k + g and the full PGS variance is 2k + 4g (g = gc = gt, h = hc = ht).
    # To use one of those scalings, replace the two fixed matrices above with the matching algebra pair:
    #   haplotypic PGS scaled to variance 1/2 at equilibrium:  k = 1/2 - g,   j = 1/2 - h
    # k <- mxAlgebra(.5 - gc, name = "k")
    # j <- mxAlgebra(.5 - hc, name = "j")
    #   full PGS standardised to variance 1 at equilibrium:    k = 1/2 - 2g,  j = 1/2 - 2h
    # k <- mxAlgebra(.5 - 2 * gc, name = "k")
    # j <- mxAlgebra(.5 - 2 * hc, name = "j")

    # 3. Covariances and Assortment
    Omega_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="Omega_p11", name="Omega_p") 
    Gamma_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="Gamma_p11", name="Gamma_p") 
    
    mu <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="mu11", name="mu") 
    ic <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.02, label="ic11", name="ic") 
    gc <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.05, label="gc11", name="gc")  
    hc <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.01, label="hc11", name="hc")  

    # 4. Vertical Transmission
    f <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="f11", name="f") 
    w <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="w11", name="w") 
    v <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="v11", name="v") 

    # 5. Scalar Algebra for Parent Trait (Equilibrium context for the trait itself)
    VF_p_Algebra <- mxAlgebra(2 * f^2 * VY_p + 2 * f^2 * tau^2 * mu, name="VF_p_Algebra")
    # Cascade: VF = 2f^2(VY + tau^2*mu) depends on tau = cov(Y, gamma~), which itself depends on VF,
    # so VF is a free parameter held to its algebra by a constraint (as VY is).
    VF    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="VF11", name="VF")
    VF_Constraint <- mxConstraint(VF == VF_p_Algebra, name='VF_Constraint')
    VY_p_Algebra <- mxAlgebra(2 * delta_p * Omega_p + 2 * a_p * Gamma_p + w * delta_p + v * a_p + VF + VE_p, name="VY_p_Algebra")
    
    Omega_p_Algebra <- mxAlgebra(2 * delta_p * gc + 2 * a_p * ic + delta_p * k + 0.5 * w, name="Omega_p_Algebra") 
    Gamma_p_Algebra <- mxAlgebra(2 * a_p * hc + 2 * delta_p * ic + a_p * j + 0.5 * v, name="Gamma_p_Algebra") 

    # 6. Scalar Algebra for Offspring Trait
    VY_o_Algebra <- mxAlgebra(2 * delta_o^2 * k + 4 * delta_o^2 * gc + 2 * a_o^2 * j + 4 * a_o^2 * hc + 8 * delta_o * ic * a_o + 2 * delta_o * w + 2 * a_o * v + VF + VE_o, name="VY_o_Algebra")

    # Identification via RDR heritability (Anchoring each trait)
    rdr_left_p  <- mxAlgebra((2*a_p^2*j + 2*delta_p^2*k) * (2*a_p^2*j + 2*delta_p^2*k + VE_p), name="rdr_left_p")
    h2mat_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_parent, name="h2mat_p")
    rdr_right_p <- mxAlgebra(h2mat_p * VY_p, name="rdr_right_p")

    rdr_left_o  <- mxAlgebra((2*a_o^2*j + 2*delta_o^2*k) * (2*a_o^2*j + 2*delta_o^2*k + VE_o), name="rdr_left_o")
    h2mat_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_offspring, name="h2mat_o")
    rdr_right_o <- mxAlgebra(h2mat_o * VY_o, name="rdr_right_o")

    # 7. Transmission and Cross-person Covariances
    gt_Algebra <- mxAlgebra(Omega_td * mu * Omega_td, name="gt_Algebra")
    ht_Algebra <- mxAlgebra(Gamma_td * mu * Gamma_td, name="ht_Algebra")

    # Constraints tying ic, w and v to the equilibrium recursion (as in the same-trait script)
    ic_Algebra <- mxAlgebra(.5 * (Gamma_td * mu * Omega_td + Omega_td * mu * Gamma_td), name="ic_Algebra")
    w_Algebra  <- mxAlgebra(2 * f * Omega_p + 2 * f * tau * mu * Omega_td, name="w_Algebra")
    v_Algebra  <- mxAlgebra(2 * f * Gamma_p + 2 * f * tau * mu * Gamma_td, name="v_Algebra")
    
    thetaNT <- mxAlgebra(2 * delta_o * gc + 2 * a_o * ic + .5 * w, name="thetaNT")
    thetaT  <- mxAlgebra(delta_o * k + thetaNT, name="thetaT")
    
    Yp_Ym   <- mxAlgebra(tau * mu * tau, name="Yp_Ym")
    Yo_Yp   <- mxAlgebra(delta_o * Omega_p + a_o * Gamma_p + f * VY_p + tau * mu * (delta_o * Omega_td + a_o * Gamma_td + f * tau), name = "Yo_Yp")
    Yp_PGSm <- mxAlgebra(tau * mu * Omega_td, name="Yp_PGSm")

    # ---- Cascade AM: latent mating phenotype gamma~ ----
    # Tilde paths = multiplier x the path into Y (ETFD Cascade MVN scripts); see the header for AM_G / AM_E.
    AM_G <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_G_free, values=AM_G_value, label="AMGenMulti", name="AM_G")
    AM_E <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_E_free, values=AM_E_value, label="AMEnvMulti", name="AM_E")
    delta_td <- mxAlgebra(delta_p * AM_G, name="delta_td")
    a_td     <- mxAlgebra(a_p * AM_G, name="a_td")
    # Shortcuts to gamma~ (PDF Part III Section 8): Omega_td = cov(gamma~, [N]T), Gamma_td = cov(gamma~, L[N]T),
    # zeta = cov(gamma~, F), tau = cov(Y, gamma~), Vgamma = var(gamma~)
    Omega_td <- mxAlgebra(2 * delta_td * gc + 2 * a_td * ic + delta_td * k + 0.5 * AM_E * w, name="Omega_td")
    Gamma_td <- mxAlgebra(2 * a_td * hc + 2 * delta_td * ic + a_td * j + 0.5 * AM_E * v, name="Gamma_td")
    zeta     <- mxAlgebra(delta_td * w + a_td * v + AM_E * VF, name="zeta")
    tau      <- mxAlgebra(2 * a_p * Gamma_td + 2 * delta_p * Omega_td + zeta + AM_E * VE_p, name="tau")
    Vgamma   <- mxAlgebra(2 * a_td * Gamma_td + 2 * delta_td * Omega_td + AM_E * zeta + AM_E^2 * VE_p, name="Vgamma")

    # 8. Expected Covariance Matrix (7x7)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_p_Algebra, Yp_Ym,      Yo_Yp,     Omega_p, Omega_p, Yp_PGSm, Yp_PGSm), 
        cbind(Yp_Ym,        VY_p_Algebra, Yo_Yp,   Yp_PGSm, Yp_PGSm, Omega_p, Omega_p), 
        cbind(Yo_Yp,        Yo_Yp,      VY_o_Algebra, thetaT, thetaNT, thetaT,  thetaNT), 
        cbind(Omega_p,      Yp_PGSm,    thetaT,    k+gc,    gc,      gt_Algebra, gt_Algebra), 
        cbind(Omega_p,      Yp_PGSm,    thetaNT,   gc,      k+gc,    gt_Algebra, gt_Algebra), 
        cbind(Yp_PGSm,      Omega_p,    thetaT,    gt_Algebra, gt_Algebra, k+gc,    gc),      
        cbind(Yp_PGSm,      Omega_p,    thetaNT,   gt_Algebra, gt_Algebra, gc,      k+gc)),   
        dimnames = list(colnames(Example_Data), colnames(Example_Data)), name="expCov")

    Means <- mxMatrix(type = "Full", nrow = 1, ncol = 7, free = TRUE, values = 0, 
        label = c("meanYp1", "meanYm1", "meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"), 
        dimnames = list(NULL, colnames(Example_Data)), name = "expMeans")

    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans", 
                                            dimnames=colnames(Example_Data))
    
    # 10. Constraints (Note: Cross-generational VY equality removed)
    Constraints <- list(
        mxConstraint(VY_p == VY_p_Algebra, name="VY_p_con"),
        mxConstraint(VY_o == VY_o_Algebra, name="VY_o_con"),
        mxConstraint(Omega_p == Omega_p_Algebra, name="Om_p_con"),
        mxConstraint(Gamma_p == Gamma_p_Algebra, name="Ga_p_con"),
        mxConstraint(gc == gt_Algebra, name="gc_eq"),
        mxConstraint(hc == ht_Algebra, name="hc_eq"),
        mxConstraint(ic == ic_Algebra, name="ic_eq"),
        mxConstraint(w == w_Algebra, name="w_eq"),
        mxConstraint(v == v_Algebra, name="v_eq"),
        mxConstraint(rdr_left_p == rdr_right_p, name="rdr_p_con"),
        mxConstraint(rdr_left_o == rdr_right_o, name="rdr_o_con")
    )

    Params <- list(
        AM_G, AM_E, delta_td, a_td, Omega_td, Gamma_td, zeta, tau, Vgamma, VF, VF_Constraint,
        VY_p, VY_o, VE_p, VE_o, delta_p, a_p, delta_o, a_o, k, j, Omega_p, Gamma_p, mu, ic, gc, hc, f, w, v,
        VY_p_Algebra, VF_p_Algebra, VY_o_Algebra, Omega_p_Algebra, Gamma_p_Algebra,
        gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
        h2mat_p, h2mat_o, rdr_left_p, rdr_right_p, rdr_left_o, rdr_right_o,
        thetaNT, thetaT, Yp_Ym, Yo_Yp, Yp_PGSm,
        CovMatrix, Means, ModelExpectations, mxFitFunctionML(), Constraints
    )

    Model1 <- mxModel("Cascade_UniSEM_DiffTrait_Observed_Corrected", Params, mxData(observed=Example_Data, type="raw"))
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1), 
                          intervals=T, silent=T, exhaustive = exhaustive, 
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)
    
    return(summary(fitModel1, verbose = TRUE))
}
