## Univariate SEM-PGS Cascade model: different traits in parents and offspring, parental phenotypes observed.
## The parental trait drives assortative mating and vertical transmission; the offspring trait is the outcome.
## The parents' latent genetic effect a_p is estimated from the parental phenotypes.
## Equilibrium: assortative mating and vertical transmission have gone on for many generations.
## Spouses assort on a latent mating phenotype gamma~ = AM_G*(delta*(T + NT) + a*(LT + LNT)) + AM_E*(F + E).
##
## Estimates (returned by summary(fit, verbose = TRUE)):
##   delta_p, a_p, delta_o, a_o  PGS and latent-genetic-score effects on the parental / offspring trait
##   f                           vertical transmission from each parent's phenotype to the offspring's environment F
##   mu                          copath between spouses' mating phenotypes gamma~
##   AM_G, AM_E                  mating-mechanism multipliers, when freed (AMGenMulti, AMEnvMulti)
##   VY_p, VY_o, VE_p, VE_o      phenotypic and residual variances of the parental and offspring trait
##   plus Omega_p, Gamma_p, gc, hc, ic, w, v, VF: covariance terms held to their model values by constraints.
##
## Input:
##   data_path   text file with a header row (read by data.table::fread) and the columns
##                 Yp1, Ym1    father's and mother's phenotype on the parental trait
##                 Yo1         offspring phenotype on the offspring trait
##                 Tp1, NTp1   father's transmitted and non-transmitted haplotypic PGS
##                 Tm1, NTm1   mother's transmitted and non-transmitted haplotypic PGS
##               The file must hold exactly these seven columns, in the order Yp1, Ym1, Yo1, Tp1, NTp1, Tm1, NTm1.
##               Each haplotypic PGS is assumed to have variance k = .5 in the base population.
##   h2_RDR_offspring  RDR heritability of the offspring trait; identifies a_o
##   AM_G_value, AM_G_free, AM_E_value, AM_E_free  multipliers on the genetic (AM_G) and the F + E (AM_E) paths into
##               gamma~: (1, 1) primary phenotypic AM, (1, 0) genetic homogamy, (0, 1) social homogamy. At least one
##               must be fixed. Default: AM_E fixed at 1, AM_G estimated (for genetic homogamy fix AM_G = 1, free AM_E).
##   feaTol, optTol, extraTries, exhaustive, jitterMean, jitterVar  optimizer settings (NPSOL, mxTryHard)

fitUniSEMPGS_Cascade_DiffTrait_ObservedYPYM_EstimatedAParent <- function(data_path, h2_RDR_offspring, AM_G_value = .5, AM_G_free = TRUE, AM_E_value = 1, AM_E_free = FALSE, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
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

    # 1. Phenotypic and Residual Variances (Separated)
    VY_p  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=2,   label="VY_parent", name="VY_p", lbound = .001) 
    VY_o  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=2,   label="VY_offspring", name="VY_o", lbound = .001) 
    VE_p  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.5,  label="VE_parent", name="VE_p", lbound = .001) 
    VE_o  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.5,  label="VE_offspring", name="VE_o", lbound = .001) 

    # 2. Genetic effects (Generation Specific)
    delta_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="delta_p11", name="delta_p") 
    a_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.8, label="a_p11", name="a_p", lbound = .001)     
    
    delta_o <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="delta_o11", name="delta_o") 
    a_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.8, label="a_o11", name="a_o", lbound = .001)     
    
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
    Omega_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.6, label="Omega_p11", name="Omega_p") 
    Gamma_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.5, label="Gamma_p11", name="Gamma_p") 
    
    mu <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="mu11", name="mu") 
    ic <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.07, label="ic11", name="ic") 
    gc <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.14, label="gc11", name="gc")  
    hc <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.04, label="hc11", name="hc")  

    # 4. Vertical Transmission
    f <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.16, label="f11", name="f") 
    w <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3,  label="w11", name="w") 
    v <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2,  label="v11", name="v") 

    # 5. Scalar Algebra for Parent Trait (Equilibrium)
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

    # Identification of a_o via RDR
    rdr_left_o  <- mxAlgebra((2*a_o^2*j + 2*delta_o^2*k) * (2*a_o^2*j + 2*delta_o^2*k + VE_o), name="rdr_left_o")
    h2mat_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_offspring, name="h2mat_o")
    rdr_right_o <- mxAlgebra(h2mat_o * VY_o, name="rdr_right_o")

    # 7. Assortment and Transmission Algebra
    gt_Algebra <- mxAlgebra(Omega_td * mu * Omega_td, name="gt_Algebra")
    ht_Algebra <- mxAlgebra(Gamma_td * mu * Gamma_td, name="ht_Algebra")

    # Constraints tying ic, w and v to the equilibrium recursion (as in the same-trait script)
    ic_Algebra <- mxAlgebra(.5 * (Gamma_td * mu * Omega_td + Omega_td * mu * Gamma_td), name="ic_Algebra")
    w_Algebra  <- mxAlgebra(2 * f * Omega_p + 2 * f * tau * mu * Omega_td, name="w_Algebra")
    v_Algebra  <- mxAlgebra(2 * f * Gamma_p + 2 * f * tau * mu * Gamma_td, name="v_Algebra")
    
    thetaNT <- mxAlgebra(2 * delta_o * gc + 2 * a_o * ic + .5 * w, name="thetaNT")
    thetaT  <- mxAlgebra(delta_o * k + thetaNT, name="thetaT")
    
    # Intergenerational Phenotypic Covariances (Between Parent Trait and Child Trait)
    Yp_Ym   <- mxAlgebra(tau * mu * tau, name="Yp_Ym")
    Yo_Yp   <- mxAlgebra(delta_o * Omega_p + a_o * Gamma_p + f * VY_p + tau * mu * (delta_o * Omega_td + a_o * Gamma_td + f * tau), name = "Yo_Yp")
    
    # Cross-person PGS-Phenotype Covariances
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
        cbind(VY_p_Algebra, Yp_Ym,      Yo_Yp,     Omega_p, Omega_p, Yp_PGSm, Yp_PGSm), # Yp
        cbind(Yp_Ym,        VY_p_Algebra, Yo_Yp,   Yp_PGSm, Yp_PGSm, Omega_p, Omega_p), # Ym
        cbind(Yo_Yp,        Yo_Yp,      VY_o_Algebra, thetaT, thetaNT, thetaT,  thetaNT), # Yo
        cbind(Omega_p,      Yp_PGSm,    thetaT,    k+gc,    gc,      gt_Algebra, gt_Algebra), # Tp
        cbind(Omega_p,      Yp_PGSm,    thetaNT,   gc,      k+gc,    gt_Algebra, gt_Algebra), # NTp
        cbind(Yp_PGSm,      Omega_p,    thetaT,    gt_Algebra, gt_Algebra, k+gc,    gc), # Tm
        cbind(Yp_PGSm,      Omega_p,    thetaNT,   gt_Algebra, gt_Algebra, gc,      k+gc)), # NTm
        dimnames=list(colnames(Example_Data), colnames(Example_Data)), name="expCov")

    # 9. Means and Expectations
    Means <- mxMatrix(type = "Full", nrow = 1, ncol = 7, free = TRUE, values = 0, 
        label = c("meanYp1", "meanYm1", "meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"), 
        dimnames = list(NULL, colnames(Example_Data)), name = "expMeans")

    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans", dimnames=colnames(Example_Data))
    
    # 10. Constraints
    Constraints <- list(
        mxConstraint(VY_p == VY_p_Algebra, name="VY_p_eq"),
        mxConstraint(VY_o == VY_o_Algebra, name="VY_o_eq"),
        mxConstraint(Omega_p == Omega_p_Algebra, name="Om_p_eq"),
        mxConstraint(Gamma_p == Gamma_p_Algebra, name="Ga_p_eq"),
        mxConstraint(gc == gt_Algebra, name="gc_eq"),
        mxConstraint(hc == ht_Algebra, name="hc_eq"),
        mxConstraint(ic == ic_Algebra, name="ic_eq"),
        mxConstraint(w == w_Algebra, name="w_eq"),
        mxConstraint(v == v_Algebra, name="v_eq"),
        mxConstraint(rdr_left_o == rdr_right_o, name="rdr_o_con")
    )

    Params <- list(
        AM_G, AM_E, delta_td, a_td, Omega_td, Gamma_td, zeta, tau, Vgamma, VF, VF_Constraint,
        VY_p, VY_o, VE_p, VE_o, delta_p, a_p, delta_o, a_o, k, j, Omega_p, Gamma_p, mu, ic, gc, hc, f, w, v,
        VY_p_Algebra, VF_p_Algebra, VY_o_Algebra, Omega_p_Algebra, Gamma_p_Algebra,
        gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
        h2mat_o, rdr_left_o, rdr_right_o, thetaNT, thetaT, Yp_Ym, Yo_Yp, Yp_PGSm,
        CovMatrix, Means, ModelExpectations, mxFitFunctionML(), Constraints
    )

    Model1 <- mxModel("Cascade_UniSEM_DiffTrait_ObservedParents_RDR", Params, mxData(observed=Example_Data, type="raw"))
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1), 
                          intervals=T, silent=T, exhaustive = exhaustive, 
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)
    
    return(summary(fitModel1, verbose = TRUE))
}
