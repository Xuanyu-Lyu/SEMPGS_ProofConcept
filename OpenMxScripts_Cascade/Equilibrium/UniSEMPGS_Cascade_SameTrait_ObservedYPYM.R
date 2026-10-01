## SEM-PGS CASCADE version of OpenMxScripts/Equilibrium/UniSEMPGS_SameTrait_ObservedYPYM.R: the same model, except that spouses
## assort on the latent mating phenotype gamma~ = delta~*(T+NT) + a~*(LT+LNT) + 1~*F + 1~*E instead of on Y
## (SEM_PGS_Cascade_model.pdf, Part III). Following the ETFD Cascade MVN scripts, each path into
## gamma~ is the path into Y times a multiplier: AM_G on the genetic paths (delta~ = AM_G*delta, a~ = AM_G*a)
## and AM_E on the non-genetic paths (the PDF's 1~ on F and E). AM_G = AM_E = 1: primary phenotypic AM;
## AM_E = 0: genetic homogamy; AM_G = 0: social homogamy. gamma~ has no scale of its own, so at least one
## multiplier must be fixed (default AM_E = 1, the analogue of the MVN script's fixed AM_U); mu is the
## copath between spouses' gamma~. By default AM_G is estimated; for genetic homogamy fix AM_G = 1 and free AM_E instead.
##
## This script is a function that fits a version univariate SEM-PGS where parent and offspring have the same trait and the parental phenotypes are observed.

fitUniSEMPGS_Cascade_SameTrait_ObservedYPYM <- function(data_path, AM_G_value = .5, AM_G_free = TRUE, AM_E_value = 1, AM_E_free = FALSE, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
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

    # Phenotypic and Residual Variance
    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=2, label="VY11", name="VY", lbound = .001) 
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.5, label="VE11", name="VE", lbound = .001) 

    # Scalar Algebra for Variances
    VY_Algebra <- mxAlgebra(2 * delta * Omega + 2 * a * Gamma + w * delta + v * a + VF + VE, name="VY_Algebra")
    VF_Algebra <- mxAlgebra(2 * f^2 * VY + 2 * f^2 * tau^2 * mu, name="VF_Algebra")
    # Cascade: VF = 2f^2(VY + tau^2*mu) depends on tau = cov(Y, gamma~), which itself depends on VF,
    # so VF is a free parameter held to its algebra by a constraint (as VY is).
    VF    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="VF11", name="VF")
    VF_Constraint <- mxConstraint(VF == VF_Algebra, name='VF_Constraint')

    VY_Constraint    <- mxConstraint(VY == VY_Algebra, name='VY_Constraint')

    # Genetic effects
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="delta11", name="delta") 
    a     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.8, label="a11", name="a", lbound = .001)     
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")     
    j     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")     
    # k and j are fixed at .5, which assumes the PGS (and, by convention, the latent LGS) is standardised
    # in the base population. Empirical PGS are usually scaled in the analysed sample, i.e. at equilibrium,
    # where the haplotypic PGS variance is k + g and the full PGS variance is 2k + 4g (g = gc = gt, h = hc = ht).
    # To use one of those scalings, replace the two fixed matrices above with the matching algebra pair
    # (and drop j_Algebra / j_constraint below, which would otherwise force gc == hc):
    #   haplotypic PGS scaled to variance 1/2 at equilibrium:  k = 1/2 - g,   j = 1/2 - h
    # k <- mxAlgebra(.5 - gc, name = "k")
    # j <- mxAlgebra(.5 - hc, name = "j")
    #   full PGS standardised to variance 1 at equilibrium:    k = 1/2 - 2g,  j = 1/2 - 2h
    # k <- mxAlgebra(.5 - 2 * gc, name = "k")
    # j <- mxAlgebra(.5 - 2 * hc, name = "j")
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.6, label="Omega11", name="Omega") 
    Gamma <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.5, label="Gamma11", name="Gamma") 

    Omega_Algebra <- mxAlgebra(2 * delta * gc + 2 * a * ic + delta * k + 0.5 * w , name="Omega_Algebra") 
    Gamma_Algebra <- mxAlgebra(2 * a * hc + 2 * delta * ic + a * j + 0.5 * v, name="Gamma_Algebra") 

    Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
    Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')
    
    adelta_Constraint_Algebra <- mxAlgebra(delta, name = "adelta_Constraint_Algebra")
    adelta_Constraint <- mxConstraint(a == delta, name = "adelta_Constraint")

    j_Algebra    <- mxAlgebra(k, name = "j_Algebra")
    j_constraint <- mxConstraint(j == j_Algebra, name = "j_constraint")

    # Assortative mating effects
    mu    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="mu11", name="mu") 
    gt    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2,  label="gt11", name="gt")  
    ht    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.05, label="ht11", name="ht")  
    gc    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.14, label="gc11", name="gc")  
    hc    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=0.04,label="hc11", name="hc")  
    
    gt_Algebra <- mxAlgebra(Omega_td * mu * Omega_td, name="gt_Algebra") 
    ht_Algebra <- mxAlgebra(Gamma_td * mu * Gamma_td, name="ht_Algebra") 
    gc_Algebra <- mxAlgebra(gt, name="gc_Algebra") # Simplified for scalar
    hc_Algebra <- mxAlgebra(ht, name="hc_Algebra") 

    ic    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.07, label="ic11",   name="ic") 
    
    ic_Algebra   <- mxAlgebra(Omega_td * mu * Gamma_td, name="ic_Algebra") 
    
    gt_constraint   <- mxConstraint(gt == gt_Algebra, name='gt_constraint')
    ht_constraint   <- mxConstraint(ht == ht_Algebra, name='ht_constraint')
    gc_constraint   <- mxConstraint(gc == gc_Algebra, name='gc_constraint')
    hc_constraint   <- mxConstraint(hc == hc_Algebra, name='hc_constraint')
    ic_constraint   <- mxConstraint(ic == ic_Algebra, name='ic_constraint')

    # Vertical transmission effects
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.16, label="f11", name="f") 
    w     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3,  label="w11", name="w") 
    v     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2,  label="v11", name="v") 
    
    w_Algebra     <- mxAlgebra(2 * f * Omega + 2 * f * tau * mu * Omega_td, name="w_Algebra")    
    v_Algebra     <- mxAlgebra(2 * f * Gamma + 2 * f * tau * mu * Gamma_td, name="v_Algebra")    
    wv_constraint_algebra <- mxAlgebra((w * sqrt(2 * delta^2 * k) / sqrt(2 * a^2 * j)), name='wv_constraint_algebra')

    v_constraint  <- mxConstraint(v == v_Algebra, name='v_constraint')
    w_constraint  <- mxConstraint(w == w_Algebra, name='w_constraint')
    wv_constraint <- mxConstraint(v == wv_constraint_algebra, name='wv_constraint')

    # Between-people covariances
    thetaNT <- mxAlgebra(2 * delta * gc + 2 * a * ic + .5 * w, name="thetaNT")
    thetaT  <- mxAlgebra(delta * k + thetaNT, name="thetaT")
    Yp_PGSm <- mxAlgebra(tau * mu * Omega_td, name="Yp_PGSm")
    Ym_PGSp <- mxAlgebra(tau * mu * Omega_td, name="Ym_PGSp")
    Yp_Ym   <- mxAlgebra(tau * mu * tau,    name="Yp_Ym")
    Ym_Yp   <- mxAlgebra(tau * mu * tau,    name="Ym_Yp")
    Yo_Yp   <- mxAlgebra(delta * Omega + a * Gamma + f * VY + tau * mu * (delta * Omega_td + a * Gamma_td + f * tau), name = "Yo_Yp")
    Yo_Ym   <- mxAlgebra(delta * Omega + a * Gamma + f * VY + tau * mu * (delta * Omega_td + a * Gamma_td + f * tau), name = "Yo_Ym")

    # ---- Cascade AM: latent mating phenotype gamma~ ----
    # Tilde paths = multiplier x the path into Y (ETFD Cascade MVN scripts); see the header for AM_G / AM_E.
    AM_G <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_G_free, values=AM_G_value, label="AMGenMulti", name="AM_G")
    AM_E <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_E_free, values=AM_E_value, label="AMEnvMulti", name="AM_E")
    delta_td <- mxAlgebra(delta * AM_G, name="delta_td")
    a_td     <- mxAlgebra(a * AM_G, name="a_td")
    # Shortcuts to gamma~ (PDF Part III Section 8): Omega_td = cov(gamma~, [N]T), Gamma_td = cov(gamma~, L[N]T),
    # zeta = cov(gamma~, F), tau = cov(Y, gamma~), Vgamma = var(gamma~)
    Omega_td <- mxAlgebra(2 * delta_td * gc + 2 * a_td * ic + delta_td * k + 0.5 * AM_E * w, name="Omega_td")
    Gamma_td <- mxAlgebra(2 * a_td * hc + 2 * delta_td * ic + a_td * j + 0.5 * AM_E * v, name="Gamma_td")
    zeta     <- mxAlgebra(delta_td * w + a_td * v + AM_E * VF, name="zeta")
    tau      <- mxAlgebra(2 * a * Gamma_td + 2 * delta * Omega_td + zeta + AM_E * VE, name="tau")
    Vgamma   <- mxAlgebra(2 * a_td * Gamma_td + 2 * delta_td * Omega_td + AM_E * zeta + AM_E^2 * VE, name="Vgamma")

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

    Means <- mxMatrix(type = "Full", nrow = 1, ncol = 7, free = TRUE, values = 0, 
        label = c("meanYp1", "meanYm1", "meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"), 
        dimnames = list(NULL, c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1")),
        name = "expMeans")

    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans", dimnames=c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1"))
    Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )
    FitFunctionML   <- mxFitFunctionML()

    Params <- list(
                AM_G, AM_E, delta_td, a_td, Omega_td, Gamma_td, zeta, tau, Vgamma, VF, VF_Constraint,
                VY, VE, delta, a, k, j, Omega, Gamma, mu, gt, ht, gc, hc, ic, f, w, v,
                VY_Algebra, VF_Algebra, Omega_Algebra, Gamma_Algebra, adelta_Constraint_Algebra, j_Algebra, gt_Algebra, ht_Algebra, gc_Algebra, hc_Algebra, ic_Algebra, w_Algebra, v_Algebra, wv_constraint_algebra,
                VY_Constraint, Gamma_Constraint, j_constraint, ht_constraint, hc_constraint, ic_constraint, v_constraint, w_constraint,
                thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Ym_Yp, Yo_Yp, Yo_Ym, 
                CovMatrix, Means, ModelExpectations, FitFunctionML)

    Model1 <- mxModel("Cascade_UniSEM_Scalar", Params, Example_Data_Mx)
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1), intervals=T, silent=T, exhaustive = exhaustive, jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)
    
    return(summary(fitModel1, verbose = TRUE))
}
