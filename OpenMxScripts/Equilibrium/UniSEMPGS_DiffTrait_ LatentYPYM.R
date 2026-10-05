## Univariate SEM-PGS model: different traits in parents and offspring, parental phenotypes not observed.
## The parental trait drives assortative mating and vertical transmission; the offspring trait is the outcome.
## Equilibrium: assortative mating and vertical transmission have gone on for many generations.
##
## Estimates (returned by summary(fit, verbose = TRUE)):
##   delta_p, a_p, delta_o, a_o  PGS and latent-genetic-score effects on the parental / offspring trait
##   f                           vertical transmission from each parent's phenotype to the offspring's environment F
##   mu                          assortative-mating copath, cov(Yp, Ym) = mu*VY^2; fixed at mu_fixed (estimated if NULL)
##   VY, VE_p, VE_o              shared phenotypic variance and the two traits' residual variances
##   plus Omega_p, Gamma_p, gc, hc, ic, w, v: covariance terms held to their model values by constraints.
##
## Input:
##   data_path   text file with a header row (read by data.table::fread) and the columns
##                 Yo1         offspring phenotype on the offspring trait
##                 Tp1, NTp1   father's transmitted and non-transmitted haplotypic PGS
##                 Tm1, NTm1   mother's transmitted and non-transmitted haplotypic PGS
##               The file should hold only these columns.
##               Each haplotypic PGS is assumed to have variance k = .5 in the base population.
##   h2_RDR_parent, h2_RDR_offspring  RDR heritability of the parental and of the offspring trait; identify a_p, a_o
##   mu_fixed    known assortative-mating copath of the parental trait. Needed for identification: with NULL, mu is
##               estimated, but the model is under-identified and can converge to a wrong point.
##   feaTol, optTol, extraTries, exhaustive, jitterMean, jitterVar  optimizer settings (NPSOL, mxTryHard)

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
    Omega_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="Omega_p11", name="Omega_p", lbound = .001) 
    Gamma_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="Gamma_p11", name="Gamma_p") 
    
    # 'mu' is free unless an external value is supplied via mu_fixed (see mu_fixed in the header).
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