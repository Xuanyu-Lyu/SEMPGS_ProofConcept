## Univariate SEM-PGS model: different traits in parents and offspring, parental phenotypes observed.
## The parental trait drives assortative mating and vertical transmission; the offspring trait is the outcome.
## RDR heritability identifies the latent genetic effect of both traits (a_p and a_o).
## Disequilibrium: vertical transmission is at equilibrium; the parents are the first assortatively mated generation.
##
## Estimates (returned by summary(fit, verbose = TRUE)):
##   delta_p, a_p, delta_o, a_o  PGS and latent-genetic-score effects on the parental / offspring trait
##   f                           vertical transmission from each parent's phenotype to the offspring's environment F
##   mu                          assortative-mating copath, cov(Yp, Ym) = mu*VY^2
##   VY_p, VY_o, VE_p, VE_o      phenotypic and residual variances of the parental and offspring trait
##   plus Omega_p, Gamma_p: covariance terms held to their model values by constraints.
##
## Input:
##   data_path   text file with a header row (read by data.table::fread) and the columns
##                 Yp1, Ym1    father's and mother's phenotype on the parental trait
##                 Yo1         offspring phenotype on the offspring trait
##                 Tp1, NTp1   father's transmitted and non-transmitted haplotypic PGS
##                 Tm1, NTm1   mother's transmitted and non-transmitted haplotypic PGS
##               The file must hold exactly these seven columns, in the order Yp1, Ym1, Yo1, Tp1, NTp1, Tm1, NTm1.
##               Each haplotypic PGS is assumed to have variance k = .5 in the base population.
##   h2_RDR_parent, h2_RDR_offspring  RDR heritability of the parental and of the offspring trait; identify a_p, a_o
##   feaTol, optTol, extraTries, exhaustive, jitterMean, jitterVar  optimizer settings (NPSOL, mxTryHard)

fitUniSEMPGS_DisEq_DiffTrait_ObservedYPYM_FixedAParent <- function(data_path, h2_RDR_parent, h2_RDR_offspring, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
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
    delta_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.25, label="delta_p11", name="delta_p")
    a_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.7, label="a_p11", name="a_p", lbound = .001)

    delta_o <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="delta_o11", name="delta_o")
    a_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.65, label="a_o11", name="a_o", lbound = .001)

    k <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    j <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")
    # k and j are fixed at .5, which assumes the PGS (and, by convention, the latent LGS) is standardised
    # in the base population. Under disequilibrium there is no within-person haplotype covariance
    # (gc = hc = 0), so the haplotypic PGS variance is k in every generation, but the offspring PGS
    # variance is 2k + 2gt because its two haplotypes come from assorted parents. If the PGS is instead
    # standardised to variance 1 in the offspring sample, replace the two fixed matrices above with
    # k <- mxAlgebra(.5 - gt, name = "k")
    # j <- mxAlgebra(.5 - ht, name = "j")
    # (scaling the haplotypic PGS to variance 1/2 leaves k = j = .5 unchanged here).

    # 3. Covariances and Assortment
    Omega_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="Omega_p11", name="Omega_p")
    Gamma_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="Gamma_p11", name="Gamma_p")

    mu <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="mu11", name="mu")

    # gc = hc = 0 (No AM in previous generation)
    gc <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=0, label="gc11", name="gc")
    hc <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=0, label="hc11", name="hc")

    # 4. Vertical Transmission
    f <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="f11", name="f")

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

    # 6. Offspring-trait variance: the offspring's own VY (w_o, v_o, VF_o: see the offspring-generation block below)
    VY_o_Algebra <- mxAlgebra(2 * delta_o^2 * k + 2 * delta_o^2 * gt + 2 * a_o^2 * j + 2 * a_o^2 * ht + 4 * delta_o * a_o * ic + 2 * delta_o * w_o + 2 * a_o * v_o + VF_o + VE_o, name="VY_o_Algebra")

    # Identification via RDR heritability (Anchoring each trait)
    rdr_left_p  <- mxAlgebra((2*a_p^2*j + 2*delta_p^2*k) * (2*a_p^2*j + 2*delta_p^2*k + VE_p), name="rdr_left_p")
    h2mat_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_parent, name="h2mat_p")
    rdr_right_p <- mxAlgebra(h2mat_p * VY_p, name="rdr_right_p")

    rdr_left_o  <- mxAlgebra((2*a_o^2*j + 2*delta_o^2*k) * (2*a_o^2*j + 2*delta_o^2*k + VE_o), name="rdr_left_o")
    h2mat_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_offspring, name="h2mat_o")
    rdr_right_o <- mxAlgebra(h2mat_o * VY_o, name="rdr_right_o")

    # 7. Transmission and Cross-person Covariances (functions of mu and current Omega_p/Gamma_p).
    # (itlo and itol -- the transmitted-vs-non-transmitted haplotype covariances -- collapse to a
    # single quantity 'ic' here because Gamma_p*mu*Omega_p == Omega_p*mu*Gamma_p for scalars; the
    # itlo/itol split only matters once Omega_p and Gamma_p become asymmetric bivariate matrices.)
    gt_Algebra <- mxAlgebra(Omega_p * mu * Omega_p, name="gt_Algebra")
    ht_Algebra <- mxAlgebra(Gamma_p * mu * Gamma_p, name="ht_Algebra")
    ic_Algebra <- mxAlgebra(Omega_p * mu * Gamma_p, name="ic_Algebra")

    gt <- mxAlgebra(gt_Algebra, name="gt")
    ht <- mxAlgebra(ht_Algebra, name="ht")
    ic <- mxAlgebra(ic_Algebra, name="ic")

    # ---- Offspring generation ----
    # The offspring's F comes from AM-mated parents, so its covariance with the parents' haplotypes (w_o, v_o)
    # and its variance (VF_o) carry the cross-mate term mu*VY. (w, v and the VF algebra above belong to the
    # parents, whose own parents mated at random.)
    w_o  <- mxAlgebra(2 * f * Omega_p + 2 * f * VY_p * mu * Omega_p, name="w_o")
    v_o  <- mxAlgebra(2 * f * Gamma_p + 2 * f * VY_p * mu * Gamma_p, name="v_o")
    VF_o <- mxAlgebra(2 * f^2 * VY_p + 2 * f^2 * VY_p^2 * mu, name="VF_o")
    # thetaNT / thetaT: covariance of Yo with ONE parental haplotype, i.e. with one data column (NTp1, Tp1, ...).
    # NOTE: in the papers (Balbona et al. 2021; Lyu et al.) theta is the SUM over the father's and the mother's
    # haplotype, thetaNT = cov(Yo, NTp + NTm) = 2*delta_o*gt + 2*a_o*ic + w_o: twice these.
    # The parents' own two haplotypes are uncorrelated (gc = hc = 0: first AM event), so per haplotype
    #   cov(Yo, NTp) = delta_o*cov(Tm, NTp) + a_o*cov(LTm, NTp) + cov(Fo, NTp) = delta_o*gt + a_o*ic + w_o/2
    # (ic here is the across-parent LGS-PGS covariance Omega*mu*Gamma)
    thetaNT <- mxAlgebra(delta_o * gt + a_o * ic + 0.5 * w_o, name="thetaNT")
    # cov(Yo, Tp) adds delta_o*var(Tp) = delta_o*k
    thetaT  <- mxAlgebra(delta_o * k + thetaNT, name="thetaT")

    Yp_Ym   <- mxAlgebra(VY_p * mu * VY_p, name="Yp_Ym")
    Yo_Yp   <- mxAlgebra((delta_o * Omega_p + a_o * Gamma_p + f * VY_p) * (1 + mu * VY_p), name = "Yo_Yp")
    Yp_PGSm <- mxAlgebra(VY_p * mu * Omega_p, name="Yp_PGSm")

    # 8. Expected Covariance Matrix (7x7)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_p_Algebra, Yp_Ym,      Yo_Yp,     Omega_p, Omega_p, Yp_PGSm, Yp_PGSm),
        cbind(Yp_Ym,        VY_p_Algebra, Yo_Yp,   Yp_PGSm, Yp_PGSm, Omega_p, Omega_p),
        cbind(Yo_Yp,        Yo_Yp,      VY_o_Algebra, thetaT, thetaNT, thetaT,  thetaNT),
        cbind(Omega_p,      Yp_PGSm,    thetaT,    k+gc,    gc,      gt,      gt),
        cbind(Omega_p,      Yp_PGSm,    thetaNT,   gc,      k+gc,    gt,      gt),
        cbind(Yp_PGSm,      Omega_p,    thetaT,    gt,      gt,      k+gc,    gc),
        cbind(Yp_PGSm,      Omega_p,    thetaNT,   gt,      gt,      gc,      k+gc)),
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
        mxConstraint(rdr_left_p == rdr_right_p, name="rdr_p_con"),
        mxConstraint(rdr_left_o == rdr_right_o, name="rdr_o_con")
    )

    Params <- list(
        w_o, v_o, VF_o,
        VY_p, VY_o, VE_p, VE_o, delta_p, a_p, delta_o, a_o, k, j, Omega_p, Gamma_p, mu, gc, hc, f,
        VY_p_Algebra, VF_p_Algebra, VY_o_Algebra, Omega_p_Algebra, Gamma_p_Algebra,
        gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
        gt, ht, ic, w, v,
        h2mat_p, h2mat_o, rdr_left_p, rdr_right_p, rdr_left_o, rdr_right_o,
        thetaNT, thetaT, Yp_Ym, Yo_Yp, Yp_PGSm,
        CovMatrix, Means, ModelExpectations, mxFitFunctionML(), Constraints
    )

    Model1 <- mxModel("UniSEM_NonEquilibrium_DiffTrait_Observed_FixedAParent", Params, mxData(observed=Example_Data, type="raw"))
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
