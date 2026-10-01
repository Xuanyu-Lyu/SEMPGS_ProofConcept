## SEM-PGS CASCADE version of OpenMxScripts/Disequilibrium/UniSEMPGS_DisEq_DiffTrait_LatentYPYM.R: the same model, except that spouses
## assort on the latent mating phenotype gamma~ = delta~*(T+NT) + a~*(LT+LNT) + 1~*F + 1~*E instead of on Y
## (SEM_PGS_Cascade_model.pdf, Part IV). Following the ETFD Cascade MVN scripts, each path into
## gamma~ is the path into Y times a multiplier: AM_G on the genetic paths (delta~ = AM_G*delta, a~ = AM_G*a)
## and AM_E on the non-genetic paths (the PDF's 1~ on F and E). AM_G = AM_E = 1: primary phenotypic AM;
## AM_E = 0: genetic homogamy; AM_G = 0: social homogamy. gamma~ has no scale of its own, so at least one
## multiplier must be fixed (default AM_E = 1, the analogue of the MVN script's fixed AM_U); mu is the
## copath between spouses' gamma~. Both multipliers are fixed by default: with latent parental phenotypes the
## mating mechanism is not identified, so it must be supplied (the hypothesized mechanism under test).
## Offspring side, derived from Yo = delta*(Tp+Tm) + a*(LTp+LTm) + Fo + Eo (see ../PDF_Errata.md): the offspring's
## haplotypes covary across parents (gt, ht, ic), their F comes from AM-mated parents (w_o, v_o, VF_o carry the
## tau*mu terms), they have their own phenotypic variance VY_o, and thetaT/thetaNT are per single-haplotype
## column (the papers' theta is the sum over the father's and the mother's haplotype). Yo_Yp follows PDF
## Part IV Section 19 (it includes f*tau^2*mu = f*cov(Ym, Yp), which the original DisEq Yo_Yp omits).
##
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

fitUniSEMPGS_Cascade_DisEq_DiffTrait_LatentYPYM <- function(data_path, h2_RDR_parent, h2_RDR_offspring, mu_fixed = NULL, AM_G_value = 1, AM_G_free = FALSE, AM_E_value = 1, AM_E_free = FALSE, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
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
    # k and j are fixed at .5, which assumes the PGS (and, by convention, the latent LGS) is standardised
    # in the base population. Under disequilibrium there is no within-person haplotype covariance
    # (gc = hc = 0), so the haplotypic PGS variance is k in every generation, but the offspring PGS
    # variance is 2k + 2gt because its two haplotypes come from assorted parents. If the PGS is instead
    # standardised to variance 1 in the offspring sample, replace the two fixed matrices above with
    # k <- mxAlgebra(.5 - gt, name = "k")
    # j <- mxAlgebra(.5 - ht, name = "j")
    # (scaling the haplotypic PGS to variance 1/2 leaves k = j = .5 unchanged here).

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

    # Offspring-trait variance: the offspring's own VY (w_o, v_o, VF_o: see the offspring-generation block below)
    VY_o_Algebra <- mxAlgebra(2 * delta_o^2 * k + 2 * delta_o^2 * gt + 2 * a_o^2 * j + 2 * a_o^2 * ht + 4 * delta_o * a_o * ic + 2 * delta_o * w_o + 2 * a_o * v_o + VF_o + VE_o, name="VY_o_Algebra")

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
    gt_Algebra <- mxAlgebra(Omega_td * mu * Omega_td, name="gt_Algebra")
    ht_Algebra <- mxAlgebra(Gamma_td * mu * Gamma_td, name="ht_Algebra")
    ic_Algebra <- mxAlgebra(Omega_td * mu * Gamma_td, name="ic_Algebra")

    gt <- mxAlgebra(gt_Algebra, name="gt")
    ht <- mxAlgebra(ht_Algebra, name="ht")
    ic <- mxAlgebra(ic_Algebra, name="ic")

    # 7. Offspring-PGS Covariances
    # ---- Offspring generation ----
    # The offspring's F comes from AM-mated parents, so its covariance with the parents' haplotypes (w_o, v_o)
    # and its variance (VF_o) carry the cross-mate tau*mu terms (PDF Part IV Section 18: w2, v2, VF2).
    # (w, v and the VF algebra above belong to the parents, whose own parents mated at random.)
    w_o  <- mxAlgebra(2 * f * Omega_p + 2 * f * tau * mu * Omega_td, name="w_o")
    v_o  <- mxAlgebra(2 * f * Gamma_p + 2 * f * tau * mu * Gamma_td, name="v_o")
    VF_o <- mxAlgebra(2 * f^2 * VY + 2 * f^2 * tau^2 * mu, name="VF_o")
    # thetaNT / thetaT: covariance of Yo with ONE parental haplotype, i.e. with one data column (NTp1, Tp1, ...).
    # NOTE: in the papers (Balbona et al. 2021, Lyu et al., SEM_PGS_Cascade_model.pdf) theta is the SUM over the
    # father's and the mother's haplotype, thetaNT = cov(Yo, NTp + NTm) = 2*delta_o*gt + 2*a_o*ic + w_o: twice these.
    # The parents' own two haplotypes are uncorrelated (gc = hc = 0: first AM event), so per haplotype
    #   cov(Yo, NTp) = delta_o*cov(Tm, NTp) + a_o*cov(LTm, NTp) + cov(Fo, NTp) = delta_o*gt + a_o*ic + w_o/2
    # (ic here is the across-parent LGS-PGS covariance Omega_td*mu*Gamma_td)
    thetaNT <- mxAlgebra(delta_o * gt + a_o * ic + 0.5 * w_o, name="thetaNT")
    # cov(Yo, Tp) adds delta_o*var(Tp) = delta_o*k
    thetaT  <- mxAlgebra(delta_o * k + thetaNT, name="thetaT")

    # ---- Cascade AM: latent mating phenotype gamma~ ----
    # Tilde paths = multiplier x the path into Y (ETFD Cascade MVN scripts); see the header for AM_G / AM_E.
    AM_G <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_G_free, values=AM_G_value, label="AMGenMulti", name="AM_G")
    AM_E <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_E_free, values=AM_E_value, label="AMEnvMulti", name="AM_E")
    delta_td <- mxAlgebra(delta_p * AM_G, name="delta_td")
    a_td     <- mxAlgebra(a_p * AM_G, name="a_td")
    # Shortcuts to gamma~ (PDF Part IV Section 15): Omega_td = cov(gamma~, [N]T), Gamma_td = cov(gamma~, L[N]T),
    # zeta = cov(gamma~, F), tau = cov(Y, gamma~), Vgamma = var(gamma~)
    Omega_td <- mxAlgebra(delta_td * k + 0.5 * AM_E * w, name="Omega_td")
    Gamma_td <- mxAlgebra(a_td * j + 0.5 * AM_E * v, name="Gamma_td")
    zeta     <- mxAlgebra(delta_td * w + a_td * v + AM_E * VF_p_Algebra, name="zeta")
    tau      <- mxAlgebra(2 * a_p * Gamma_td + 2 * delta_p * Omega_td + zeta + AM_E * VE_p, name="tau")
    Vgamma   <- mxAlgebra(2 * a_td * Gamma_td + 2 * delta_td * Omega_td + AM_E * zeta + AM_E^2 * VE_p, name="Vgamma")

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
        w_o, v_o, VF_o,
        AM_G, AM_E, delta_td, a_td, Omega_td, Gamma_td, zeta, tau, Vgamma,
        VY, VE_p, VE_o, delta_p, a_p, delta_o, a_o, k, j, Omega_p, Gamma_p, mu, gc, hc, f,
        VY_p_Algebra, VF_p_Algebra, VY_o_Algebra, Omega_p_Algebra, Gamma_p_Algebra,
        gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
        gt, ht, ic, w, v,
        h2mat_p, h2mat_o, rdr_left_p, rdr_right_p, rdr_left_o, rdr_right_o,
        thetaNT, thetaT, CovMatrix, Means, ModelExpectations, mxFitFunctionML(), Constraints
    )

    Model1 <- mxModel("Cascade_UniSEM_NonEquilibrium_DiffTrait_Latent", Params, mxData(observed=Example_Data, type="raw"))
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
