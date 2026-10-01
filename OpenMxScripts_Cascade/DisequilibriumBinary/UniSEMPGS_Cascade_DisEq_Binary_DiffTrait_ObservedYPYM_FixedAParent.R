## SEM-PGS CASCADE version of OpenMxScripts/DisequilibriumBinary/UniSEMPGS_DisEq_Binary_DiffTrait_ObservedYPYM_FixedAParent.R: the same model, except that spouses
## assort on the latent mating phenotype gamma~ = delta~*(T+NT) + a~*(LT+LNT) + 1~*F + 1~*E instead of on Y
## (SEM_PGS_Cascade_model.pdf, Part IV). Following the ETFD Cascade MVN scripts, each path into
## gamma~ is the path into Y times a multiplier: AM_G on the genetic paths (delta~ = AM_G*delta, a~ = AM_G*a)
## and AM_E on the non-genetic paths (the PDF's 1~ on F and E). AM_G = AM_E = 1: primary phenotypic AM;
## AM_E = 0: genetic homogamy; AM_G = 0: social homogamy. gamma~ has no scale of its own, so at least one
## multiplier must be fixed (default AM_E = 1, the analogue of the MVN script's fixed AM_U); mu is the
## copath between spouses' gamma~. By default AM_G is estimated; for genetic homogamy fix AM_G = 1 and free AM_E instead.
## Offspring side, derived from Yo = delta*(Tp+Tm) + a*(LTp+LTm) + Fo + Eo (see ../PDF_Errata.md): the offspring's
## haplotypes covary across parents (gt, ht, ic), their F comes from AM-mated parents (w_o, v_o, VF_o carry the
## tau*mu terms), they have their own phenotypic variance VY_o, and thetaT/thetaNT are per single-haplotype
## column (the papers' theta is the sum over the father's and the mother's haplotype). Yo_Yp follows PDF
## Part IV Section 19 (it includes f*tau^2*mu = f*cov(Ym, Yp), which the original DisEq Yo_Yp omits).
##
## Binary/liability-threshold version of UniSEMPGS_DisEq_DiffTrait_ObservedYPYM_FixedAParent.R:
## parent and offspring express DIFFERENT binary traits, the parental phenotypes ARE observed, and
## the trait is NOT in equilibrium (gc = hc = ic = 0, one generation of assortative mating). Yp1, Ym1
## and Yo1 each get their own free threshold. RDR identifies BOTH generations' 'a' (parental and
## offspring). Each trait's liability variance (VY_p, VY_o) is independently fixed to 1.
## Optional covariates (covars) enter as definition variables on the PGS
## means and on the thresholds, as in r2_omx_partial().

fitUniSEMPGS_Cascade_DisEq_Binary_DiffTrait_ObservedYPYM_FixedAParent <- function(data_path, h2_RDR_parent, h2_RDR_offspring, covars = NULL, AM_G_value = .5, AM_G_free = TRUE, AM_E_value = 1, AM_E_free = FALSE, feaTol = 1e-6, optTol = 1e-8, jitterMean = .2, jitterVar = .05, extraTries = 30, exhaustive = F){
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
    obs_vars    <- c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1")   # extra columns (e.g. covariates) are allowed
    binary_vars <- c("Yp1", "Ym1", "Yo1")
    for (v in binary_vars) Example_Data[[v]] <- mxFactor(Example_Data[[v]], levels = c(0, 1))
    if (length(covars) > 0) {
        missing_cov <- setdiff(covars, colnames(Example_Data))
        if (length(missing_cov) > 0) stop("Covariate column(s) not found in data: ", paste(missing_cov, collapse = ", "))
        Example_Data <- Example_Data[complete.cases(Example_Data[, covars, with = FALSE]), ]
    }

    # 1. Phenotypic and Residual Variances (Independently estimated)
    # Each trait's liability variance is FIXED at 1 (not free), independently for parent/offspring.
    VY_p  <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1, label="VY_p11", name="VY_p")
    VY_o  <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1, label="VY_o11", name="VY_o")
    VE_p  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4,  label="VE_p11", name="VE_p", lbound = .001)
    VE_o  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4,  label="VE_o11", name="VE_o", lbound = .001)

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
    gt_Algebra <- mxAlgebra(Omega_td * mu * Omega_td, name="gt_Algebra")
    ht_Algebra <- mxAlgebra(Gamma_td * mu * Gamma_td, name="ht_Algebra")
    ic_Algebra <- mxAlgebra(Omega_td * mu * Gamma_td, name="ic_Algebra")

    gt <- mxAlgebra(gt_Algebra, name="gt")
    ht <- mxAlgebra(ht_Algebra, name="ht")
    ic <- mxAlgebra(ic_Algebra, name="ic")

    # ---- Offspring generation ----
    # The offspring's F comes from AM-mated parents, so its covariance with the parents' haplotypes (w_o, v_o)
    # and its variance (VF_o) carry the cross-mate tau*mu terms (PDF Part IV Section 18: w2, v2, VF2).
    # (w, v and the VF algebra above belong to the parents, whose own parents mated at random.)
    w_o  <- mxAlgebra(2 * f * Omega_p + 2 * f * tau * mu * Omega_td, name="w_o")
    v_o  <- mxAlgebra(2 * f * Gamma_p + 2 * f * tau * mu * Gamma_td, name="v_o")
    VF_o <- mxAlgebra(2 * f^2 * VY_p + 2 * f^2 * tau^2 * mu, name="VF_o")
    # thetaNT / thetaT: covariance of Yo with ONE parental haplotype, i.e. with one data column (NTp1, Tp1, ...).
    # NOTE: in the papers (Balbona et al. 2021, Lyu et al., SEM_PGS_Cascade_model.pdf) theta is the SUM over the
    # father's and the mother's haplotype, thetaNT = cov(Yo, NTp + NTm) = 2*delta_o*gt + 2*a_o*ic + w_o: twice these.
    # The parents' own two haplotypes are uncorrelated (gc = hc = 0: first AM event), so per haplotype
    #   cov(Yo, NTp) = delta_o*cov(Tm, NTp) + a_o*cov(LTm, NTp) + cov(Fo, NTp) = delta_o*gt + a_o*ic + w_o/2
    # (ic here is the across-parent LGS-PGS covariance Omega_td*mu*Gamma_td)
    thetaNT <- mxAlgebra(delta_o * gt + a_o * ic + 0.5 * w_o, name="thetaNT")
    # cov(Yo, Tp) adds delta_o*var(Tp) = delta_o*k
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
    # Shortcuts to gamma~ (PDF Part IV Section 15): Omega_td = cov(gamma~, [N]T), Gamma_td = cov(gamma~, L[N]T),
    # zeta = cov(gamma~, F), tau = cov(Y, gamma~), Vgamma = var(gamma~)
    Omega_td <- mxAlgebra(delta_td * k + 0.5 * AM_E * w, name="Omega_td")
    Gamma_td <- mxAlgebra(a_td * j + 0.5 * AM_E * v, name="Gamma_td")
    zeta     <- mxAlgebra(delta_td * w + a_td * v + AM_E * VF_p_Algebra, name="zeta")
    tau      <- mxAlgebra(2 * a_p * Gamma_td + 2 * delta_p * Omega_td + zeta + AM_E * VE_p, name="tau")
    Vgamma   <- mxAlgebra(2 * a_td * Gamma_td + 2 * delta_td * Omega_td + AM_E * zeta + AM_E^2 * VE_p, name="Vgamma")

    # 8. Expected Covariance Matrix (7x7)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_p_Algebra, Yp_Ym,      Yo_Yp,     Omega_p, Omega_p, Yp_PGSm, Yp_PGSm),
        cbind(Yp_Ym,        VY_p_Algebra, Yo_Yp,   Yp_PGSm, Yp_PGSm, Omega_p, Omega_p),
        cbind(Yo_Yp,        Yo_Yp,      VY_o_Algebra, thetaT, thetaNT, thetaT,  thetaNT),
        cbind(Omega_p,      Yp_PGSm,    thetaT,    k+gc,    gc,      gt,      gt),
        cbind(Omega_p,      Yp_PGSm,    thetaNT,   gc,      k+gc,    gt,      gt),
        cbind(Yp_PGSm,      Omega_p,    thetaT,    gt,      gt,      k+gc,    gc),
        cbind(Yp_PGSm,      Omega_p,    thetaNT,   gt,      gt,      gc,      k+gc)),
        dimnames = list(obs_vars, obs_vars), name="expCov")

    # Means: fixed at 0 for the binary phenotypes (identification requires fixing either the mean
    # or the threshold; we fix the mean and let the threshold float), free for the continuous PGS.
    # Each observed phenotype gets its own free threshold.
    # Covariates (if any) enter as definition variables, as in r2_omx_partial(): every PGS mean
    # becomes mean + sum(b * covar) and every threshold becomes t0 + sum(g * covar).
    Means0 <- mxMatrix(type = "Full", nrow = 1, ncol = 7, free = c(F, F, F, T, T, T, T), values = 0,
        label = c("meanYp1", "meanYm1", "meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"),
        dimnames = list(NULL, obs_vars), name = "Means0")
    Th0 <- mxMatrix(type = "Full", nrow = 1, ncol = 3, free = TRUE, values = 0,
        labels = c("thresh_Yp11", "thresh_Ym11", "thresh_Yo11"), name = "Th0")
    if (length(covars) > 0) {
        nCov   <- length(covars)
        isCont <- !(obs_vars %in% binary_vars)
        defCov <- mxMatrix(type = "Full", nrow = 1, ncol = nCov, free = FALSE, labels = paste0("data.", covars), name = "defCov")
        bCov   <- mxMatrix(type = "Full", nrow = nCov, ncol = length(obs_vars), free = matrix(isCont, nCov, length(obs_vars), byrow = TRUE),
            values = 0, labels = outer(covars, obs_vars, function(cv, vr) ifelse(vr %in% binary_vars, NA, paste0("b_", vr, "_", cv))), name = "bCov")
        gCov   <- mxMatrix(type = "Full", nrow = nCov, ncol = length(binary_vars), free = TRUE, values = 0,
            labels = outer(covars, binary_vars, function(cv, vr) paste0("g_", vr, "_", cv)), name = "gCov")
        Means   <- mxAlgebra(Means0 + defCov %*% bCov, dimnames = list(NULL, obs_vars), name = "expMeans")
        Th      <- mxAlgebra(Th0 + defCov %*% gCov, name = "Th")
        CovObjs <- list(Means0, Th0, defCov, bCov, gCov)
    } else {
        Means   <- mxAlgebra(Means0, dimnames = list(NULL, obs_vars), name = "expMeans")
        Th      <- mxAlgebra(Th0, name = "Th")
        CovObjs <- list(Means0, Th0)
    }

    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans", dimnames=obs_vars,
                                            thresholds="Th", threshnames=binary_vars)

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
        AM_G, AM_E, delta_td, a_td, Omega_td, Gamma_td, zeta, tau, Vgamma,
        VY_p, VY_o, VE_p, VE_o, delta_p, a_p, delta_o, a_o, k, j, Omega_p, Gamma_p, mu, gc, hc, f,
        VY_p_Algebra, VF_p_Algebra, VY_o_Algebra, Omega_p_Algebra, Gamma_p_Algebra,
        gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
        gt, ht, ic, w, v,
        h2mat_p, h2mat_o, rdr_left_p, rdr_right_p, rdr_left_o, rdr_right_o,
        thetaNT, thetaT, Yp_Ym, Yo_Yp, Yp_PGSm,
        CovMatrix, CovObjs, Means, Th, ModelExpectations, mxFitFunctionML(jointConditionOn = "continuous"), Constraints
    )

    Model1 <- mxModel("Cascade_UniSEM_Binary_NonEquilibrium_DiffTrait_Observed_FixedAParent", Params, mxData(observed=Example_Data, type="raw"))
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
