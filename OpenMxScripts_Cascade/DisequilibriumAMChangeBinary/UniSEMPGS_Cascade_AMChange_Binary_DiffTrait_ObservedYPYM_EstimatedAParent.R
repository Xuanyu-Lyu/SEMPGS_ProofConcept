## Univariate SEM-PGS Cascade model: different traits in parents and offspring, parental phenotypes observed.
## The parental trait drives assortative mating and vertical transmission; the offspring trait is the outcome.
## The parents' latent genetic effect a_p is estimated from the parental phenotypes.
## Change in assortment strength (two-regime disequilibrium, SEM_PGS_Cascade_model.pdf Part VII): the grandparents
## (_pp) and every generation before them mated with copath mu_pp long enough to reach equilibrium; the parents (_p)
## are an ordinary member of that population and mate once with their own copath mu_p; the offspring (_o) are the
## first generation formed under mu_p. mu_pp = 0 gives the Disequilibrium/ model, mu_pp = mu_p the Equilibrium/ one.
## Spouses assort on a latent mating phenotype gamma~ = AM_G*(delta*(T + NT) + a*(LT + LNT)) + AM_E*(F + E).
## Binary trait: fit as a liability-threshold model with both traits' liability variances fixed at 1.
##
## Estimates (returned by summary(fit, verbose = TRUE)):
##   delta_p, a_p, delta_o, a_o  PGS and latent-genetic-score effects on the parental / offspring trait
##   f                           vertical transmission from each parent's phenotype to the offspring's environment F
##   mu_pp                       copath between the grandparents' (and every earlier generation's) mating phenotypes
##   mu_p                        copath between the parents' mating phenotypes gamma~
##   AM_G, AM_E                  mating-mechanism multipliers, when freed (AMGenMulti, AMEnvMulti)
##   VE_p, VE_o                  residual variances of the parental and offspring trait
##   thresholds                  liability thresholds (thresh_Yp, thresh_Ym, thresh_Yo)
##   plus Omega_p, Gamma_p, gc_p, hc_p, ic_p, w_p, v_p, VF_p: the parents' covariance terms, held to their model
##   values by constraints.
##
## Input:
##   data_path   text file with a header row (read by data.table::fread) and the columns
##                 Yp1, Ym1    father's and mother's phenotype on the parental trait (0/1)
##                 Yo1         offspring phenotype on the offspring trait (0/1)
##                 Tp1, NTp1   father's transmitted and non-transmitted haplotypic PGS
##                 Tm1, NTm1   mother's transmitted and non-transmitted haplotypic PGS
##               Other columns (e.g. covariates) may be present.
##               Each haplotypic PGS is assumed to have variance k = .5 in the base population.
##   h2_RDR_offspring  RDR heritability of the offspring trait; identifies a_o
##   covars      optional covariates, each applied only to the person it belongs to: a named list from data variables
##               to that person's covariate columns, e.g. list(Yp1 = "age_p", Ym1 = "age_m", Yo1 = c("age_o", "sex_o"),
##               Tp1 = "PC1_p", NTp1 = "PC1_p", Tm1 = "PC1_m", NTm1 = "PC1_m"). A covariate shifts the threshold or PGS
##               mean of only the variables it is listed under; variables not listed get no covariates.
##   AM_G_value, AM_G_free, AM_E_value, AM_E_free  multipliers on the genetic (AM_G) and the F + E (AM_E) paths into
##               gamma~: (1, 1) primary phenotypic AM, (1, 0) genetic homogamy, (0, 1) social homogamy. At least one
##               must be fixed. Default: AM_E fixed at 1, AM_G estimated (for genetic homogamy fix AM_G = 1, free AM_E).
##               The mechanism is assumed to be the same in both regimes; only the strength (mu_pp vs mu_p) changes.
##   feaTol, optTol, extraTries, exhaustive, jitterMean, jitterVar  optimizer settings (NPSOL, mxTryHard)

fitUniSEMPGS_Cascade_AMChange_Binary_DiffTrait_ObservedYPYM_EstimatedAParent <- function(data_path, h2_RDR_offspring, covars = NULL, AM_G_value = .5, AM_G_free = TRUE, AM_E_value = 1, AM_E_free = FALSE, feaTol = 1e-6, optTol = 1e-8, jitterMean = .2, jitterVar = .05, extraTries = 30, exhaustive = F){
    library(OpenMx)
    library(data.table)
    library(stringr)

    mxOption(NULL,"Calculate Hessian","Yes")
    mxOption(NULL,"Standard Errors","Yes")
    mxOption(NULL,"Default optimizer","NPSOL")
    mxOption(NULL,"Feasibility tolerance",as.character(feaTol))
    mxOption(NULL,"Optimality tolerance",as.character(optTol))
    mxOption(NULL,"Number of Threads", value = parallel::detectCores())

    obs_vars     <- c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1")   # extra columns (e.g. covariates) are allowed
    binary_vars  <- c("Yp1", "Ym1", "Yo1")
    Example_Data <- fread(data_path, header = T)
    for (v in binary_vars) Example_Data[[v]] <- mxFactor(Example_Data[[v]], levels = c(0, 1))
    # covars: a named list from data variables to that person's own covariate columns (see the header); each covariate
    # adjusts only the variables it is listed under
    if (length(covars) > 0 && (!is.list(covars) || is.null(names(covars)) || any(names(covars) == "")))
        stop("covars must be a named list from data variables to that person's covariate columns, e.g. ",
             "list(Yo1 = c(\"age_o\", \"sex_o\"), Tp1 = \"PC1_p\", NTp1 = \"PC1_p\")")
    badVar <- setdiff(names(covars), obs_vars)
    if (length(badVar) > 0) stop("covars names must be data variables (", paste(obs_vars, collapse = ", "), "): ",
                                 paste(badVar, collapse = ", "))
    covCols <- unique(unlist(covars, use.names = FALSE))
    if (length(covCols) > 0) {
        missing_cov <- setdiff(covCols, colnames(Example_Data))
        if (length(missing_cov) > 0) stop("Covariate column(s) not found in data: ", paste(missing_cov, collapse = ", "))
        Example_Data <- Example_Data[complete.cases(Example_Data[, covCols, with = FALSE]), ]
    }

    # Starting values on the liability scale, where the parents' VY_p is fixed at 1 (the binary phenotypes give no
    # phenotypic moments to start from); gc_p starts at the sample cov(Tp, NTp).
    VY_s <- 1; VYo_s <- 1; Omega_s <- .15; delta_s <- .2
    gc_s <- mean(c(cov(Example_Data$Tp1, Example_Data$NTp1, use = "complete.obs"),
                   cov(Example_Data$Tm1, Example_Data$NTm1, use = "complete.obs")))

    # ---- Structural parameters (the same in every generation) ----
    delta_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=Omega_s, label="delta_p11", name="delta_p")
    a_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4*sqrt(VY_s), label="a_p11", name="a_p", lbound = .001)
    delta_o <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=delta_s, label="delta_o11", name="delta_o")
    a_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4*sqrt(VYo_s), label="a_o11", name="a_o", lbound = .001)
    f       <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="f11", name="f")
    VE_p    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=VY_s/3, label="VE_p11", name="VE_p", lbound = .001)
    VE_o    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=VYo_s/3, label="VE_o11", name="VE_o", lbound = .001)
    k       <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    j       <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")
    # k and j are fixed at .5, which assumes the PGS (and, by convention, the latent LGS) is standardised in the
    # base population, before any AM. Each parental haplotypic PGS then has variance k + gc_p (gc_p is inherited from
    # the grandparents' assortment), and the offspring PGS has variance 2(k + gc_p) + 2gt_p. If the haplotypic PGS is
    # instead scaled to variance 1/2 in the parents, replace the two fixed matrices above with
    # k <- mxAlgebra(.5 - gc_p, name = "k")
    # j <- mxAlgebra(.5 - hc_p, name = "j")

    # ---- Cascade AM: latent mating phenotype gamma~ (on the parental trait) ----
    # Tilde paths = multiplier x the path into Y (ETFD Cascade MVN scripts); see the header for AM_G / AM_E.
    AM_G <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_G_free, values=AM_G_value, label="AMGenMulti", name="AM_G")
    AM_E <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_E_free, values=AM_E_value, label="AMEnvMulti", name="AM_E")
    delta_td <- mxAlgebra(delta_p * AM_G, name="delta_td")
    a_td     <- mxAlgebra(a_p * AM_G, name="a_td")

    # ---- Grandparents (_pp) and parents (_p) ----
    # The parents are an ordinary member of the grandparents' mu_pp equilibrium (PDF Part VII Sections 34-35), so their
    # own quantities are the Equilibrium/ fixed point at mu_pp: free parameters held to their algebra by constraints.
    mu_pp   <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1/VY_s, label="mu_pp11", name="mu_pp")
    # Liability variance is FIXED at 1 (not free): identifies the scale of the binary/threshold model.
    VY_p    <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1, label="VY_p11", name="VY_p")
    VF_p    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.05*VY_s, label="VF_p11", name="VF_p")
    Omega_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=Omega_s, label="Omega_p11", name="Omega_p")
    Gamma_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=Omega_s, label="Gamma_p11", name="Gamma_p")
    w_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1*sqrt(VY_s), label="w_p11", name="w_p")
    v_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1*sqrt(VY_s), label="v_p11", name="v_p")
    gc_p    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=gc_s, label="gc_p11", name="gc_p")
    hc_p    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=gc_s, label="hc_p11", name="hc_p")
    ic_p    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=gc_s, label="ic_p11", name="ic_p")

    VY_p_Algebra    <- mxAlgebra(2 * delta_p * Omega_p + 2 * a_p * Gamma_p + delta_p * w_p + a_p * v_p + VF_p + VE_p, name="VY_p_Algebra")
    VF_p_Algebra    <- mxAlgebra(2 * f^2 * VY_p + 2 * f^2 * tau_p^2 * mu_pp, name="VF_p_Algebra")
    Omega_p_Algebra <- mxAlgebra(2 * delta_p * gc_p + 2 * a_p * ic_p + delta_p * k + 0.5 * w_p, name="Omega_p_Algebra")
    Gamma_p_Algebra <- mxAlgebra(2 * a_p * hc_p + 2 * delta_p * ic_p + a_p * j + 0.5 * v_p, name="Gamma_p_Algebra")
    w_p_Algebra     <- mxAlgebra(2 * f * Omega_p + 2 * f * tau_p * mu_pp * Omega_td_p, name="w_p_Algebra")
    v_p_Algebra     <- mxAlgebra(2 * f * Gamma_p + 2 * f * tau_p * mu_pp * Gamma_td_p, name="v_p_Algebra")
    # The grandparents' assortment made their two gametes covary (trans: gt_pp, ht_pp, it_pp); at equilibrium that is
    # each parent's within-person (cis) covariance between their own two haplotypes.
    gt_pp <- mxAlgebra(Omega_td_p * mu_pp * Omega_td_p, name="gt_pp")
    ht_pp <- mxAlgebra(Gamma_td_p * mu_pp * Gamma_td_p, name="ht_pp")
    it_pp <- mxAlgebra(Omega_td_p * mu_pp * Gamma_td_p, name="it_pp")

    # The parents' shortcuts to gamma~ (PDF Part III Section 8): Omega_td_p = cov(gamma~, [N]T),
    # Gamma_td_p = cov(gamma~, L[N]T), zeta_p = cov(gamma~, F), tau_p = cov(Y, gamma~), Vgamma_p = var(gamma~)
    Omega_td_p <- mxAlgebra(2 * delta_td * gc_p + 2 * a_td * ic_p + delta_td * k + 0.5 * AM_E * w_p, name="Omega_td_p")
    Gamma_td_p <- mxAlgebra(2 * a_td * hc_p + 2 * delta_td * ic_p + a_td * j + 0.5 * AM_E * v_p, name="Gamma_td_p")
    zeta_p     <- mxAlgebra(delta_td * w_p + a_td * v_p + AM_E * VF_p, name="zeta_p")
    tau_p      <- mxAlgebra(2 * a_p * Gamma_td_p + 2 * delta_p * Omega_td_p + zeta_p + AM_E * VE_p, name="tau_p")
    Vgamma_p   <- mxAlgebra(2 * a_td * Gamma_td_p + 2 * delta_td * Omega_td_p + AM_E * zeta_p + AM_E^2 * VE_p, name="Vgamma_p")

    # ---- The parents' own mating, with copath mu_p (PDF Part VII Sections 36-37) ----
    mu_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1/VY_s, label="mu_p11", name="mu_p")
    # trans covariances between the father's and the mother's haplotypes
    gt_p <- mxAlgebra(Omega_td_p * mu_p * Omega_td_p, name="gt_p")
    ht_p <- mxAlgebra(Gamma_td_p * mu_p * Gamma_td_p, name="ht_p")
    it_p <- mxAlgebra(Omega_td_p * mu_p * Gamma_td_p, name="it_p")
    Yp_Ym   <- mxAlgebra(tau_p * mu_p * tau_p, name="Yp_Ym")
    Yp_PGSm <- mxAlgebra(tau_p * mu_p * Omega_td_p, name="Yp_PGSm")      # = cov(Ym, [N]Tp)

    # ---- Offspring (_o) ----
    # The offspring's F comes from parents who mated with copath mu_p, so its covariance with the parents' haplotypes
    # (w_o, v_o) and its variance (VF_o) carry the cross-mate tau_p*mu_p terms.
    w_o  <- mxAlgebra(2 * f * Omega_p + 2 * f * tau_p * mu_p * Omega_td_p, name="w_o")
    v_o  <- mxAlgebra(2 * f * Gamma_p + 2 * f * tau_p * mu_p * Gamma_td_p, name="v_o")
    VF_o <- mxAlgebra(2 * f^2 * VY_p + 2 * f^2 * tau_p^2 * mu_p, name="VF_o")
    # thetaNT / thetaT: covariance of Yo with ONE parental haplotype, i.e. with one data column (NTp1, Tp1, ...).
    # NOTE: in the papers theta is the SUM over the father's and the mother's haplotype, cov(Yo, NTp + NTm): twice these.
    # A parent's gamete carries the parent's inherited cis covariance (cov(Tp, NTp) = gc_p) and covaries with the
    # spouse's by gt_p, so per haplotype
    #   cov(Yo, NTp) = delta_o*(gc_p + gt_p) + a_o*(ic_p + it_p) + w_o/2
    # (PDF Part VII Section 38 uses 2*delta*gt_p + 2*a*it_p instead, which holds only when mu_pp = mu_p.)
    thetaNT <- mxAlgebra(delta_o * (gc_p + gt_p) + a_o * (ic_p + it_p) + 0.5 * w_o, name="thetaNT")
    # cov(Yo, Tp) adds delta_o*k
    thetaT  <- mxAlgebra(delta_o * k + thetaNT, name="thetaT")
    # The offspring trait's variance: each offspring haplotype has variance k + gc_p, the two covary by gt_p (and ht_p,
    # it_p across the LGS), and the offspring's F has variance VF_o and covariances w_o, v_o with them
    # the offspring trait's liability variance is FIXED at 1 as well
    VY_o         <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1, label="VY_o11", name="VY_o")
    VY_o_Algebra <- mxAlgebra(2 * delta_o^2 * (k + gc_p + gt_p) + 2 * a_o^2 * (j + hc_p + ht_p) + 4 * delta_o * a_o * (ic_p + it_p) + 2 * delta_o * w_o + 2 * a_o * v_o + VF_o + VE_o, name="VY_o_Algebra")
    # Offspring-parent covariance (PDF Part VII Section 38, with the offspring trait's delta_o, a_o)
    Yo_Yp <- mxAlgebra(delta_o * Omega_p + a_o * Gamma_p + f * VY_p + tau_p * mu_p * (delta_o * Omega_td_p + a_o * Gamma_td_p + f * tau_p), name="Yo_Yp")

    # Identification of a_o via RDR
    rdr_left_o  <- mxAlgebra((2*a_o^2*j + 2*delta_o^2*k) * (2*a_o^2*j + 2*delta_o^2*k + VE_o), name="rdr_left_o")
    h2mat_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_offspring, name="h2mat_o")
    rdr_right_o <- mxAlgebra(h2mat_o * VY_o, name="rdr_right_o")

    Constraints <- list(
        mxConstraint(VY_p == VY_p_Algebra, name="VY_p_Constraint"),
        mxConstraint(VY_o == VY_o_Algebra, name="VY_o_Constraint"),
        mxConstraint(VF_p == VF_p_Algebra, name="VF_p_Constraint"),
        mxConstraint(Omega_p == Omega_p_Algebra, name="Omega_p_Constraint"),
        mxConstraint(Gamma_p == Gamma_p_Algebra, name="Gamma_p_Constraint"),
        mxConstraint(w_p == w_p_Algebra, name="w_p_Constraint"),
        mxConstraint(v_p == v_p_Algebra, name="v_p_Constraint"),
        mxConstraint(gc_p == gt_pp, name="gc_p_Constraint"),
        mxConstraint(hc_p == ht_pp, name="hc_p_Constraint"),
        mxConstraint(ic_p == it_pp, name="ic_p_Constraint"),
        mxConstraint(rdr_left_o == rdr_right_o, name="rdr_o_con")
    )

    # Expected Covariance Matrix (7x7)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_p,    Yp_Ym,   Yo_Yp,   Omega_p, Omega_p, Yp_PGSm, Yp_PGSm), # Yp
        cbind(Yp_Ym,   VY_p,    Yo_Yp,   Yp_PGSm, Yp_PGSm, Omega_p, Omega_p), # Ym
        cbind(Yo_Yp,   Yo_Yp,   VY_o,    thetaT,  thetaNT, thetaT,  thetaNT), # Yo
        cbind(Omega_p, Yp_PGSm, thetaT,  k+gc_p,  gc_p,    gt_p,    gt_p),    # Tp
        cbind(Omega_p, Yp_PGSm, thetaNT, gc_p,    k+gc_p,  gt_p,    gt_p),    # NTp
        cbind(Yp_PGSm, Omega_p, thetaT,  gt_p,    gt_p,    k+gc_p,  gc_p),    # Tm
        cbind(Yp_PGSm, Omega_p, thetaNT, gt_p,    gt_p,    gc_p,    k+gc_p)), # NTm
        dimnames=list(obs_vars, obs_vars), name="expCov")

    # Means: fixed at 0 for the binary phenotypes (identification requires fixing either the mean
    # or the threshold; we fix the mean and let the threshold float), free for the continuous PGS.
    # Each observed phenotype gets its own free threshold.
    # Covariates (if any) enter as definition variables, as in r2_omx_partial(): a PGS mean becomes mean + sum(b * covar)
    # and a threshold becomes t0 + sum(g * covar), over the covariates covars lists for that variable; every other
    # covariate-variable effect is fixed at 0.
    Means0 <- mxMatrix(type = "Full", nrow = 1, ncol = 7, free = c(F, F, F, T, T, T, T), values = 0,
        label = c("meanYp1", "meanYm1", "meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"),
        dimnames = list(NULL, obs_vars), name = "Means0")
    Th0 <- mxMatrix(type = "Full", nrow = 1, ncol = 3, free = TRUE, values = 0,
        labels = c("thresh_Yp11", "thresh_Ym11", "thresh_Yo11"), name = "Th0")
    if (length(covCols) > 0) {
        nCov   <- length(covCols)
        uses   <- matrix(sapply(obs_vars, function(vr) covCols %in% covars[[vr]]), nrow = nCov,
                         dimnames = list(covCols, obs_vars))                 # does covariate i adjust variable v?
        bFree  <- uses & matrix(!(obs_vars %in% binary_vars), nCov, length(obs_vars), byrow = TRUE)
        gFree  <- uses[, binary_vars, drop = FALSE]
        defCov <- mxMatrix(type = "Full", nrow = 1, ncol = nCov, free = FALSE, labels = paste0("data.", covCols), name = "defCov")
        bCov   <- mxMatrix(type = "Full", nrow = nCov, ncol = length(obs_vars), free = bFree, values = 0,
            labels = ifelse(bFree, outer(covCols, obs_vars, function(cv, vr) paste0("b_", vr, "_", cv)), NA), name = "bCov")
        gCov   <- mxMatrix(type = "Full", nrow = nCov, ncol = length(binary_vars), free = gFree, values = 0,
            labels = ifelse(gFree, outer(covCols, binary_vars, function(cv, vr) paste0("g_", vr, "_", cv)), NA), name = "gCov")
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

    Params <- list(
        delta_p, a_p, delta_o, a_o, f, VE_p, VE_o, k, j, AM_G, AM_E, delta_td, a_td,
        mu_pp, VY_p, VF_p, Omega_p, Gamma_p, w_p, v_p, gc_p, hc_p, ic_p,
        VY_p_Algebra, VF_p_Algebra, Omega_p_Algebra, Gamma_p_Algebra, w_p_Algebra, v_p_Algebra,
        gt_pp, ht_pp, it_pp, Omega_td_p, Gamma_td_p, zeta_p, tau_p, Vgamma_p,
        mu_p, gt_p, ht_p, it_p, Yp_Ym, Yp_PGSm,
        w_o, v_o, VF_o, thetaNT, thetaT, VY_o, VY_o_Algebra, Yo_Yp,
        h2mat_o, rdr_left_o, rdr_right_o,
        Constraints, CovMatrix, CovObjs, Means, Th, ModelExpectations, mxFitFunctionML(jointConditionOn = "continuous")
    )

    Model1 <- mxModel("Cascade_UniSEM_Binary_AMChange_DiffTrait_Observed_RDR", Params, mxData(observed=Example_Data, type="raw"))
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
