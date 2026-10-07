## Univariate SEM-PGS Cascade model: the same trait in parents and offspring, parental phenotypes observed.
## Disequilibrium: vertical transmission is at equilibrium; the parents are the first assortatively mated generation.
## Spouses assort on a latent mating phenotype gamma~ = AM_G*(delta*(T + NT) + a*(LT + LNT)) + AM_E*(F + E).
## Binary trait: fit as a liability-threshold model with the parents' liability variance fixed at 1.
##
## Estimates (returned by summary(fit, verbose = TRUE)):
##   delta, a    effects of the haplotypic PGS and of the latent genetic score on the phenotype
##   f           vertical transmission from each parent's phenotype to the offspring's environment F
##   mu          copath between spouses' mating phenotypes gamma~
##   AM_G, AM_E  mating-mechanism multipliers, when freed (AMGenMulti, AMEnvMulti)
##   VE          residual variance
##   thresholds  liability thresholds (thresh_Yp, thresh_Ym, thresh_Yo)
##   plus Omega, Gamma: covariance terms held to their model values by constraints.
##
## Input:
##   data_path   text file with a header row (read by data.table::fread) and the columns
##                 Yp1, Ym1    father's and mother's phenotype (0/1)
##                 Yo1         offspring phenotype (0/1)
##                 Tp1, NTp1   father's transmitted and non-transmitted haplotypic PGS
##                 Tm1, NTm1   mother's transmitted and non-transmitted haplotypic PGS
##               Other columns (e.g. covariates) may be present.
##               Each haplotypic PGS is assumed to have variance k = .5 in the base population.
##   covars      optional covariates, each applied only to the person it belongs to: a named list from data variables
##               to that person's covariate columns, e.g. list(Yp1 = "age_p", Ym1 = "age_m", Yo1 = c("age_o", "sex_o"),
##               Tp1 = "PC1_p", NTp1 = "PC1_p", Tm1 = "PC1_m", NTm1 = "PC1_m"). A covariate shifts the threshold or PGS
##               mean of only the variables it is listed under; variables not listed get no covariates.
##   AM_G_value, AM_G_free, AM_E_value, AM_E_free  multipliers on the genetic (AM_G) and the F + E (AM_E) paths into
##               gamma~: (1, 1) primary phenotypic AM, (1, 0) genetic homogamy, (0, 1) social homogamy. At least one
##               must be fixed. Default: AM_E fixed at 1, AM_G estimated (for genetic homogamy fix AM_G = 1, free AM_E).
##   feaTol, optTol, extraTries, exhaustive, jitterMean, jitterVar  optimizer settings (NPSOL, mxTryHard)

fitUniSEMPGS_Cascade_DisEq_Binary_SameTrait_ObservedYPYM <- function(data_path, covars = NULL, AM_G_value = .5, AM_G_free = TRUE, AM_E_value = 1, AM_E_free = FALSE, feaTol = 1e-6, optTol = 1e-8, jitterMean = .2, jitterVar = .05, extraTries = 30, exhaustive = F){
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

    # --- Variances ---
    # Liability variance is FIXED at 1 (not free): identifies the scale of the binary/threshold model.
    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1, label="VY11", name="VY")
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="VE11", name="VE", lbound = .001)

    # Equation for VY from screenshot: 2*Omega*delta + 2*Gamma*a + delta*w + a*v + VF + VE
    VY_Algebra <- mxAlgebra(2 * Omega * delta + 2 * Gamma * a + delta * w + a * v + VF_Algebra + VE, name="VY_Algebra")
    # Equation for VF: f*VY*f (assuming father/mother equality)
    VF_Algebra <- mxAlgebra(2 * f^2 * VY, name="VF_Algebra")

    VY_Constraint <- mxConstraint(VY == VY_Algebra, name='VY_Constraint')

    # --- Genetic effects & Path Coefficients ---
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="delta11", name="delta")
    a     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="a11", name="a", lbound = .001)
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    j     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")

    # Omega and Gamma Equations from screenshot
    # Omega = delta*k + 0.5*w
    # Gamma = a*j + 0.5*v
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="Omega11", name="Omega")
    Gamma <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="Gamma11", name="Gamma")

    Omega_Algebra <- mxAlgebra(delta * k + 0.5 * w , name="Omega_Algebra")
    Gamma_Algebra <- mxAlgebra(a * j + 0.5 * v, name="Gamma_Algebra")

    Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
    Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')

    # --- Assortative Mating (Non-Equilibrium) ---
    mu    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="mu11", name="mu")

    # From screenshot: gc = hc = ic = 0 (No AM in previous generation)
    gc    <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=0, label="gc11", name="gc")
    hc    <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=0, label="hc11", name="hc")

    # gt, ht, ic are functions of mu and current Omega/Gamma. (itlo and itol -- the
    # transmitted-vs-non-transmitted haplotype covariances -- collapse to a single quantity 'ic'
    # here because Gamma*mu*Omega == Omega*mu*Gamma for scalars; the itlo/itol split only matters
    # once Omega and Gamma become asymmetric matrices in the bivariate model.)
    gt_Algebra <- mxAlgebra(Omega_td * mu * Omega_td, name="gt_Algebra")
    ht_Algebra <- mxAlgebra(Gamma_td * mu * Gamma_td, name="ht_Algebra")
    ic_Algebra <- mxAlgebra(Omega_td * mu * Gamma_td, name="ic_Algebra")

    gt <- mxAlgebra(gt_Algebra, name="gt")
    ht <- mxAlgebra(ht_Algebra, name="ht")
    ic <- mxAlgebra(ic_Algebra, name="ic")

    # --- Vertical Transmission ---
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="f11", name="f")

    # Equation: w = f*Omega + f*Omega (since parent traits/effects are equal)
    w_Algebra     <- mxAlgebra(2 * f * Omega, name="w_Algebra")
    v_Algebra     <- mxAlgebra(2 * f * Gamma, name="v_Algebra")

    w <- mxAlgebra(w_Algebra, name="w")
    v <- mxAlgebra(v_Algebra, name="v")

    # ---- Offspring generation ----
    # The offspring's F comes from AM-mated parents, so its covariance with the parents' haplotypes (w_o, v_o)
    # and its variance (VF_o) carry the cross-mate tau*mu terms (PDF Part IV Section 18: w2, v2, VF2).
    # (w, v and the VF algebra above belong to the parents, whose own parents mated at random.)
    w_o  <- mxAlgebra(2 * f * Omega + 2 * f * tau * mu * Omega_td, name="w_o")
    v_o  <- mxAlgebra(2 * f * Gamma + 2 * f * tau * mu * Gamma_td, name="v_o")
    VF_o <- mxAlgebra(2 * f^2 * VY + 2 * f^2 * tau^2 * mu, name="VF_o")
    # thetaNT / thetaT: covariance of Yo with ONE parental haplotype, i.e. with one data column (NTp1, Tp1, ...).
    # NOTE: in the papers (Balbona et al. 2021, Lyu et al., SEM_PGS_Cascade_model.pdf) theta is the SUM over the
    # father's and the mother's haplotype, thetaNT = cov(Yo, NTp + NTm) = 2*delta*gt + 2*a*ic + w_o: twice these.
    # The parents' own two haplotypes are uncorrelated (gc = hc = 0: first AM event), so per haplotype
    #   cov(Yo, NTp) = delta*cov(Tm, NTp) + a*cov(LTm, NTp) + cov(Fo, NTp) = delta*gt + a*ic + w_o/2
    # (ic here is the across-parent LGS-PGS covariance Omega_td*mu*Gamma_td)
    thetaNT <- mxAlgebra(delta * gt + a * ic + 0.5 * w_o, name="thetaNT")
    # cov(Yo, Tp) adds delta*var(Tp) = delta*k
    thetaT  <- mxAlgebra(delta * k + thetaNT, name="thetaT")
    # The offspring's own phenotypic variance: their two haplotypes covary by gt / ht / ic across parents but carry
    # no within-haplotype AM covariance yet, and their F has variance VF_o and covariances w_o, v_o with them
    VY_o_Algebra <- mxAlgebra(2 * delta^2 * k + 2 * delta^2 * gt + 2 * a^2 * j + 2 * a^2 * ht + 4 * delta * a * ic + 2 * delta * w_o + 2 * a * v_o + VF_o + VE, name="VY_o_Algebra")

    # Cross-parental expectations
    Yp_PGSm <- mxAlgebra(tau * mu * Omega_td, name="Yp_PGSm")
    Ym_PGSp <- mxAlgebra(tau * mu * Omega_td, name="Ym_PGSp")
    Yp_Ym   <- mxAlgebra(tau * mu * tau,    name="Yp_Ym")

    # Offspring-parent covariance (PDF Part IV Section 19): cov(Yo, Yp) = delta*[cov(Tp,Yp) + cov(Tm,Yp)]
    # + a*[cov(LTp,Yp) + cov(LTm,Yp)] + f*[VY + cov(Ym,Yp)], with the cross-mate terms tau*mu*(...)
    Yo_Yp   <- mxAlgebra(delta * Omega + a * Gamma + f * VY + tau * mu * (delta * Omega_td + a * Gamma_td + f * tau), name = "Yo_Yp")

    # ---- Cascade AM: latent mating phenotype gamma~ ----
    # Tilde paths = multiplier x the path into Y (ETFD Cascade MVN scripts); see the header for AM_G / AM_E.
    AM_G <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_G_free, values=AM_G_value, label="AMGenMulti", name="AM_G")
    AM_E <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_E_free, values=AM_E_value, label="AMEnvMulti", name="AM_E")
    delta_td <- mxAlgebra(delta * AM_G, name="delta_td")
    a_td     <- mxAlgebra(a * AM_G, name="a_td")
    # Shortcuts to gamma~ (PDF Part IV Section 15): Omega_td = cov(gamma~, [N]T), Gamma_td = cov(gamma~, L[N]T),
    # zeta = cov(gamma~, F), tau = cov(Y, gamma~), Vgamma = var(gamma~)
    Omega_td <- mxAlgebra(delta_td * k + 0.5 * AM_E * w, name="Omega_td")
    Gamma_td <- mxAlgebra(a_td * j + 0.5 * AM_E * v, name="Gamma_td")
    zeta     <- mxAlgebra(delta_td * w + a_td * v + AM_E * VF_Algebra, name="zeta")
    tau      <- mxAlgebra(2 * a * Gamma_td + 2 * delta * Omega_td + zeta + AM_E * VE, name="tau")
    Vgamma   <- mxAlgebra(2 * a_td * Gamma_td + 2 * delta_td * Omega_td + AM_E * zeta + AM_E^2 * VE, name="Vgamma")

    # --- Expected Covariance Matrix ---
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY,      Yp_Ym,   Yo_Yp,   Omega,   Omega,   Yp_PGSm, Yp_PGSm),
        cbind(Yp_Ym,   VY,      Yo_Yp,   Ym_PGSp, Ym_PGSp, Omega,   Omega),
        cbind(Yo_Yp,   Yo_Yp,   VY_o_Algebra, thetaT,  thetaNT, thetaT,  thetaNT),
        cbind(Omega,   Ym_PGSp, thetaT,  k+gc,    gc,      gt,      gt),
        cbind(Omega,   Ym_PGSp, thetaNT, gc,      k+gc,    gt,      gt),
        cbind(Yp_PGSm, Omega,   thetaT,  gt,      gt,      k+gc,    gc),
        cbind(Yp_PGSm, Omega,   thetaNT, gt,      gt,      gc,      k+gc)),
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
    Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )
    FitFunctionML   <- mxFitFunctionML(jointConditionOn = "continuous")

    Params <- list(
                w_o, v_o, VF_o, VY_o_Algebra,
                AM_G, AM_E, delta_td, a_td, Omega_td, Gamma_td, zeta, tau, Vgamma,
                VY, VE, delta, a, k, j, Omega, Gamma, mu, gc, hc, f,
                VY_Algebra, VF_Algebra, Omega_Algebra, Gamma_Algebra,
                gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
                gt, ht, ic, w, v,
                VY_Constraint, Omega_Constraint, Gamma_Constraint,
                thetaNT, thetaT, Yp_PGSm, Ym_PGSp, Yp_Ym, Yo_Yp,
                CovMatrix, CovObjs, Means, Th, ModelExpectations, FitFunctionML)

    Model1 <- mxModel("Cascade_UniSEM_Binary_NonEquilibrium", Params, Example_Data_Mx)
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1), intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
