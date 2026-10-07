## Univariate SEM-PGS model: different traits in parents and offspring, parental phenotypes observed.
## The parental trait drives assortative mating and vertical transmission; the offspring trait is the outcome.
## The parents' latent genetic effect a_p is estimated from the parental phenotypes.
## Disequilibrium: vertical transmission is at equilibrium; the parents are the first assortatively mated generation.
## Binary trait: fit as a liability-threshold model with both traits' liability variances fixed at 1.
##
## Estimates (returned by summary(fit, verbose = TRUE)):
##   delta_p, a_p, delta_o, a_o  PGS and latent-genetic-score effects on the parental / offspring trait
##   f                           vertical transmission from each parent's phenotype to the offspring's environment F
##   mu                          assortative-mating copath, cov(Yp, Ym) = mu*VY^2
##   VE_p, VE_o                  residual variances of the parental and offspring trait
##   thresholds                  liability thresholds (thresh_Yp, thresh_Ym, thresh_Yo)
##   plus Omega_p, Gamma_p: covariance terms held to their model values by constraints.
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
##   feaTol, optTol, extraTries, exhaustive, jitterMean, jitterVar  optimizer settings (NPSOL, mxTryHard)

fitUniSEMPGS_DisEq_Binary_DiffTrait_ObservedYPYM_EstimatedAParent <- function(data_path, h2_RDR_offspring, covars = NULL, feaTol = 1e-6, optTol = 1e-8, jitterMean = .2, jitterVar = .05, extraTries = 30, exhaustive = F){
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

    # 1. Phenotypic and Residual Variances (Separated)
    # Each trait's liability variance is FIXED at 1 (not free), independently for parent/offspring.
    VY_p  <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1,   label="VY_parent", name="VY_p")
    VY_o  <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=1,   label="VY_offspring", name="VY_o")
    VE_p  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4,  label="VE_parent", name="VE_p", lbound = .001)
    VE_o  <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4,  label="VE_offspring", name="VE_o", lbound = .001)

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

    # Identification of a_o via RDR
    rdr_left_o  <- mxAlgebra((2*a_o^2*j + 2*delta_o^2*k) * (2*a_o^2*j + 2*delta_o^2*k + VE_o), name="rdr_left_o")
    h2mat_o     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR_offspring, name="h2mat_o")
    rdr_right_o <- mxAlgebra(h2mat_o * VY_o, name="rdr_right_o")

    # 7. Assortment and Transmission Algebra (functions of mu and current Omega_p/Gamma_p). (itlo
    # and itol -- the transmitted-vs-non-transmitted haplotype covariances -- collapse to a single
    # quantity 'ic' here because Gamma_p*mu*Omega_p == Omega_p*mu*Gamma_p for scalars; the
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

    # Intergenerational Phenotypic Covariances (Between Parent Trait and Child Trait)
    Yp_Ym   <- mxAlgebra(VY_p * mu * VY_p, name="Yp_Ym")
    Yo_Yp   <- mxAlgebra((delta_o * Omega_p + a_o * Gamma_p + f * VY_p) * (1 + mu * VY_p), name = "Yo_Yp")

    # Cross-person PGS-Phenotype Covariances
    Yp_PGSm <- mxAlgebra(VY_p * mu * Omega_p, name="Yp_PGSm")

    # 8. Expected Covariance Matrix (7x7)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_p_Algebra, Yp_Ym,      Yo_Yp,     Omega_p, Omega_p, Yp_PGSm, Yp_PGSm), # Yp
        cbind(Yp_Ym,        VY_p_Algebra, Yo_Yp,   Yp_PGSm, Yp_PGSm, Omega_p, Omega_p), # Ym
        cbind(Yo_Yp,        Yo_Yp,      VY_o_Algebra, thetaT, thetaNT, thetaT,  thetaNT), # Yo
        cbind(Omega_p,      Yp_PGSm,    thetaT,    k+gc,    gc,      gt,      gt), # Tp
        cbind(Omega_p,      Yp_PGSm,    thetaNT,   gc,      k+gc,    gt,      gt), # NTp
        cbind(Yp_PGSm,      Omega_p,    thetaT,    gt,      gt,      k+gc,    gc), # Tm
        cbind(Yp_PGSm,      Omega_p,    thetaNT,   gt,      gt,      gc,      k+gc)), # NTm
        dimnames=list(obs_vars, obs_vars), name="expCov")

    # 9. Means and Expectations
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

    # 10. Constraints
    Constraints <- list(
        mxConstraint(VY_p == VY_p_Algebra, name="VY_p_eq"),
        mxConstraint(VY_o == VY_o_Algebra, name="VY_o_eq"),
        mxConstraint(Omega_p == Omega_p_Algebra, name="Om_p_eq"),
        mxConstraint(Gamma_p == Gamma_p_Algebra, name="Ga_p_eq"),
        mxConstraint(rdr_left_o == rdr_right_o, name="rdr_o_con")
    )

    Params <- list(
        w_o, v_o, VF_o,
        VY_p, VY_o, VE_p, VE_o, delta_p, a_p, delta_o, a_o, k, j, Omega_p, Gamma_p, mu, gc, hc, f,
        VY_p_Algebra, VF_p_Algebra, VY_o_Algebra, Omega_p_Algebra, Gamma_p_Algebra,
        gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
        gt, ht, ic, w, v,
        h2mat_o, rdr_left_o, rdr_right_o, thetaNT, thetaT, Yp_Ym, Yo_Yp, Yp_PGSm,
        CovMatrix, CovObjs, Means, Th, ModelExpectations, mxFitFunctionML(jointConditionOn = "continuous"), Constraints
    )

    Model1 <- mxModel("UniSEM_Binary_NonEquilibrium_DiffTrait_Observed_RDR", Params, mxData(observed=Example_Data, type="raw"))
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
