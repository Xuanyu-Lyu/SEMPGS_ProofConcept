## Model 1 of Balbona, Kim & Keller (2021), as a univariate SEM-PGS Cascade model.
## The PGS is assumed to capture all of the trait's heritability (no latent genetic score); only the offspring
## phenotype and the parents' PGS are observed.
## Change in assortment strength (two-regime disequilibrium, SEM_PGS_Cascade_model.pdf Part VII): the grandparents
## (_pp) and every generation before them mated with copath mu_pp long enough to reach equilibrium; the parents (_p)
## are an ordinary member of that population and mate once with their own copath mu_p; the offspring (_o) are the
## first generation formed under mu_p. mu_pp = 0 gives the Disequilibrium/ model, mu_pp = mu_p the Equilibrium/ one.
## Spouses assort on a latent mating phenotype gamma~ = AM_G*delta*(T + NT) + AM_E*(F + E).
##
## Estimates (returned by summary(fit, verbose = TRUE)):
##   delta       effect of the haplotypic PGS on the phenotype
##   f           vertical transmission from each parent's phenotype to the offspring's environment F
##   mu_pp       copath between the grandparents' (and every earlier generation's) mating phenotypes gamma~,
##               identified by the parents' cis covariance gc_p = cov(Tp, NTp)
##   mu_p        copath between the parents' mating phenotypes gamma~, identified by gt_p = cov(Tp, Tm)
##   VY_p, VE    the parents' phenotypic variance and the residual variance
##   plus Omega_p, gc_p, w_p, VF_p: the parents' covariance terms, held to their model values by constraints.
## The five moments VY_o, thetaT, thetaNT, gc_p and gt_p just-identify delta, f, mu_pp, mu_p and VE.
##
## Input:
##   data_path   text file with a header row (read by data.table::fread) and the columns
##                 Yo1         offspring phenotype
##                 Tp1, NTp1   father's transmitted and non-transmitted haplotypic PGS
##                 Tm1, NTm1   mother's transmitted and non-transmitted haplotypic PGS
##               Other columns are ignored.
##               Each haplotypic PGS is assumed to have variance k = .5 in the base population.
##   AM_G_value, AM_G_free, AM_E_value, AM_E_free  multipliers on the genetic (AM_G) and the F + E (AM_E) paths into
##               gamma~: (1, 1) primary phenotypic AM, (1, 0) genetic homogamy, (0, 1) social homogamy. At least one
##               must be fixed. Default: both fixed at 1. With latent parents the mechanism is not identified, so
##               fix both at the mechanism to be tested. It is assumed to be the same in both regimes.
##   feaTol, optTol, extraTries, exhaustive, jitterMean, jitterVar  optimizer settings (NPSOL, mxTryHard)

fitUniSEMPGS_Cascade_AMChange_Model1 <- function(data_path, AM_G_value = 1, AM_G_free = FALSE, AM_E_value = 1, AM_E_free = FALSE, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
    library(OpenMx)
    library(data.table)
    library(stringr)

    mxOption(NULL,"Calculate Hessian","Yes")
    mxOption(NULL,"Standard Errors","Yes")
    mxOption(NULL,"Default optimizer","NPSOL")
    mxOption(NULL,"Feasibility tolerance",as.character(feaTol))
    mxOption(NULL,"Optimality tolerance",as.character(optTol))
    mxOption(NULL,"Number of Threads", value = parallel::detectCores())

    obs_vars     <- c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")
    Example_Data <- fread(data_path, header = T)[, obs_vars, with = FALSE]

    # Starting values from the sample moments, which put the first attempt near the solution: thetaT - thetaNT =
    # delta*k, Omega_p is close to thetaT, the parents' VY_p is close to the offspring's variance, gc_p = cov(Tp, NTp).
    S <- cov(Example_Data, use = "pairwise.complete.obs")
    thetaT_s  <- mean(S["Yo1", c("Tp1", "Tm1")])
    thetaNT_s <- mean(S["Yo1", c("NTp1", "NTm1")])
    delta_s   <- (thetaT_s - thetaNT_s) / .5               # k = .5
    VY_s      <- S["Yo1", "Yo1"]
    gc_s      <- mean(c(S["Tp1", "NTp1"], S["Tm1", "NTm1"]))

    # ---- Structural parameters (the same in every generation) ----
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=delta_s, label="delta11", name="delta")
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="f11", name="f")
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=VY_s/2, label="VE11", name="VE", lbound = .001)
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    # k is fixed at .5, which assumes the PGS is standardised in the base population, before any AM. Each parental
    # haplotypic PGS then has variance k + gc_p (gc_p is inherited from the grandparents' assortment), and the offspring
    # PGS has variance 2(k + gc_p) + 2gt_p. If the haplotypic PGS is instead scaled to variance 1/2 in the parents,
    # replace the fixed matrix above with
    # k <- mxAlgebra(.5 - gc_p, name = "k")

    # ---- Cascade AM: latent mating phenotype gamma~ ----
    # Tilde paths = multiplier x the path into Y (ETFD Cascade MVN scripts); see the header for AM_G / AM_E.
    AM_G <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_G_free, values=AM_G_value, label="AMGenMulti", name="AM_G")
    AM_E <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_E_free, values=AM_E_value, label="AMEnvMulti", name="AM_E")
    delta_td <- mxAlgebra(delta * AM_G, name="delta_td")

    # ---- Grandparents (_pp) and parents (_p) ----
    # The parents are an ordinary member of the grandparents' mu_pp equilibrium (PDF Part VII Sections 34-35), so their
    # own quantities are the Equilibrium/ fixed point at mu_pp: free parameters held to their algebra by constraints.
    mu_pp   <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1/VY_s, label="mu_pp11", name="mu_pp")
    VY_p    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=VY_s, label="VY_p11", name="VY_p", lbound = .001)
    VF_p    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.05*VY_s, label="VF_p11", name="VF_p")
    Omega_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=thetaT_s, label="Omega_p11", name="Omega_p")
    w_p     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1*sqrt(VY_s), label="w_p11", name="w_p")
    gc_p    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=gc_s, label="gc_p11", name="gc_p")

    VY_p_Algebra    <- mxAlgebra(2 * delta * Omega_p + delta * w_p + VF_p + VE, name="VY_p_Algebra")
    VF_p_Algebra    <- mxAlgebra(2 * f^2 * VY_p + 2 * f^2 * tau_p^2 * mu_pp, name="VF_p_Algebra")
    Omega_p_Algebra <- mxAlgebra(2 * delta * gc_p + delta * k + 0.5 * w_p, name="Omega_p_Algebra")
    w_p_Algebra     <- mxAlgebra(2 * f * Omega_p + 2 * f * tau_p * mu_pp * Omega_td_p, name="w_p_Algebra")
    # The grandparents' assortment made their two gametes covary (trans: gt_pp); at equilibrium that is each parent's
    # within-person (cis) covariance between their own two haplotypes.
    gt_pp <- mxAlgebra(Omega_td_p * mu_pp * Omega_td_p, name="gt_pp")

    # The parents' shortcuts to gamma~ (PDF Part III Section 8, a = 0): Omega_td_p = cov(gamma~, [N]T),
    # zeta_p = cov(gamma~, F), tau_p = cov(Y, gamma~), Vgamma_p = var(gamma~)
    Omega_td_p <- mxAlgebra(2 * delta_td * gc_p + delta_td * k + 0.5 * AM_E * w_p, name="Omega_td_p")
    zeta_p     <- mxAlgebra(delta_td * w_p + AM_E * VF_p, name="zeta_p")
    tau_p      <- mxAlgebra(2 * delta * Omega_td_p + zeta_p + AM_E * VE, name="tau_p")
    Vgamma_p   <- mxAlgebra(2 * delta_td * Omega_td_p + AM_E * zeta_p + AM_E^2 * VE, name="Vgamma_p")

    # ---- The parents' own mating, with copath mu_p (PDF Part VII Section 36) ----
    mu_p <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1/VY_s, label="mu_p11", name="mu_p")
    gt_p <- mxAlgebra(Omega_td_p * mu_p * Omega_td_p, name="gt_p")   # trans covariance of the parents' haplotypes

    # ---- Offspring (_o) ----
    # The offspring's F comes from parents who mated with copath mu_p, so its covariance with the parents' haplotypes
    # (w_o) and its variance (VF_o) carry the cross-mate tau_p*mu_p terms.
    w_o  <- mxAlgebra(2 * f * Omega_p + 2 * f * tau_p * mu_p * Omega_td_p, name="w_o")
    VF_o <- mxAlgebra(2 * f^2 * VY_p + 2 * f^2 * tau_p^2 * mu_p, name="VF_o")
    # thetaNT / thetaT: covariance of Yo with ONE parental haplotype, i.e. with one data column (NTp1, Tp1, ...).
    # NOTE: in the papers theta is the SUM over the father's and the mother's haplotype, cov(Yo, NTp + NTm): twice these.
    # A parent's gamete carries the parent's inherited cis covariance (cov(Tp, NTp) = gc_p) and covaries with the
    # spouse's by gt_p, so per haplotype
    #   cov(Yo, NTp) = delta*[cov(Tp, NTp) + cov(Tm, NTp)] + cov(Fo, NTp) = delta*(gc_p + gt_p) + w_o/2
    # (PDF Part VII Section 38 uses 2*delta*gt_p instead, which holds only when mu_pp = mu_p.)
    thetaNT <- mxAlgebra(delta * (gc_p + gt_p) + 0.5 * w_o, name="thetaNT")
    # cov(Yo, Tp) adds delta*k
    thetaT  <- mxAlgebra(delta * k + thetaNT, name="thetaT")
    # The offspring's own phenotypic variance: each of their haplotypes has variance k + gc_p, the two covary by gt_p,
    # and their F has variance VF_o and covariance w_o with them
    VY_o_Algebra <- mxAlgebra(2 * delta^2 * (k + gc_p + gt_p) + 2 * delta * w_o + VF_o + VE, name="VY_o_Algebra")

    Constraints <- list(
        mxConstraint(VY_p == VY_p_Algebra, name="VY_p_Constraint"),
        mxConstraint(VF_p == VF_p_Algebra, name="VF_p_Constraint"),
        mxConstraint(Omega_p == Omega_p_Algebra, name="Omega_p_Constraint"),
        mxConstraint(w_p == w_p_Algebra, name="w_p_Constraint"),
        mxConstraint(gc_p == gt_pp, name="gc_p_Constraint")
    )

    # Expected covariances matrix (5x5: Yo1, Tp1, NTp1, Tm1, NTm1)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_o_Algebra, thetaT,  thetaNT, thetaT,  thetaNT),  # Yo1
        cbind(thetaT,       k+gc_p,  gc_p,    gt_p,    gt_p),     # Tp1
        cbind(thetaNT,      gc_p,    k+gc_p,  gt_p,    gt_p),     # NTp1
        cbind(thetaT,       gt_p,    gt_p,    k+gc_p,  gc_p),     # Tm1
        cbind(thetaNT,      gt_p,    gt_p,    gc_p,    k+gc_p)),  # NTm1
        dimnames = list(obs_vars, obs_vars), name="expCov")

    Means <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = TRUE, values = 0,
        label = c("meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"),
        dimnames = list(NULL, obs_vars), name = "expMeans")

    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans", dimnames=obs_vars)
    Example_Data_Mx   <- mxData(observed=Example_Data, type="raw")
    FitFunctionML     <- mxFitFunctionML()

    Params <- list(
                delta, f, VE, k, AM_G, AM_E, delta_td,
                mu_pp, VY_p, VF_p, Omega_p, w_p, gc_p,
                VY_p_Algebra, VF_p_Algebra, Omega_p_Algebra, w_p_Algebra, gt_pp,
                Omega_td_p, zeta_p, tau_p, Vgamma_p,
                mu_p, gt_p,
                w_o, VF_o, thetaNT, thetaT, VY_o_Algebra,
                Constraints, CovMatrix, Means, ModelExpectations, FitFunctionML)

    Model1 <- mxModel("Cascade_UniSEM_AMChange_Model1", Params, Example_Data_Mx)
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
