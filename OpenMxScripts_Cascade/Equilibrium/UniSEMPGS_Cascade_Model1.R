## Model 1 of Balbona, Kim & Keller (2021), as a univariate SEM-PGS Cascade model.
## The PGS is assumed to capture all of the trait's heritability (no latent genetic score); only the offspring
## phenotype and the parents' PGS are observed.
## Equilibrium: assortative mating and vertical transmission have gone on for many generations.
## Spouses assort on a latent mating phenotype gamma~ = AM_G*delta*(T + NT) + AM_E*(F + E).
##
## Estimates (returned by summary(fit, verbose = TRUE)):
##   delta       effect of the haplotypic PGS on the phenotype
##   f           vertical transmission from each parent's phenotype to the offspring's environment F
##   mu          copath between spouses' mating phenotypes gamma~
##   VY, VE      phenotypic and residual variance
##   plus Omega, gt, gc, w, VF: covariance terms held to their model values by constraints.
##
## Input:
##   data_path   text file with a header row (read by data.table::fread) and the columns
##                 Yo1         offspring phenotype
##                 Tp1, NTp1   father's transmitted and non-transmitted haplotypic PGS
##                 Tm1, NTm1   mother's transmitted and non-transmitted haplotypic PGS
##               The file should hold only these columns.
##               Each haplotypic PGS is assumed to have variance k = .5 in the base population.
##   AM_G_value, AM_G_free, AM_E_value, AM_E_free  multipliers on the genetic (AM_G) and the F + E (AM_E) paths into
##               gamma~: (1, 1) primary phenotypic AM, (1, 0) genetic homogamy, (0, 1) social homogamy. At least one
##               must be fixed. Default: both fixed at 1. With latent parents the mechanism is not identified, so
##               fix both at the mechanism to be tested.
##   feaTol, optTol, extraTries, exhaustive, jitterMean, jitterVar  optimizer settings (NPSOL, mxTryHard)

fitUniSEMPGS_Cascade_Model1 <- function(data_path, AM_G_value = 1, AM_G_free = FALSE, AM_E_value = 1, AM_E_free = FALSE, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
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
    # Ensure data only contains: Yo1, Tp1, NTp1, Tm1, NTm1

    # Starting values from the sample moments. Model 1 is just-identified, so these put the first attempt near the
    # solution: thetaT - thetaNT = delta*k, Omega = thetaT, w = 2*(thetaNT - 2*delta*g). (With fixed starting values
    # the first attempt could stop at a far-off point with status 1, which mxTryHard accepts unless exhaustive = T.)
    S <- cov(Example_Data[, c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")], use = "pairwise.complete.obs")
    thetaT_s  <- mean(S["Yo1", c("Tp1", "Tm1")])
    thetaNT_s <- mean(S["Yo1", c("NTp1", "NTm1")])
    H         <- S[-1, -1]
    g_s       <- max(mean(H[upper.tri(H)]), 0)             # the six haplotype covariances (gc twice, gt four times)
    delta_s   <- (thetaT_s - thetaNT_s) / .5               # k = .5

    # Phenotypic and Residual Variance
    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=S["Yo1", "Yo1"], label="VY11", name="VY", lbound = .001)
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=S["Yo1", "Yo1"] / 2, label="VE11", name="VE", lbound = .001)

    # Scalar Algebra for Variances
    VY_Algebra <- mxAlgebra(2 * delta * Omega + w * delta + VF + VE, name="VY_Algebra")
    VF_Algebra <- mxAlgebra(2 * f^2 * VY + 2 * f^2 * tau^2 * mu, name="VF_Algebra")
    # Cascade: VF = 2f^2(VY + tau^2*mu) depends on tau = cov(Y, gamma~), which itself depends on VF,
    # so VF is a free parameter held to its algebra by a constraint (as VY is).
    VF    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.02 * S["Yo1", "Yo1"], label="VF11", name="VF")
    VF_Constraint <- mxConstraint(VF == VF_Algebra, name='VF_Constraint')
    VY_Constraint <- mxConstraint(VY == VY_Algebra, name='VY_Constraint')

    # Genetic effects
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=delta_s, label="delta11", name="delta")
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    # k is fixed at .5, which assumes the PGS is standardised in the base population. Empirical PGS are usually
    # scaled in the analysed sample, i.e. at equilibrium, where the haplotypic PGS variance is k + g and the full
    # PGS variance is 2k + 4g (g = gc = gt). To use one of those scalings, replace the fixed matrix above with
    #   haplotypic PGS scaled to variance 1/2 at equilibrium:  k <- mxAlgebra(.5 - gc, name = "k")
    #   full PGS standardised to variance 1 at equilibrium:    k <- mxAlgebra(.5 - 2 * gc, name = "k")
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=thetaT_s, label="Omega11", name="Omega")

    Omega_Algebra <- mxAlgebra(2 * delta * gc + delta * k + 0.5 * w , name="Omega_Algebra")
    Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')

    # Assortative mating effects; at equilibrium the cis (gc) and trans (gt) covariances are equal
    mu    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1,  label="mu11", name="mu")
    gt    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=g_s, label="gt11", name="gt")
    gc    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=g_s, label="gc11", name="gc")

    gt_Algebra <- mxAlgebra(Omega_td * mu * Omega_td, name="gt_Algebra")
    gc_Algebra <- mxAlgebra(gt, name="gc_Algebra")

    gt_constraint   <- mxConstraint(gt == gt_Algebra, name='gt_constraint')
    gc_constraint   <- mxConstraint(gc == gc_Algebra, name='gc_constraint')

    # Vertical transmission effects
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.1, label="f11", name="f")
    w     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=2 * (thetaNT_s - 2 * delta_s * g_s), label="w11", name="w")

    w_Algebra     <- mxAlgebra(2 * f * Omega + 2 * f * tau * mu * Omega_td, name="w_Algebra")
    w_constraint  <- mxConstraint(w == w_Algebra, name='w_constraint')

    # Between-people covariances (Excluding Parent Phenotypes)
    # thetaNT / thetaT: covariance of Yo with ONE parental haplotype, i.e. with one data column (NTp1, Tp1, ...).
    # In the papers theta is the SUM over the father's and the mother's haplotype: twice these.
    thetaNT <- mxAlgebra(2 * delta * gc + .5 * w, name="thetaNT")
    thetaT  <- mxAlgebra(delta * k + thetaNT, name="thetaT")

    # ---- Cascade AM: latent mating phenotype gamma~ ----
    # Tilde paths = multiplier x the path into Y (ETFD Cascade MVN scripts); see the header for AM_G / AM_E.
    AM_G <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_G_free, values=AM_G_value, label="AMGenMulti", name="AM_G")
    AM_E <- mxMatrix(type="Full", nrow=1, ncol=1, free=AM_E_free, values=AM_E_value, label="AMEnvMulti", name="AM_E")
    delta_td <- mxAlgebra(delta * AM_G, name="delta_td")
    # Shortcuts to gamma~ (PDF Part III Section 8, a = 0): Omega_td = cov(gamma~, [N]T), zeta = cov(gamma~, F),
    # tau = cov(Y, gamma~), Vgamma = var(gamma~)
    Omega_td <- mxAlgebra(2 * delta_td * gc + delta_td * k + 0.5 * AM_E * w, name="Omega_td")
    zeta     <- mxAlgebra(delta_td * w + AM_E * VF, name="zeta")
    tau      <- mxAlgebra(2 * delta * Omega_td + zeta + AM_E * VE, name="tau")
    Vgamma   <- mxAlgebra(2 * delta_td * Omega_td + AM_E * zeta + AM_E^2 * VE, name="Vgamma")

    # Expected covariances matrix (5x5: Yo1, Tp1, NTp1, Tm1, NTm1)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_Algebra, thetaT,  thetaNT, thetaT,  thetaNT), # Yo1
        cbind(thetaT,     k+gc,    gc,      gt,      gt),      # Tp1
        cbind(thetaNT,    gc,      k+gc,    gt,      gt),      # NTp1
        cbind(thetaT,     gt,      gt,      k+gc,    gc),      # Tm1
        cbind(thetaNT,    gt,      gt,      gc,      k+gc)),   # NTm1
        dimnames = list(c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1"),
                        c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")), name="expCov")

    Means <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = TRUE, values = 0,
        label = c("meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1"),
        dimnames = list(NULL, c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")),
        name = "expMeans")

    ModelExpectations <- mxExpectationNormal(covariance="expCov", means="expMeans",
                                            dimnames=c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1"))
    Example_Data_Mx   <- mxData(observed=Example_Data, type="raw")
    FitFunctionML     <- mxFitFunctionML()

    Params <- list(
                AM_G, AM_E, delta_td, Omega_td, zeta, tau, Vgamma, VF, VF_Constraint,
                VY, VE, delta, k, Omega, mu, gt, gc, f, w,
                VY_Algebra, VF_Algebra, Omega_Algebra, gt_Algebra, gc_Algebra, w_Algebra,
                VY_Constraint, Omega_Constraint, gt_constraint, gc_constraint, w_constraint,
                thetaNT, thetaT,
                CovMatrix, Means, ModelExpectations, FitFunctionML)

    Model1 <- mxModel("Cascade_UniSEM_Model1", Params, Example_Data_Mx)
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
