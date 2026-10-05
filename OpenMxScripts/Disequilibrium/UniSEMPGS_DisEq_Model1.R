## This script is a function that fits Model 1 of Balbona, Kim & Keller (2021, Behav Genet 51:264-278): univariate
## SEM-PGS with vertical transmission and primary phenotypic assortative mating, where the PGS is assumed to explain
## ALL of the trait's heritability (no latent genetic score: a = 0, so no j, Gamma, v, h, i and no RDR constraint).
## Only the offspring phenotype and the four parental haplotypic PGS are used; the parental phenotypes are latent.
## This script is for trait that is not in equilibrium: VT is at its random-mating equilibrium and the parents mated
## assortatively for the first time (gc = 0), as in the paper's Supplement.
## It is the a = 0 special case of UniSEMPGS_DisEq_SameTrait_LatentYPYM.R, including its corrected offspring
## generation (../ConceptProof.md Section 6): thetaT/thetaNT are per single-haplotype column, and the offspring's F
## carries the AM feedback (w_o, VF_o), so the offspring have their own variance VY_o.
## (reference: perform_SEM_model1_dis() in ReferenceCode/VT_SEM_functions.R)
## Estimates are unbiased only to the degree that the PGS captures the heritability; see the LatentYPYM scripts.

fitUniSEMPGS_DisEq_Model1 <- function(data_path, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
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

    # Starting values from the sample moments, which put the first attempt near the solution: thetaT - thetaNT =
    # delta*k, Omega is close to thetaT (= Omega + delta*gt + f*Omega*mu*VY), and the parents' VY is close to the
    # offspring's variance. (With fixed starting values the first attempt can stop at a far-off point with status 1,
    # which mxTryHard accepts unless exhaustive = T.)
    S <- cov(Example_Data[, c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")], use = "pairwise.complete.obs")
    thetaT_s  <- mean(S["Yo1", c("Tp1", "Tm1")])
    thetaNT_s <- mean(S["Yo1", c("NTp1", "NTm1")])
    delta_s   <- (thetaT_s - thetaNT_s) / .5               # k = .5

    # --- Variances (parents) ---
    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=S["Yo1", "Yo1"], label="VY11", name="VY", lbound = .001)
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=S["Yo1", "Yo1"] / 2, label="VE11", name="VE", lbound = .001)

    # Equation for VY: 2*Omega*delta + delta*w + VF + VE
    VY_Algebra <- mxAlgebra(2 * Omega * delta + delta * w + VF_Algebra + VE, name="VY_Algebra")
    # Equation for VF: f*VY*f (assuming father/mother equality)
    VF_Algebra <- mxAlgebra(2 * f^2 * VY, name="VF_Algebra")

    VY_Constraint <- mxConstraint(VY == VY_Algebra, name='VY_Constraint')

    # --- Genetic effects & Path Coefficients ---
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=delta_s, label="delta11", name="delta")
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    # k is fixed at .5, which assumes the PGS is standardised in the base population. Under disequilibrium there is
    # no within-person haplotype covariance (gc = 0), so the haplotypic PGS variance is k in every generation, but
    # the offspring PGS variance is 2k + 2gt because its two haplotypes come from assorted parents. If the PGS is
    # instead standardised to variance 1 in the offspring sample, replace the fixed matrix above with
    # k <- mxAlgebra(.5 - gt, name = "k")
    # (scaling the haplotypic PGS to variance 1/2 leaves k = .5 unchanged here).

    # Omega = delta*k + 0.5*w
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=thetaT_s, label="Omega11", name="Omega")

    Omega_Algebra <- mxAlgebra(delta * k + 0.5 * w , name="Omega_Algebra")
    Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')

    # --- Assortative Mating (Non-Equilibrium) ---
    mu    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="mu11", name="mu")

    # gc = 0 (No AM in previous generation)
    gc    <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=0, label="gc11", name="gc")

    # gt is a function of mu and the current Omega
    gt_Algebra <- mxAlgebra(Omega * mu * Omega, name="gt_Algebra")
    gt <- mxAlgebra(gt_Algebra, name="gt")

    # --- Vertical Transmission ---
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="f11", name="f")

    # Equation: w = f*Omega + f*Omega (since parent traits/effects are equal)
    w_Algebra     <- mxAlgebra(2 * f * Omega, name="w_Algebra")
    w <- mxAlgebra(w_Algebra, name="w")

    # ---- Offspring generation ----
    # The offspring's F comes from AM-mated parents, so its covariance with the parents' haplotypes (w_o) and its
    # variance (VF_o) carry the cross-mate term mu*VY. (w and the VF algebra above belong to the parents, whose own
    # parents mated at random.)
    w_o  <- mxAlgebra(2 * f * Omega + 2 * f * VY * mu * Omega, name="w_o")
    VF_o <- mxAlgebra(2 * f^2 * VY + 2 * f^2 * VY^2 * mu, name="VF_o")
    # thetaNT / thetaT: covariance of Yo with ONE parental haplotype, i.e. with one data column (NTp1, Tp1, ...).
    # NOTE: in the papers (Balbona et al. 2021; Lyu et al.) theta is the SUM over the father's and the mother's
    # haplotype, thetaNT = cov(Yo, NTp + NTm) = 2*delta*gt + w_o: twice these.
    # The parents' own two haplotypes are uncorrelated (gc = 0: first AM event), so per haplotype
    #   cov(Yo, NTp) = delta*cov(Tm, NTp) + cov(Fo, NTp) = delta*gt + w_o/2
    thetaNT <- mxAlgebra(delta * gt + 0.5 * w_o, name="thetaNT")
    # cov(Yo, Tp) adds delta*var(Tp) = delta*k
    thetaT  <- mxAlgebra(delta * k + thetaNT, name="thetaT")
    # The offspring's own phenotypic variance: their two haplotypes covary by gt across parents but carry no
    # within-haplotype AM covariance yet, and their F has variance VF_o and covariance w_o with them
    VY_o_Algebra <- mxAlgebra(2 * delta^2 * k + 2 * delta^2 * gt + 2 * delta * w_o + VF_o + VE, name="VY_o_Algebra")

    # Expected covariances matrix (5x5: Yo1, Tp1, NTp1, Tm1, NTm1)
    CovMatrix <- mxAlgebra(rbind(
        cbind(VY_o_Algebra, thetaT,  thetaNT, thetaT,  thetaNT), # Yo1
        cbind(thetaT,  k+gc,    gc,      gt,      gt),      # Tp1
        cbind(thetaNT, gc,      k+gc,    gt,      gt),      # NTp1
        cbind(thetaT,  gt,      gt,      k+gc,    gc),      # Tm1
        cbind(thetaNT, gt,      gt,      gc,      k+gc)),   # NTm1
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
                w_o, VF_o, VY_o_Algebra,
                VY, VE, delta, k, Omega, mu, gc, f,
                VY_Algebra, VF_Algebra, Omega_Algebra, gt_Algebra, w_Algebra,
                gt, w,
                VY_Constraint, Omega_Constraint,
                thetaNT, thetaT,
                CovMatrix, Means, ModelExpectations, FitFunctionML)

    Model1 <- mxModel("UniSEM_NonEquilibrium_Model1", Params, Example_Data_Mx)
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
