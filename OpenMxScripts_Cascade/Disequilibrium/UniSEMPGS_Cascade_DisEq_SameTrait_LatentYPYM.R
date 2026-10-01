## SEM-PGS CASCADE version of OpenMxScripts/Disequilibrium/UniSEMPGS_DisEq_SameTrait_LatentYPYM.R: the same model, except that spouses
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
## This script is a function that fits a version univariate SEM-PGS where parent and offspring have the same trait but the parental phenotypes are latent.
## This script is for trait that is not in equilibrium
## To identify this model, the optimal solution is to use RDR to find a for both parents and offspring.

fitUniSEMPGS_Cascade_DisEq_SameTrait_LatentParents <- function(data_path, h2_RDR, AM_G_value = 1, AM_G_free = FALSE, AM_E_value = 1, AM_E_free = FALSE, feaTol = 1e-6, optTol = 1e-8, jitterMean = .5, jitterVar = .1, extraTries = 30, exhaustive = F){
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

    # --- Variances ---
    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=1.5, label="VY11", name="VY", lbound = .001)
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.5, label="VE11", name="VE", lbound = .001)

    # Equation for VY: 2*Omega*delta + 2*Gamma*a + delta*w + a*v + VF + VE
    VY_Algebra <- mxAlgebra(2 * Omega * delta + 2 * Gamma * a + delta * w + a * v + VF_Algebra + VE, name="VY_Algebra")
    # Equation for VF: f*VY*f (assuming father/mother equality)
    VF_Algebra <- mxAlgebra(2 * f^2 * VY, name="VF_Algebra")

    VY_Constraint <- mxConstraint(VY == VY_Algebra, name='VY_Constraint')

    # --- Genetic effects & Path Coefficients ---
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.25, label="delta11", name="delta")
    a     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.7, label="a11", name="a", lbound = .001)
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k11", name="k")
    j     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="j11", name="j")
    # k and j are fixed at .5, which assumes the PGS (and, by convention, the latent LGS) is standardised
    # in the base population. Under disequilibrium there is no within-person haplotype covariance
    # (gc = hc = 0), so the haplotypic PGS variance is k in every generation, but the offspring PGS
    # variance is 2k + 2gt because its two haplotypes come from assorted parents. If the PGS is instead
    # standardised to variance 1 in the offspring sample, replace the two fixed matrices above with
    # k <- mxAlgebra(.5 - gt, name = "k")
    # j <- mxAlgebra(.5 - ht, name = "j")
    # (scaling the haplotypic PGS to variance 1/2 leaves k = j = .5 unchanged here).

    # Omega = delta*k + 0.5*w
    # Gamma = a*j + 0.5*v
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="Omega11", name="Omega")
    Gamma <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="Gamma11", name="Gamma")

    Omega_Algebra <- mxAlgebra(delta * k + 0.5 * w , name="Omega_Algebra")
    Gamma_Algebra <- mxAlgebra(a * j + 0.5 * v, name="Gamma_Algebra")

    Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
    Gamma_Constraint <- mxConstraint(Gamma == Gamma_Algebra, name='Gamma_Constraint')

    # RDR Constraint to identify 'a'
    # Current notation: sigma2 is VY, e^2 is VE
    rdr_left  <- mxAlgebra((2*a^2*j + 2*delta^2*k) * (2*a^2*j + 2*delta^2*k + VE), name="rdr_left")
    h2mat     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=h2_RDR, name="h2mat")
    rdr_right <- mxAlgebra(h2mat * VY, name="rdr_right")
    rdrCon    <- mxConstraint(rdr_left == rdr_right, name="rdrCon")

    # --- Assortative Mating (Non-Equilibrium) ---
    mu    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="mu11", name="mu")

    # gc = hc = ic = 0 (No AM in previous generation)
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
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.15, label="f11", name="f")

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
                w_o, v_o, VF_o, VY_o_Algebra,
                AM_G, AM_E, delta_td, a_td, Omega_td, Gamma_td, zeta, tau, Vgamma,
                VY, VE, delta, a, k, j, Omega, Gamma, mu, gc, hc, f,
                VY_Algebra, VF_Algebra, Omega_Algebra, Gamma_Algebra,
                gt_Algebra, ht_Algebra, ic_Algebra, w_Algebra, v_Algebra,
                gt, ht, ic, w, v,
                VY_Constraint, Omega_Constraint, Gamma_Constraint,
                h2mat, rdr_left, rdr_right, rdrCon, # Identification for 'a'
                thetaNT, thetaT,
                CovMatrix, Means, ModelExpectations, FitFunctionML)

    Model1 <- mxModel("Cascade_UniSEM_NonEquilibrium_LatentParents", Params, Example_Data_Mx)
    fitModel1 <- mxTryHard(Model1, extraTries = extraTries, OKstatuscodes = c(0,1),
                          intervals=T, silent=T, exhaustive = exhaustive,
                          jitterDistrib = "rnorm", loc=jitterMean, scale = jitterVar)

    return(summary(fitModel1, verbose = TRUE))
}
