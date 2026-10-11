## Univariate (scalar) iterative math for the SEM-PGS Cascade model with a CHANGE IN ASSORTMENT STRENGTH
## (SEM_PGS_Cascade_model.pdf, Part VII: two-regime disequilibrium).
##
## Generations: _pp = the grandparents and every generation before them, _p = the observed parents, _o = their
## offspring. The grandparents mated on the latent mating phenotype
##     gamma~ = delta~*(T + NT) + a~*(LT + LNT) + 1~*F + 1~*E,  delta~ = AM_G*delta, a~ = AM_G*a, 1~ = AM_E
## with spousal gamma~ correlation `am_pp` for long enough to reach the Part III equilibrium. The parents are
## an ordinary member of that equilibrium population. They then mate once, on the same gamma~ (same AM_G, AM_E),
## with spousal correlation `am_p`. mu_pp = am_pp/Vgamma_p and mu_p = am_p/Vgamma_p.
##
## Phase 1 is the equilibrium recursion of Equilibrium/00-UniIterativeMath_Cascade.R at am_pp (copied, so this
## file stands alone). Its fixed point gives the parents' own quantities (_p), including their within-person
## (cis) haplotype covariances gc_p, hc_p, ic_p. These equal the covariances the grandparents' mating created
## between their two gametes: gt_pp = Omega~_p^2*mu_pp, ht_pp, it_pp.
##
## Phase 2 is the parents' one mating event at mu_p: gt_p = Omega~_p^2*mu_p, ht_p, it_p. The offspring side is
## derived directly from Yo = delta*(Tp + Tm) + a*(LTp + LTm) + Fo + Eo with Fo = f*(Yp + Ym), where Tp is the
## father's transmitted gamete (each locus from either of his haplotypes at random):
##  * var(Tp) = k + gc_p and cov(Tp, NTp) = gc_p: a parent's gamete carries the cis covariance the parent
##    inherited (half his within-haplotype LD plus half the covariance between his two haplotypes, both gc_p
##    at the mu_pp equilibrium). cov(Tp, Tm) = gt_p comes from the parents' own mating;
##  * w_o = cov(Fo, Tp + Tm) = 2f(Omega_p + tau_p*mu_p*Omega~_p), v_o likewise, VF_o = 2f^2(VY_p + tau_p^2*mu_p);
##  * per single-haplotype data column:
##      thetaNT = cov(Yo, NTp) = delta*(gc_p + gt_p) + a*(ic_p + it_p) + w_o/2,   thetaT = delta*k + thetaNT
##    (in the papers theta is the sum over the father's and the mother's haplotype: twice these);
##  * VY_o = 2delta^2(k + gc_p + gt_p) + 2a^2(j + hc_p + ht_p) + 4*delta*a*(ic_p + it_p)
##           + 2*delta*w_o + 2*a*v_o + VF_o + VE;
##  * Yo_Yp = delta*Omega_p + a*Gamma_p + f*VY_p + tau_p*mu_p*(delta*Omega~_p + a*Gamma~_p + f*tau_p) (PDF 38).
## am_pp = 0 gives gc_p = hc_p = ic_p = 0 and the Disequilibrium/ math exactly; am_pp = am_p gives the Part III
## equilibrium (gc_p = gt_p, ...). `pdf_partVII` reports the PDF's Section 38 generation-2 formulas for
## comparison only: they are exact only when am_pp = am_p (see ../PDF_Errata.md).

uniIterativeMath_Cascade_AMChange <- function(delta, a, f, VE, am_pp, am_p, AM_G = 1, AM_E = 1, k = .5, j = .5, gens = 200, tol = 1e-12){

    delta_td <- AM_G*delta
    a_td     <- AM_G*a

    ## ---- Phase 1: the grandparental equilibrium at am_pp (the parents' own quantities) ----
    # generation 0: a random-mating founder population with no VT history (the GeneEvolve starting state)
    gc <- hc <- ic <- 0
    w  <- v  <- VF <- 0

    traj <- vector("list", gens + 1)
    converged <- FALSE
    prev <- NULL
    for (g in 0:gens){
        Omega    <- 2*delta*gc + delta*k + .5*w + 2*a*ic
        Gamma    <- 2*a*hc + 2*delta*ic + a*j + .5*v
        VY       <- 2*delta*Omega + 2*a*Gamma + w*delta + v*a + VF + VE
        Omega_td <- 2*delta_td*gc + delta_td*k + .5*AM_E*w + 2*a_td*ic
        Gamma_td <- 2*a_td*hc + 2*delta_td*ic + a_td*j + .5*AM_E*v
        zeta     <- delta_td*w + a_td*v + AM_E*VF
        tau      <- 2*a*Gamma_td + 2*delta*Omega_td + zeta + AM_E*VE
        Vgamma   <- 2*a_td*Gamma_td + 2*delta_td*Omega_td + AM_E*zeta + AM_E^2*VE

        mu <- am_pp / Vgamma
        gt <- Omega_td^2*mu
        ht <- Gamma_td^2*mu
        it <- Omega_td*Gamma_td*mu

        traj[[g + 1]] <- c(gen = g, VY = VY, Omega = Omega, Gamma = Gamma, Omega_td = Omega_td, Gamma_td = Gamma_td,
                           zeta = zeta, tau = tau, Vgamma = Vgamma, mu = mu, gc = gc, hc = hc, ic = ic,
                           gt = gt, ht = ht, it = it, w = w, v = v, VF = VF)

        gc <- gt; hc <- ht; ic <- it
        w  <- 2*f*(Omega + tau*mu*Omega_td)
        v  <- 2*f*(Gamma + tau*mu*Gamma_td)
        VF <- 2*f^2*(VY + tau^2*mu)

        cur <- c(Omega, Gamma, VY, mu, gt, ht, w, v, VF, ic, tau, Vgamma)
        if (!is.null(prev) && max(abs(cur - prev)) < tol){ converged <- TRUE; break }
        prev <- cur
    }
    trajectory <- as.data.frame(do.call(rbind, traj[seq_len(g + 1)]))

    # the parents (_p): an ordinary member of the am_pp equilibrium; their cis covariances are the grandparents'
    # trans covariances (_pp)
    mu_pp <- mu
    gt_pp <- gt; ht_pp <- ht; it_pp <- it
    gc_p <- gc; hc_p <- hc; ic_p <- ic
    w_p <- w; v_p <- v; VF_p <- VF
    Omega_p <- Omega; Gamma_p <- Gamma; VY_p <- VY
    Omega_td_p <- Omega_td; Gamma_td_p <- Gamma_td; zeta_p <- zeta; tau_p <- tau; Vgamma_p <- Vgamma

    ## ---- Phase 2: the parents' one mating event at am_p (PDF Part VII, Sections 36-37) ----
    mu_p <- am_p / Vgamma_p
    gt_p <- Omega_td_p^2*mu_p
    ht_p <- Gamma_td_p^2*mu_p
    it_p <- Omega_td_p*Gamma_td_p*mu_p

    # offspring generation (same expressions as the AMChange Cascade fitting scripts)
    w_o  <- 2*f*(Omega_p + tau_p*mu_p*Omega_td_p)
    v_o  <- 2*f*(Gamma_p + tau_p*mu_p*Gamma_td_p)
    VF_o <- 2*f^2*(VY_p + tau_p^2*mu_p)
    thetaNT <- delta*(gc_p + gt_p) + a*(ic_p + it_p) + .5*w_o     # cov(Yo, NTp): ONE parental haplotype
    thetaT  <- delta*k + thetaNT                                   # cov(Yo, Tp)
    VY_o <- 2*delta^2*(k + gc_p + gt_p) + 2*a^2*(j + hc_p + ht_p) + 4*delta*a*(ic_p + it_p) +
            2*delta*w_o + 2*a*v_o + VF_o + VE
    Yp_PGSm <- tau_p*mu_p*Omega_td_p
    Yp_LGSm <- tau_p*mu_p*Gamma_td_p
    Yp_Ym   <- tau_p*mu_p*tau_p
    Yp_Fm   <- tau_p*mu_p*zeta_p
    Fp_Fm   <- zeta_p*mu_p*zeta_p
    Yo_Yp   <- delta*Omega_p + a*Gamma_p + f*VY_p + tau_p*mu_p*(delta*Omega_td_p + a*Gamma_td_p + f*tau_p)
    Y_Fo    <- f*(VY_p + tau_p^2*mu_p)                             # cov(Y*, Fo)

    # PDF Part VII, Section 38 as printed (NOT used by the fitting scripts; reported for comparison):
    # Omega2 / Gamma2 / thetaNT2 / VY2 use the parents' trans covariance g2 = gt_p twice and drop the
    # inherited gc_p. thetaNT2 and thetaT2 are the papers' two-haplotype sums.
    Omega2 <- delta*k + 2*delta*gt_p + 2*a*it_p + .5*w_o
    Gamma2 <- a*j + 2*a*ht_p + 2*delta*it_p + .5*v_o
    pdf_partVII <- list(w2 = w_o, v2 = v_o, VF2 = VF_o, Omega2 = Omega2, Gamma2 = Gamma2,
                        VY2 = 2*a*Gamma2 + 2*delta*Omega2 + a*v_o + delta*w_o + VF_o + VE,
                        thetaNT2 = 4*delta*gt_p + 4*a*it_p + w_o,
                        thetaT2  = 2*delta*k + 4*delta*gt_p + 4*a*it_p + w_o)

    list(delta = delta, a = a, f = f, VE = VE, am_pp = am_pp, am_p = am_p, AM_G = AM_G, AM_E = AM_E, k = k, j = j,
         delta_td = delta_td, a_td = a_td,
         mu_pp = mu_pp, gt_pp = gt_pp, ht_pp = ht_pp, it_pp = it_pp,
         VY_p = VY_p, VF_p = VF_p, Omega_p = Omega_p, Gamma_p = Gamma_p, w_p = w_p, v_p = v_p,
         gc_p = gc_p, hc_p = hc_p, ic_p = ic_p,
         Omega_td_p = Omega_td_p, Gamma_td_p = Gamma_td_p, zeta_p = zeta_p, tau_p = tau_p, Vgamma_p = Vgamma_p,
         mu_p = mu_p, gt_p = gt_p, ht_p = ht_p, it_p = it_p,
         w_o = w_o, v_o = v_o, VF_o = VF_o, VY_o = VY_o,
         thetaNT = thetaNT, thetaT = thetaT,
         Yp_PGSm = Yp_PGSm, Yp_LGSm = Yp_LGSm, Yp_Ym = Yp_Ym, Yp_Fm = Yp_Fm, Fp_Fm = Fp_Fm, Yo_Yp = Yo_Yp,
         Y_Fo = Y_Fo,
         pdf_partVII = pdf_partVII,
         trajectory = trajectory, converged = converged, generations = g)
}

## Assert that the converged quantities satisfy the AMChange Cascade model constraints
## (the mxConstraint / mxAlgebra algebra of the AMChange Cascade fitting scripts).
checkCascadeAMChangeConsistency <- function(q, tol = 1e-8){
    with(q, {
        stopifnot(
            # parents: the Part III equilibrium at mu_pp
            abs(Omega_p - (2*delta*gc_p + 2*a*ic_p + delta*k + .5*w_p)) < tol,
            abs(Gamma_p - (2*a*hc_p + 2*delta*ic_p + a*j + .5*v_p))     < tol,
            abs(VY_p - (2*delta*Omega_p + 2*a*Gamma_p + w_p*delta + v_p*a + VF_p + VE)) < tol,
            abs(Omega_td_p - (2*delta_td*gc_p + 2*a_td*ic_p + delta_td*k + .5*AM_E*w_p)) < tol,
            abs(Gamma_td_p - (2*a_td*hc_p + 2*delta_td*ic_p + a_td*j + .5*AM_E*v_p))     < tol,
            abs(zeta_p - (delta_td*w_p + a_td*v_p + AM_E*VF_p)) < tol,
            abs(tau_p - (2*a*Gamma_td_p + 2*delta*Omega_td_p + zeta_p + AM_E*VE)) < tol,
            abs(tau_p - (2*a_td*Gamma_p + 2*delta_td*Omega_p + AM_E*(a*v_p + delta*w_p + VF_p + VE))) < tol,
            abs(Vgamma_p - (2*a_td*Gamma_td_p + 2*delta_td*Omega_td_p + AM_E*zeta_p + AM_E^2*VE)) < tol,
            abs(VF_p - (2*f^2*VY_p + 2*f^2*tau_p^2*mu_pp)) < tol,
            abs(w_p - (2*f*Omega_p + 2*f*tau_p*mu_pp*Omega_td_p)) < tol,
            abs(v_p - (2*f*Gamma_p + 2*f*tau_p*mu_pp*Gamma_td_p)) < tol,
            abs(gc_p - Omega_td_p*mu_pp*Omega_td_p) < tol,
            abs(hc_p - Gamma_td_p*mu_pp*Gamma_td_p) < tol,
            abs(ic_p - Omega_td_p*mu_pp*Gamma_td_p) < tol,
            abs(mu_pp - am_pp/Vgamma_p) < tol,
            # the parents' mating at mu_p and the offspring
            abs(mu_p - am_p/Vgamma_p) < tol,
            abs(gt_p - Omega_td_p*mu_p*Omega_td_p) < tol,
            abs(ht_p - Gamma_td_p*mu_p*Gamma_td_p) < tol,
            abs(it_p - Omega_td_p*mu_p*Gamma_td_p) < tol,
            abs(w_o - (2*f*Omega_p + 2*f*tau_p*mu_p*Omega_td_p)) < tol,
            abs(v_o - (2*f*Gamma_p + 2*f*tau_p*mu_p*Gamma_td_p)) < tol,
            abs(VF_o - (2*f^2*VY_p + 2*f^2*tau_p^2*mu_p)) < tol,
            abs(thetaNT - (delta*gc_p + delta*gt_p + a*ic_p + a*it_p + .5*w_o)) < tol,
            abs(thetaT - (delta*k + thetaNT)) < tol,
            abs(VY_o - (2*delta^2*k + 2*delta^2*gc_p + 2*delta^2*gt_p + 2*a^2*j + 2*a^2*hc_p + 2*a^2*ht_p +
                        4*delta*a*ic_p + 4*delta*a*it_p + 2*delta*w_o + 2*a*v_o + VF_o + VE)) < tol
        )
    })
    invisible(TRUE)
}

## Offspring-trait quantities for the DiffTrait designs; same algebra as the AMChange Cascade DiffTrait scripts
## (the offspring trait loads delta_o / a_o on the same haplotypes and receives the same Fo). Supply VE_o,
## or VY_o_target to solve VE_o so that the offspring-trait variance equals the target.
diffTraitOffspring_Cascade_AMChange <- function(parent, delta_o, a_o, VE_o = NULL, VY_o_target = NULL){
    VYo_terms <- with(parent, 2*delta_o^2*(k + gc_p + gt_p) + 2*a_o^2*(j + hc_p + ht_p) +
                              4*delta_o*a_o*(ic_p + it_p) + 2*delta_o*w_o + 2*a_o*v_o + VF_o)
    if (is.null(VE_o)) VE_o <- VY_o_target - VYo_terms
    stopifnot(VE_o > 0)
    thetaNT_o <- with(parent, delta_o*(gc_p + gt_p) + a_o*(ic_p + it_p) + .5*w_o)     # per single haplotype
    thetaT_o  <- delta_o*parent$k + thetaNT_o
    Yo_Yp_o   <- with(parent, delta_o*Omega_p + a_o*Gamma_p + f*VY_p +
                              tau_p*mu_p*(delta_o*Omega_td_p + a_o*Gamma_td_p + f*tau_p))
    list(delta_o = delta_o, a_o = a_o, VE_o = VE_o, VY_o = VYo_terms + VE_o,
         VYo_terms = VYo_terms, thetaNT_o = thetaNT_o, thetaT_o = thetaT_o, Yo_Yp_o = Yo_Yp_o)
}
