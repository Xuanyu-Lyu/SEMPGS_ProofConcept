## Univariate (scalar) iterative math for the SEM-PGS model in DISEQUILIBRIUM.
## Matches the disequilibrium model equations: vertical transmission is iterated to its
## random-mating equilibrium (mu = 0, so gc = hc = ic = 0 and w = 2*f*Omega with no mu terms),
## then assortative mating happens for the FIRST time in the observed parents' generation.
## This exactly satisfies the DisEq scripts' self-referential constraints
## (Omega = delta*k + .5*w with w = 2*f*Omega; VF = 2*f^2*VY).
##
## Offspring generation (the children of that first AM event), derived from
## Yo = delta*(Tp + Tm) + a*(LTp + LTm) + Fo + Eo with Fo = f*(Yp + Ym):
##  * the parents' own two haplotypes are uncorrelated (gc = hc = ic = 0), so the offspring's two haplotypes
##    covary only ACROSS parents (gt, ht, it) and each still has variance k (j);
##  * Fo comes from AM-mated parents: w_o = cov(Fo, Tp + Tm) = 2f*Omega*(1 + mu*VY), v_o likewise, and
##    VF_o = 2f^2*VY*(1 + mu*VY) (the parents' own w, v, VF have no mu term);
##  * thetaNT / thetaT are per single-haplotype column (the data columns NTp1, Tp1, ...):
##    cov(Yo, NTp) = delta*gt + a*it + w_o/2 and cov(Yo, Tp) = delta*k + cov(Yo, NTp). In the papers theta is
##    the SUM over the father's and the mother's haplotype, cov(Yo, NTp + NTm): twice these;
##  * the offspring's own variance VY_off = 2delta^2(k + gt) + 2a^2(j + ht) + 4*delta*a*it + 2*delta*w_o
##    + 2*a*v_o + VF_o + VE;
##  * cov(Yo, Yp) = (delta*Omega + a*Gamma + f*VY)*(1 + mu*VY), including f*cov(Ym, Yp) = f*mu*VY^2.
## (Before Sep 2026 these used the papers' two-haplotype theta on single-haplotype columns, the parents' VY
## for the offspring, the parents' w/VF, and omitted f*mu*VY^2 from Yo_Yp; see ../ConceptProof.md, Section 6.)

uniIterativeMath_DisEq <- function(delta, a, f, VE, am, k = .5, j = .5, gens = 100, tol = 1e-12){

    ## ---- Phase 1: VT-only recursion (mu = 0) to convergence ----
    VY <- 2*delta^2*k + 2*a^2*j + VE
    w  <- v  <- 0
    VF <- 2*f^2*VY
    Omega <- delta*k
    Gamma <- a*j

    converged <- FALSE
    it_used <- 1
    for (it in 2:gens){
        prev <- c(Omega, Gamma, VY, w, v, VF)

        Omega <- delta*k + .5*w
        Gamma <- a*j + .5*v
        VY    <- 2*delta*Omega + 2*a*Gamma + w*delta + v*a + VF + VE
        w     <- 2*f*Omega
        v     <- 2*f*Gamma
        VF    <- 2*f^2*VY

        it_used <- it
        if (max(abs(c(Omega, Gamma, VY, w, v, VF) - prev)) < tol){
            converged <- TRUE
            break
        }
    }

    # closed-form fixed points as a cross-check
    Omega_closed <- delta*k / (1 - f)
    Gamma_closed <- a*j / (1 - f)
    VY_closed    <- (2*delta*Omega_closed + 2*a*Gamma_closed +
                     2*f*Omega_closed*delta + 2*f*Gamma_closed*a + VE) / (1 - 2*f^2)
    stopifnot(abs(Omega - Omega_closed) < 1e-8,
              abs(Gamma - Gamma_closed) < 1e-8,
              abs(VY - VY_closed) < 1e-8)

    ## ---- Phase 2: one generation of assortative mating ----
    mu   <- am / VY
    gt   <- Omega^2 * mu
    ht   <- Gamma^2 * mu
    itlo <- Gamma*mu*Omega
    itol <- Omega*mu*Gamma
    gc <- hc <- ic <- 0                   # no AM in any previous generation

    # offspring generation (same expressions as the DisEq fitting scripts; see the header)
    w_o  <- 2*f*Omega + 2*f*VY*mu*Omega
    v_o  <- 2*f*Gamma + 2*f*VY*mu*Gamma
    VF_o <- 2*f^2*VY + 2*f^2*mu*VY^2
    thetaNT <- delta*gt + a*itol + .5*w_o            # cov(Yo, NTp): ONE parental haplotype
    thetaT  <- delta*k + thetaNT                     # cov(Yo, Tp)
    VY_off  <- 2*delta^2*(k + gt) + 2*a^2*(j + ht) + 4*delta*a*itol + 2*delta*w_o + 2*a*v_o + VF_o + VE
    Yp_PGSm <- VY*mu*Omega
    Yp_Ym   <- mu*VY^2
    Yo_Yp   <- (delta*Omega + a*Gamma + f*VY)*(1 + mu*VY)

    list(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j,
         VY = VY, mu = mu, Omega = Omega, Gamma = Gamma,
         gt = gt, ht = ht, gc = gc, hc = hc, itlo = itlo, itol = itol, ic = ic,
         w = w, v = v, VF = VF, w_o = w_o, v_o = v_o, VF_o = VF_o, VY_off = VY_off,
         thetaNT = thetaNT, thetaT = thetaT,
         Yp_PGSm = Yp_PGSm, Yp_Ym = Yp_Ym, Yo_Yp = Yo_Yp,
         converged = converged, generations = it_used)
}

## Assert that the converged quantities satisfy the DisEq model constraints
## (the mxConstraint algebra of the disequilibrium fitting scripts).
checkDisEqConsistency <- function(q, tol = 1e-8){
    with(q, {
        stopifnot(
            abs(Omega - (delta*k + .5*w)) < tol,
            abs(Gamma - (a*j + .5*v))     < tol,
            abs(VY - (2*delta*Omega + 2*a*Gamma + w*delta + v*a + VF + VE)) < tol,
            abs(VF - 2*f^2*VY) < tol,
            abs(w - 2*f*Omega) < tol,
            abs(v - 2*f*Gamma) < tol,
            abs(gt - Omega^2*mu) < tol,
            abs(ht - Gamma^2*mu) < tol,
            gc == 0, hc == 0, ic == 0,
            abs(mu - am/VY) < tol,
            abs(w_o - (2*f*Omega + 2*f*VY*mu*Omega)) < tol,
            abs(VF_o - (2*f^2*VY + 2*f^2*mu*VY^2)) < tol,
            abs(thetaNT - (delta*gt + a*itol + .5*w_o)) < tol,
            abs(thetaT - (delta*k + thetaNT)) < tol,
            abs(VY_off - (2*delta^2*k + 2*delta^2*gt + 2*a^2*j + 2*a^2*ht + 4*delta*a*itol +
                          2*delta*w_o + 2*a*v_o + VF_o + VE)) < tol
        )
    })
    invisible(TRUE)
}

## Offspring-trait quantities for the DiffTrait designs (same algebra as the DisEq DiffTrait scripts): the offspring
## trait loads delta_o / a_o on the same haplotypes and receives the same Fo. Supply VE_o, or VY_o_target to solve
## VE_o so that the offspring-trait variance equals the target (the DiffTrait_Latent models need VY_o == VY_p).
diffTraitOffspring_DisEq <- function(parent, delta_o, a_o, VE_o = NULL, VY_o_target = NULL){
    VYo_terms <- with(parent, 2*delta_o^2*(k + gt) + 2*a_o^2*(j + ht) + 4*delta_o*a_o*itol +
                              2*delta_o*w_o + 2*a_o*v_o + VF_o)
    if (is.null(VE_o)) VE_o <- VY_o_target - VYo_terms
    stopifnot(VE_o > 0)
    thetaNT_o <- with(parent, delta_o*gt + a_o*itol + .5*w_o)     # per single haplotype
    thetaT_o  <- delta_o*parent$k + thetaNT_o
    Yo_Yp_o   <- with(parent, (delta_o*Omega + a_o*Gamma + f*VY)*(1 + mu*VY))
    list(delta_o = delta_o, a_o = a_o, VE_o = VE_o, VY_o = VYo_terms + VE_o,
         VYo_terms = VYo_terms, thetaNT_o = thetaNT_o, thetaT_o = thetaT_o, Yo_Yp_o = Yo_Yp_o)
}
