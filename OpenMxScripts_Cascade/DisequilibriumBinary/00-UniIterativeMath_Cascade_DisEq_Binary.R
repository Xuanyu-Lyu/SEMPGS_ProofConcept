## Univariate (scalar) iterative math for the SEM-PGS Cascade model in DISEQUILIBRIUM, standardized to a
## liability scale (parental VY == 1) for use with the binary/liability-threshold Cascade DisEq scripts.
##
## The base recursion (uniIterativeMath_Cascade_DisEq / checkCascadeDisEqConsistency /
## diffTraitOffspring_Cascade_DisEq) is copied verbatim from
## Disequilibrium/00-UniIterativeMath_Cascade_DisEq.R. The parental VT-only equilibrium does not depend
## on the mating mechanism, so the rescaling s = 1/sqrt(VY) is the same as in the original
## DisequilibriumBinary/00-UniIterativeMath_DisEq_Binary.R; the one AM event on gamma~ is homogeneous
## under it (mu scales by 1/s^2, gt/ht/itlo/itol unchanged).

uniIterativeMath_Cascade_DisEq <- function(delta, a, f, VE, am, AM_G = 1, AM_E = 1, k = .5, j = .5, gens = 100, tol = 1e-12){

    delta_td <- AM_G*delta
    a_td     <- AM_G*a

    ## ---- Phase 1: VT-only recursion (mu = 0) to convergence (identical to the original DisEq math) ----
    VY <- 2*delta^2*k + 2*a^2*j + VE
    w  <- v  <- 0
    VF <- 0
    Omega <- delta*k
    Gamma <- a*j

    traj <- list(c(gen = 0, VY = VY, Omega = Omega, Gamma = Gamma, w = w, v = v, VF = VF))
    converged <- FALSE
    it_used <- 0
    for (it in 1:gens){
        prev <- c(Omega, Gamma, VY, w, v, VF)

        w     <- 2*f*Omega
        v     <- 2*f*Gamma
        VF    <- 2*f^2*VY
        Omega <- delta*k + .5*w
        Gamma <- a*j + .5*v
        VY    <- 2*delta*Omega + 2*a*Gamma + w*delta + v*a + VF + VE
        traj[[it + 1]] <- c(gen = it, VY = VY, Omega = Omega, Gamma = Gamma, w = w, v = v, VF = VF)

        it_used <- it
        if (max(abs(c(Omega, Gamma, VY, w, v, VF) - prev)) < tol){
            converged <- TRUE
            break
        }
    }
    # re-close w, v, VF on the converged Omega/Gamma/VY (they satisfy the fixed point to < tol)
    w <- 2*f*Omega; v <- 2*f*Gamma; VF <- 2*f^2*VY

    # closed-form fixed points as a cross-check
    Omega_closed <- delta*k / (1 - f)
    Gamma_closed <- a*j / (1 - f)
    VY_closed    <- (2*delta*Omega_closed + 2*a*Gamma_closed +
                     2*f*Omega_closed*delta + 2*f*Gamma_closed*a + VE) / (1 - 2*f^2)
    stopifnot(abs(Omega - Omega_closed) < 1e-8,
              abs(Gamma - Gamma_closed) < 1e-8,
              abs(VY - VY_closed) < 1e-8)

    ## ---- Phase 2: the parents' one mating event, on gamma~ (PDF Part IV, Sections 15-16) ----
    Omega_td <- delta_td*k + .5*AM_E*w                # gc = hc = ic = 0 in the parents
    Gamma_td <- a_td*j + .5*AM_E*v
    zeta     <- delta_td*w + a_td*v + AM_E*VF
    tau      <- 2*a*Gamma_td + 2*delta*Omega_td + zeta + AM_E*VE
    Vgamma   <- 2*a_td*Gamma_td + 2*delta_td*Omega_td + AM_E*zeta + AM_E^2*VE

    mu   <- am / Vgamma
    gt   <- Omega_td^2 * mu
    ht   <- Gamma_td^2 * mu
    itlo <- Gamma_td*mu*Omega_td
    itol <- Omega_td*mu*Gamma_td
    gc <- hc <- ic <- 0                   # no AM in any previous generation

    # offspring generation (same expressions as the DisEq Cascade fitting scripts)
    w_o  <- 2*f*(Omega + tau*mu*Omega_td)
    v_o  <- 2*f*(Gamma + tau*mu*Gamma_td)
    VF_o <- 2*f^2*(VY + tau^2*mu)
    thetaNT <- delta*gt + a*itol + .5*w_o          # cov(Yo, NTp): ONE parental haplotype
    thetaT  <- delta*k + thetaNT                   # cov(Yo, Tp)
    VY_off  <- 2*delta^2*(k + gt) + 2*a^2*(j + ht) + 4*delta*a*itol + 2*delta*w_o + 2*a*v_o + VF_o + VE
    Yp_PGSm <- tau*mu*Omega_td
    Yp_LGSm <- tau*mu*Gamma_td
    Yp_Ym   <- tau*mu*tau
    Yp_Fm   <- tau*mu*zeta
    Fp_Fm   <- zeta*mu*zeta
    Yo_Yp   <- delta*Omega + a*Gamma + f*VY + tau*mu*(delta*Omega_td + a*Gamma_td + f*tau)

    # PDF Part IV, Sections 18-19: the offspring generation's own quantities as written in the PDF
    # (NOT used by the fitting scripts; reported so the simulation can adjudicate between the forms)
    w2  <- 2*f*(Omega + tau*mu*Omega_td)
    v2  <- 2*f*(Gamma + tau*mu*Gamma_td)
    VF2 <- 2*f^2*(VY + tau^2*mu)
    Omega2 <- delta*k + 2*delta*gt + 2*a*itol + .5*w2
    Gamma2 <- a*j + 2*a*ht + 2*delta*itol + .5*v2
    pdf_partIV <- list(w2 = w2, v2 = v2, VF2 = VF2, Omega2 = Omega2, Gamma2 = Gamma2,
                       VY2 = 2*a*Gamma2 + 2*delta*Omega2 + a*v2 + delta*w2 + VF2 + VE,
                       thetaNT2 = 4*delta*gt + 4*a*itol + w2,          # two-haplotype form
                       thetaT2  = 2*delta*k + 4*delta*gt + 4*a*itol + w2,
                       Y_Fo = f*(VY + tau^2*mu))

    list(delta = delta, a = a, f = f, VE = VE, am = am, AM_G = AM_G, AM_E = AM_E, k = k, j = j,
         delta_td = delta_td, a_td = a_td,
         VY = VY, VF = VF, mu = mu, Omega = Omega, Gamma = Gamma,
         Omega_td = Omega_td, Gamma_td = Gamma_td, zeta = zeta, tau = tau, Vgamma = Vgamma,
         gt = gt, ht = ht, gc = gc, hc = hc, itlo = itlo, itol = itol, ic = ic,
         w = w, v = v, w_o = w_o, v_o = v_o, VF_o = VF_o, VY_off = VY_off,
         thetaNT = thetaNT, thetaT = thetaT,
         Yp_PGSm = Yp_PGSm, Yp_LGSm = Yp_LGSm, Yp_Ym = Yp_Ym, Yp_Fm = Yp_Fm, Fp_Fm = Fp_Fm, Yo_Yp = Yo_Yp,
         Y_Fo = f*(VY + tau^2*mu),
         pdf_partIV = pdf_partIV,
         trajectory = as.data.frame(do.call(rbind, traj)),
         converged = converged, generations = it_used)
}

## Assert that the converged quantities satisfy the DisEq Cascade model constraints
## (the mxConstraint / mxAlgebra algebra of the disequilibrium Cascade fitting scripts).
checkCascadeDisEqConsistency <- function(q, tol = 1e-8){
    with(q, {
        stopifnot(
            abs(Omega - (delta*k + .5*w)) < tol,
            abs(Gamma - (a*j + .5*v))     < tol,
            abs(VY - (2*delta*Omega + 2*a*Gamma + w*delta + v*a + VF + VE)) < tol,
            abs(VF - 2*f^2*VY) < tol,
            abs(w - 2*f*Omega) < tol,
            abs(v - 2*f*Gamma) < tol,
            abs(Omega_td - (delta_td*k + .5*AM_E*w)) < tol,
            abs(Gamma_td - (a_td*j + .5*AM_E*v))     < tol,
            abs(tau - (2*a*Gamma_td + 2*delta*Omega_td + zeta + AM_E*VE)) < tol,
            abs(tau - (2*a_td*Gamma + 2*delta_td*Omega + AM_E*(a*v + delta*w + VF + VE))) < tol,
            abs(gt - Omega_td^2*mu) < tol,
            abs(ht - Gamma_td^2*mu) < tol,
            gc == 0, hc == 0, ic == 0,
            abs(mu - am/Vgamma) < tol,
            abs(w_o - (2*f*Omega + 2*f*tau*mu*Omega_td)) < tol,
            abs(VF_o - (2*f^2*VY + 2*f^2*tau^2*mu)) < tol,
            abs(thetaNT - (delta*gt + a*itol + .5*w_o)) < tol,
            abs(thetaT - (delta*k + thetaNT)) < tol,
            abs(VY_off - (2*delta^2*k + 2*delta^2*gt + 2*a^2*j + 2*a^2*ht + 4*delta*a*itol +
                          2*delta*w_o + 2*a*v_o + VF_o + VE)) < tol
        )
    })
    invisible(TRUE)
}

## Offspring-trait quantities for the DiffTrait designs; same algebra as the DisEq Cascade DiffTrait scripts
## (the offspring trait loads delta_o / a_o on the same haplotypes and receives the same Fo). Supply VE_o,
## or VY_o_target to solve VE_o so that the offspring-trait variance equals the target.
diffTraitOffspring_Cascade_DisEq <- function(parent, delta_o, a_o, VE_o = NULL, VY_o_target = NULL){
    VYo_terms <- with(parent, 2*delta_o^2*(k + gt) + 2*a_o^2*(j + ht) + 4*delta_o*a_o*itol +
                              2*delta_o*w_o + 2*a_o*v_o + VF_o)
    if (is.null(VE_o)) VE_o <- VY_o_target - VYo_terms
    stopifnot(VE_o > 0)
    thetaNT_o <- with(parent, delta_o*gt + a_o*itol + .5*w_o)     # per single haplotype
    thetaT_o  <- delta_o*parent$k + thetaNT_o
    Yo_Yp_o   <- with(parent, delta_o*Omega + a_o*Gamma + f*VY + tau*mu*(delta_o*Omega_td + a_o*Gamma_td + f*tau))
    list(delta_o = delta_o, a_o = a_o, VE_o = VE_o, VY_o = VYo_terms + VE_o,
         VYo_terms = VYo_terms, thetaNT_o = thetaNT_o, thetaT_o = thetaT_o, Yo_Yp_o = Yo_Yp_o)
}

## ---------------------------------------------------------------------------
## Binary/liability-standardized wrapper: same-trait / parent-trait disequilibrium
## point re-expressed on a scale where VY == 1 exactly.
## ---------------------------------------------------------------------------
uniIterativeMath_Cascade_DisEq_Binary <- function(delta, a, f, VE, am, AM_G = 1, AM_E = 1, k = .5, j = .5, gens = 100, tol = 1e-12){
    raw <- uniIterativeMath_Cascade_DisEq(delta = delta, a = a, f = f, VE = VE, am = am, AM_G = AM_G, AM_E = AM_E,
                                          k = k, j = j, gens = gens, tol = tol)
    stopifnot(raw$converged)
    s <- 1 / sqrt(raw$VY)
    std <- uniIterativeMath_Cascade_DisEq(delta = s*delta, a = s*a, f = f, VE = s^2*VE, am = am, AM_G = AM_G, AM_E = AM_E,
                                          k = k, j = j, gens = gens, tol = tol)
    stopifnot(std$converged, abs(std$VY - 1) < 1e-8)
    std$scale <- s
    std
}

## Given an already-standardized parent-trait disequilibrium point (VY_p == 1) and a choice of offspring
## delta_o/a_o, solve for the VE_o that makes the offspring's implied liability variance exactly 1 too
## (disequilibrium Cascade DiffTrait formula: the offspring's own variance, see diffTraitOffspring_Cascade_DisEq).
solve_VYo_Cascade_DisEq_Binary <- function(delta_o, a_o, parent, target = 1){
    diffTraitOffspring_Cascade_DisEq(parent, delta_o = delta_o, a_o = a_o, VY_o_target = target)
}
