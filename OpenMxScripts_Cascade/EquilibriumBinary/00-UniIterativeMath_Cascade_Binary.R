## Univariate (scalar) iterative math for the SEM-PGS Cascade model at equilibrium, standardized to a
## liability scale (VY == 1) for use with the binary/liability-threshold Cascade fitting scripts.
##
## The base recursion (uniIterativeMath_Cascade / checkCascadeEquilibriumConsistency /
## diffTraitOffspring_Cascade) is copied verbatim from Equilibrium/00-UniIterativeMath_Cascade.R.
## On top of it, uniIterativeMath_Cascade_Binary() uses the same homogeneity as the original
## EquilibriumBinary/00-UniIterativeMath_Binary.R: with f, am, AM_G and AM_E held fixed, the Cascade
## fixed point is homogeneous under (delta, a, VE) -> (s*delta, s*a, s^2*VE), s = 1/sqrt(VY):
## gamma~ scales by s along with Y, so mu = am/Vgamma scales by 1/s^2 and every PGS-PGS covariance
## (gt/ht/gc/hc/ic = Omega~^2*mu, ...) is unchanged. Re-running the recursion with the rescaled inputs
## therefore lands exactly on VY == 1.

uniIterativeMath_Cascade <- function(delta, a, f, VE, am, AM_G = 1, AM_E = 1, k = .5, j = .5, gens = 200, tol = 1e-12){

    delta_td <- AM_G*delta
    a_td     <- AM_G*a

    # generation 0: no AM-induced haplotype covariance, no VT
    gc <- hc <- ic <- 0
    w  <- v  <- VF <- 0

    traj <- vector("list", gens + 1)
    converged <- FALSE
    prev <- NULL
    for (g in 0:gens){
        # this generation's own quantities (PDF Part III, Sections 8 and 11)
        Omega    <- 2*delta*gc + delta*k + .5*w + 2*a*ic
        Gamma    <- 2*a*hc + 2*delta*ic + a*j + .5*v
        VY       <- 2*delta*Omega + 2*a*Gamma + w*delta + v*a + VF + VE
        Omega_td <- 2*delta_td*gc + delta_td*k + .5*AM_E*w + 2*a_td*ic   # cov(gamma~, [N]T)
        Gamma_td <- 2*a_td*hc + 2*delta_td*ic + a_td*j + .5*AM_E*v       # cov(gamma~, L[N]T)
        zeta     <- delta_td*w + a_td*v + AM_E*VF                         # cov(gamma~, F)
        tau      <- 2*a*Gamma_td + 2*delta*Omega_td + zeta + AM_E*VE      # cov(Y, gamma~)
        Vgamma   <- 2*a_td*Gamma_td + 2*delta_td*Omega_td + AM_E*zeta + AM_E^2*VE

        # this generation mates on gamma~ (Section 9)
        mu <- am / Vgamma
        gt <- Omega_td^2*mu
        ht <- Gamma_td^2*mu
        it <- Omega_td*Gamma_td*mu

        traj[[g + 1]] <- c(gen = g, VY = VY, Omega = Omega, Gamma = Gamma, Omega_td = Omega_td, Gamma_td = Gamma_td,
                           zeta = zeta, tau = tau, Vgamma = Vgamma, mu = mu, gc = gc, hc = hc, ic = ic,
                           gt = gt, ht = ht, it = it, w = w, v = v, VF = VF,
                           Yp_Ym = tau^2*mu, Yp_PGSm = tau*mu*Omega_td)

        # offspring generation's inputs: cis = parents' trans (equilibrium closure), VT through Y (Section 10)
        gc <- gt; hc <- ht; ic <- it
        w  <- 2*f*(Omega + tau*mu*Omega_td)
        v  <- 2*f*(Gamma + tau*mu*Gamma_td)
        VF <- 2*f^2*(VY + tau^2*mu)

        cur <- c(Omega, Gamma, VY, mu, gt, ht, w, v, VF, ic, tau, Vgamma)
        if (!is.null(prev) && max(abs(cur - prev)) < tol){ converged <- TRUE; break }
        prev <- cur
    }
    trajectory <- as.data.frame(do.call(rbind, traj[seq_len(g + 1)]))

    # derived covariance pieces (same expressions as the equilibrium Cascade fitting scripts)
    thetaNT <- 2*delta*gc + 2*a*ic + .5*w
    thetaT  <- delta*k + thetaNT
    Yp_PGSm <- tau*mu*Omega_td                  # cov(Yp, [N]Tm) = tau*mu*Omega~
    Yp_LGSm <- tau*mu*Gamma_td                  # cov(Yp, L[N]Tm)
    Yp_Ym   <- tau*mu*tau                       # cov(Yp, Ym)
    Yp_Fm   <- tau*mu*zeta                      # cov(Yp, Fm)
    Fp_Fm   <- zeta*mu*zeta                     # cov(Fp, Fm)
    Yo_Yp   <- delta*Omega + a*Gamma + f*VY + tau*mu*(delta*Omega_td + a*Gamma_td + f*tau)
    Y_Fo    <- f*(VY + tau^2*mu)                # cov(Y*, Fo)

    list(delta = delta, a = a, f = f, VE = VE, am = am, AM_G = AM_G, AM_E = AM_E, k = k, j = j,
         delta_td = delta_td, a_td = a_td,
         VY = VY, VF = VF, mu = mu, Omega = Omega, Gamma = Gamma,
         Omega_td = Omega_td, Gamma_td = Gamma_td, zeta = zeta, tau = tau, Vgamma = Vgamma,
         gt = gt, ht = ht, gc = gc, hc = hc, itlo = it, itol = it, ic = ic,
         w = w, v = v,
         thetaNT = thetaNT, thetaT = thetaT,
         Yp_PGSm = Yp_PGSm, Yp_LGSm = Yp_LGSm, Yp_Ym = Yp_Ym, Yp_Fm = Yp_Fm, Fp_Fm = Fp_Fm,
         Yo_Yp = Yo_Yp, Y_Fo = Y_Fo,
         trajectory = trajectory, converged = converged, generations = g)
}

## Assert that the converged quantities satisfy the equilibrium Cascade model constraints
## (the mxConstraint algebra of the equilibrium Cascade fitting scripts).
checkCascadeEquilibriumConsistency <- function(q, tol = 1e-8){
    with(q, {
        stopifnot(
            abs(Omega - (2*delta*gc + delta*k + .5*w + 2*a*ic)) < tol,
            abs(Gamma - (2*a*hc + 2*delta*ic + a*j + .5*v))     < tol,
            abs(VY - (2*delta*Omega + 2*a*Gamma + w*delta + v*a + VF + VE)) < tol,
            abs(Omega_td - (2*delta_td*gc + delta_td*k + .5*AM_E*w + 2*a_td*ic)) < tol,
            abs(Gamma_td - (2*a_td*hc + 2*delta_td*ic + a_td*j + .5*AM_E*v))     < tol,
            abs(zeta - (delta_td*w + a_td*v + AM_E*VF)) < tol,
            abs(tau - (2*a*Gamma_td + 2*delta*Omega_td + zeta + AM_E*VE)) < tol,
            # tau's second form (PDF Part III, Section 8) must agree with the first
            abs(tau - (2*a_td*Gamma + 2*delta_td*Omega + AM_E*(a*v + delta*w + VF + VE))) < tol,
            abs(Vgamma - (2*a_td*Gamma_td + 2*delta_td*Omega_td + AM_E*zeta + AM_E^2*VE)) < tol,
            abs(VF - (2*f^2*VY + 2*f^2*tau^2*mu)) < tol,
            abs(w - (2*f*Omega + 2*f*tau*mu*Omega_td)) < tol,
            abs(v - (2*f*Gamma + 2*f*tau*mu*Gamma_td)) < tol,
            abs(gt - Omega_td^2*mu) < tol,
            abs(gc - gt) < tol,
            abs(ht - Gamma_td^2*mu) < tol,
            abs(hc - ht) < tol,
            abs(ic - Omega_td*mu*Gamma_td) < tol,
            abs(mu - am/Vgamma) < tol
        )
    })
    invisible(TRUE)
}

## Offspring-trait quantities for the DiffTrait designs (parent trait = the trait spouses assort and
## transmit on; offspring trait loads delta_o/a_o on the same PGS/LGS and receives the same F).
## Same algebra as the equilibrium Cascade DiffTrait scripts. Supply VE_o, or VY_o_target to solve VE_o
## so that the offspring-trait variance equals the target (the DiffTrait_Latent models need VY_o == VY_p).
diffTraitOffspring_Cascade <- function(parent, delta_o, a_o, VE_o = NULL, VY_o_target = NULL){
    VYo_terms <- with(parent,
        2*delta_o^2*k + 4*delta_o^2*gc + 2*a_o^2*j + 4*a_o^2*hc +
        8*delta_o*ic*a_o + 2*delta_o*w + 2*a_o*v + VF)
    if (is.null(VE_o)) VE_o <- VY_o_target - VYo_terms
    stopifnot(VE_o > 0)
    thetaNT_o <- with(parent, 2*delta_o*gc + 2*a_o*ic + .5*w)
    thetaT_o  <- parent$k*delta_o + thetaNT_o
    Yo_Yp_o   <- with(parent, delta_o*Omega + a_o*Gamma + f*VY + tau*mu*(delta_o*Omega_td + a_o*Gamma_td + f*tau))
    list(delta_o = delta_o, a_o = a_o, VE_o = VE_o, VY_o = VYo_terms + VE_o,
         VYo_terms = VYo_terms, thetaNT_o = thetaNT_o, thetaT_o = thetaT_o, Yo_Yp_o = Yo_Yp_o)
}

## ---------------------------------------------------------------------------
## Binary/liability-standardized wrapper: same-trait / parent-trait equilibrium
## point re-expressed on a scale where VY == 1 exactly.
## ---------------------------------------------------------------------------
uniIterativeMath_Cascade_Binary <- function(delta, a, f, VE, am, AM_G = 1, AM_E = 1, k = .5, j = .5, gens = 200, tol = 1e-12){
    raw <- uniIterativeMath_Cascade(delta = delta, a = a, f = f, VE = VE, am = am, AM_G = AM_G, AM_E = AM_E,
                                    k = k, j = j, gens = gens, tol = tol)
    stopifnot(raw$converged)
    s <- 1 / sqrt(raw$VY)
    std <- uniIterativeMath_Cascade(delta = s*delta, a = s*a, f = f, VE = s^2*VE, am = am, AM_G = AM_G, AM_E = AM_E,
                                    k = k, j = j, gens = gens, tol = tol)
    stopifnot(std$converged, abs(std$VY - 1) < 1e-8)
    std$scale <- s
    std
}

## Given an already-standardized parent-trait equilibrium point (VY_p == 1) and a choice of offspring
## delta_o/a_o, solve for the VE_o that makes the offspring's implied liability variance exactly 1 too
## (equilibrium Cascade DiffTrait formula; VF is the Cascade VF = 2f^2(VY + tau^2*mu)).
solve_VYo_Cascade_Equilibrium_Binary <- function(delta_o, a_o, parent, target = 1){
    diffTraitOffspring_Cascade(parent, delta_o = delta_o, a_o = a_o, VY_o_target = target)
}
