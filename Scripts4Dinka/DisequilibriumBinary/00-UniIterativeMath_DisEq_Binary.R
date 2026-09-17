## Univariate (scalar) iterative math for the SEM-PGS model in DISEQUILIBRIUM, standardized
## to a liability scale (VY == 1) for use with the binary/liability-threshold fitting scripts.
##
## The base recursion (uniIterativeMath_DisEq / checkDisEqConsistency) is copied verbatim from
## Disequilibrium/00-UniIterativeMath_DisEq.R. uniIterativeMath_DisEq_Binary() re-runs it with
## (delta, a, VE) -> (s*delta, s*a, s^2*VE), s = 1/sqrt(VY), f and am held fixed, which drives the
## fixed point to VY == 1 exactly (see EquilibriumBinary/00-UniIterativeMath_Binary.R for the
## derivation and numerical verification; the same homogeneity holds for the VT-then-AM recursion).

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

    # derived covariance pieces (same expressions as the DisEq fitting scripts)
    thetaNT <- a*itlo + a*itol + 2*delta*gt + w
    thetaT  <- 2*delta*k + thetaNT
    Yp_PGSm <- VY*mu*Omega
    Yp_Ym   <- mu*VY^2
    Yo_Yp   <- delta*Omega + a*Gamma + f*VY + (delta*Omega + a*Gamma)*mu*VY

    list(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j,
         VY = VY, mu = mu, Omega = Omega, Gamma = Gamma,
         gt = gt, ht = ht, gc = gc, hc = hc, itlo = itlo, itol = itol, ic = ic,
         w = w, v = v, VF = VF,
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
            abs(mu - am/VY) < tol
        )
    })
    invisible(TRUE)
}

## ---------------------------------------------------------------------------
## Binary/liability-standardized wrapper: same-trait / parent-trait disequilibrium
## point re-expressed on a scale where VY == 1 exactly.
## ---------------------------------------------------------------------------
uniIterativeMath_DisEq_Binary <- function(delta, a, f, VE, am, k = .5, j = .5, gens = 100, tol = 1e-12){
    raw <- uniIterativeMath_DisEq(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j, gens = gens, tol = tol)
    stopifnot(raw$converged)
    s <- 1 / sqrt(raw$VY)
    std <- uniIterativeMath_DisEq(delta = s*delta, a = s*a, f = f, VE = s^2*VE, am = am, k = k, j = j, gens = gens, tol = tol)
    stopifnot(std$converged, abs(std$VY - 1) < 1e-8)
    std$scale <- s
    std
}

## Given an already-standardized parent-trait disequilibrium point (VY_p == 1) and a choice of
## offspring generation's delta_o/a_o, solve for the VE_o that makes the offspring's implied
## liability variance equal exactly 1 too (disequilibrium DiffTrait formula: gc=hc=ic=0).
solve_VYo_DisEq_Binary <- function(delta_o, a_o, k, j, parent, target = 1){
    VYo_terms <- with(parent, 2*delta_o^2*k + 2*a_o^2*j + 2*delta_o*w + 2*a_o*v + 2*f^2*VY)
    VE_o <- target - VYo_terms
    stopifnot(VE_o > 0)
    thetaNT_o <- with(parent, a_o*itlo + a_o*itol + 2*delta_o*gt + w)
    thetaT_o  <- 2*delta_o*k + thetaNT_o
    Yo_Yp_o   <- with(parent, delta_o*Omega + a_o*Gamma + f*VY + (delta_o*Omega + a_o*Gamma)*mu*VY)
    list(delta_o = delta_o, a_o = a_o, VE_o = VE_o, VY_o = target,
         VYo_terms = VYo_terms, thetaNT_o = thetaNT_o, thetaT_o = thetaT_o, Yo_Yp_o = Yo_Yp_o)
}
