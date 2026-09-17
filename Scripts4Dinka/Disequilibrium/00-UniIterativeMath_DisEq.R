## Univariate (scalar) iterative math for the SEM-PGS model in DISEQUILIBRIUM.
## Matches the disequilibrium model equations: vertical transmission is iterated to its
## random-mating equilibrium (mu = 0, so gc = hc = ic = 0 and w = 2*f*Omega with no mu terms),
## then assortative mating happens for the FIRST time in the observed parents' generation.
## This exactly satisfies the DisEq scripts' self-referential constraints
## (Omega = delta*k + .5*w with w = 2*f*Omega; VF = 2*f^2*VY).

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
