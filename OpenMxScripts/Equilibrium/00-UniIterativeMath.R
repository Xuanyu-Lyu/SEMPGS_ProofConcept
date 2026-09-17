## Univariate (scalar) iterative math for the SEM-PGS model at equilibrium.
## Scalar translation of the bivariate recursion in PaperScripts/09-fitWithMVN.R
## (fix-the-mate-correlation convention; scalar mu is symmetric so all transposes drop).
## Iterates the AM + VT recursion from a random-mating base until convergence and
## returns every equilibrium quantity plus the derived covariance pieces used in expCov.

uniIterativeMath <- function(delta, a, f, VE, am, k = .5, j = .5, gens = 100, tol = 1e-12){

    # t0: random-mating base population, no VT history
    VY <- 2*delta^2*k + 2*a^2*j + VE
    mu <- am / VY                          # scalar mu = solve(VY) %*% (am*VY) %*% solve(VY)
    gc <- hc <- ic <- 0
    gt <- ht <- itlo <- itol <- 0
    w  <- v  <- 0
    VF <- 2*f^2*VY + 2*f^2*mu*VY^2
    Omega <- delta*k
    Gamma <- a*j

    converged <- FALSE
    it_used <- 1
    for (it in 2:gens){
        prev <- c(Omega, Gamma, VY, mu, gt, ht, w, v, VF, ic)

        # same update order as 09-fitWithMVN.R:123-177
        Omega <- 2*delta*gc + delta*k + .5*w + 2*a*ic
        Gamma <- 2*a*hc + 2*delta*ic + a*j + .5*v
        VY    <- 2*delta*Omega + 2*a*Gamma + w*delta + v*a + VF + VE
        mu    <- am / VY                   # mate correlation held fixed across generations
        gt    <- Omega^2 * mu;  gc <- gt
        ht    <- Gamma^2 * mu;  hc <- ht
        w     <- 2*f*Omega + 2*f*VY*mu*Omega
        v     <- 2*f*Gamma + 2*f*VY*mu*Gamma
        VF    <- 2*f^2*VY + 2*f^2*mu*VY^2
        itlo  <- Gamma*mu*Omega
        itol  <- Omega*mu*Gamma
        ic    <- .5*(itlo + itol)

        it_used <- it
        if (max(abs(c(Omega, Gamma, VY, mu, gt, ht, w, v, VF, ic) - prev)) < tol){
            converged <- TRUE
            break
        }
    }

    # derived covariance pieces (same expressions as the equilibrium fitting scripts)
    thetaNT <- 2*delta*gc + 2*a*ic + .5*w
    thetaT  <- delta*k + thetaNT
    Yp_PGSm <- VY*mu*Omega
    Yp_Ym   <- mu*VY^2
    Yo_Yp   <- delta*Omega + a*Gamma + (delta*Omega + a*Gamma)*mu*VY + f*VY + f*VY^2*mu

    list(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j,
         VY = VY, mu = mu, Omega = Omega, Gamma = Gamma,
         gt = gt, ht = ht, gc = gc, hc = hc, itlo = itlo, itol = itol, ic = ic,
         w = w, v = v, VF = VF,
         thetaNT = thetaNT, thetaT = thetaT,
         Yp_PGSm = Yp_PGSm, Yp_Ym = Yp_Ym, Yo_Yp = Yo_Yp,
         converged = converged, generations = it_used)
}

## Assert that the converged quantities satisfy the equilibrium model constraints
## (the mxConstraint algebra of the equilibrium fitting scripts).
checkEquilibriumConsistency <- function(q, tol = 1e-8){
    with(q, {
        stopifnot(
            abs(Omega - (2*delta*gc + delta*k + .5*w + 2*a*ic)) < tol,
            abs(Gamma - (2*a*hc + 2*delta*ic + a*j + .5*v))     < tol,
            abs(VY - (2*delta*Omega + 2*a*Gamma + w*delta + v*a + VF + VE)) < tol,
            abs(VF - (2*f^2*VY + 2*f^2*mu*VY^2)) < tol,
            abs(w - (2*f*Omega + 2*f*VY*mu*Omega)) < tol,
            abs(v - (2*f*Gamma + 2*f*VY*mu*Gamma)) < tol,
            abs(gt - Omega^2*mu) < tol,
            abs(gc - gt) < tol,
            abs(ht - Gamma^2*mu) < tol,
            abs(ic - .5*(itlo + itol)) < tol,
            abs(mu - am/VY) < tol
        )
    })
    invisible(TRUE)
}
