## Univariate (scalar) iterative math for the SEM-PGS model at equilibrium, standardized
## to a liability scale (VY == 1) for use with the binary/liability-threshold fitting scripts.
##
## The base recursion (uniIterativeMath / checkEquilibriumConsistency) is copied verbatim from
## Equilibrium/00-UniIterativeMath.R. On top of it, uniIterativeMath_Binary() exploits the fact
## that the equilibrium fixed point is homogeneous under the rescaling
##     (delta, a, VE) -> (s*delta, s*a, s^2*VE),   s = 1/sqrt(VY)
## with f and am held fixed: re-running the SAME recursion with these rescaled inputs converges
## to a fixed point with VY == 1 exactly, while Omega/Gamma/w/v/thetaT/thetaNT/Yp_PGSm scale by s,
## Yo_Yp/Yp_Ym scale by s^2 (Yp_Ym becomes exactly `am`), mu becomes exactly `am`, and the
## PGS-PGS covariances gt/ht/gc/hc/ic/itlo/itol are unchanged (they never depended on the Y scale).
## This was verified numerically to match a direct algebraic rescaling of every derived quantity.

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

## ---------------------------------------------------------------------------
## Binary/liability-standardized wrapper: same-trait / parent-trait equilibrium
## point re-expressed on a scale where VY == 1 exactly.
## ---------------------------------------------------------------------------
uniIterativeMath_Binary <- function(delta, a, f, VE, am, k = .5, j = .5, gens = 100, tol = 1e-12){
    raw <- uniIterativeMath(delta = delta, a = a, f = f, VE = VE, am = am, k = k, j = j, gens = gens, tol = tol)
    stopifnot(raw$converged)
    s <- 1 / sqrt(raw$VY)
    std <- uniIterativeMath(delta = s*delta, a = s*a, f = f, VE = s^2*VE, am = am, k = k, j = j, gens = gens, tol = tol)
    stopifnot(std$converged, abs(std$VY - 1) < 1e-8)
    std$scale <- s
    std
}

## Given an already-standardized parent-trait equilibrium point (VY_p == 1) and a choice of
## offspring generation's delta_o/a_o, solve for the VE_o that makes the offspring's implied
## liability variance equal exactly 1 too (equilibrium DiffTrait formula).
solve_VYo_Equilibrium_Binary <- function(delta_o, a_o, k, j, parent, target = 1){
    VYo_terms <- with(parent,
        2*delta_o^2*k + 4*delta_o^2*gc + 2*a_o^2*j + 4*a_o^2*hc +
        8*delta_o*ic*a_o + 2*delta_o*w + 2*a_o*v + 2*f^2*VY*(1 + VY*mu))
    VE_o <- target - VYo_terms
    stopifnot(VE_o > 0)
    thetaNT_o <- with(parent, 2*delta_o*gc + 2*a_o*ic + .5*w)
    thetaT_o  <- delta_o*k + thetaNT_o
    Yo_Yp_o   <- with(parent, delta_o*Omega + a_o*Gamma + (delta_o*Omega + a_o*Gamma)*mu*VY + f*VY + f*VY^2*mu)
    list(delta_o = delta_o, a_o = a_o, VE_o = VE_o, VY_o = target,
         VYo_terms = VYo_terms, thetaNT_o = thetaNT_o, thetaT_o = thetaT_o, Yo_Yp_o = Yo_Yp_o)
}
