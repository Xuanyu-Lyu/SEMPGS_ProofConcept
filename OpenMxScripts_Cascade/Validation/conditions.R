## Shared conditions for validating the SEM-PGS Cascade scripts (sourced by 00_truth.R and 02_fit_models.R).
##
## Generating values: one base parameter set, three mating mechanisms, two regimes. Within every
## regime x mechanism the base (delta, a, VE) are rescaled so the parental phenotypic variance is exactly
## 1 (the binary iterative-math wrappers), which lets the continuous and the binary (liability-threshold)
## models be fit to the same simulated populations. The offspring trait of the DiffTrait designs loads
## delta_o / a_o on the same PGS / LGS, receives the same F, and gets VE_o solved so its variance is 1
## too (the DiffTrait_Latent models need VY_o == VY_p).

cascadeValidationSetup <- function(cascadeDir){
    source(file.path(cascadeDir, "EquilibriumBinary", "00-UniIterativeMath_Cascade_Binary.R"), local = FALSE)
    source(file.path(cascadeDir, "DisequilibriumBinary", "00-UniIterativeMath_Cascade_DisEq_Binary.R"), local = FALSE)
    invisible(TRUE)
}

BASE <- list(delta = sqrt(.2), a = sqrt(.3), VE = .5, f = .15, am = .4, k = .5, j = .5,
             delta_o = .3, a_o = .45)
## The original OpenMxScripts test values (delta = .2, a = sqrt(.6), VE = .36). Used only for the
## exact-covariance coding check when BASE implies a non-positive-definite covariance matrix.
BASE_ORIGINAL_TEST <- modifyList(BASE, list(delta = .2, a = sqrt(.6), VE = .36))
MECHANISMS <- list(primary = c(AM_G = 1, AM_E = 1),     # gamma~ = Y
                   genetic = c(AM_G = 1, AM_E = 0),     # gamma~ = A (PGS + LGS parts)
                   social  = c(AM_G = 0, AM_E = 1))     # gamma~ = F + E
REGIMES <- c("Eq", "DisEq")
PREV <- c(parent = .25, offspring_diff = .40)           # binary prevalences (thresholds on the VY = 1 scale)
THRESH <- c(parent = qnorm(1 - PREV[["parent"]]), offspring_diff = qnorm(1 - PREV[["offspring_diff"]]))

## Iterative-math truth for one regime x mechanism: parent-trait fixed point `q` (VY == 1) and offspring trait `o`.
conditionTruth <- function(regime, mech, base = BASE){
    m <- MECHANISMS[[mech]]
    args <- list(delta = base$delta, a = base$a, f = base$f, VE = base$VE, am = base$am,
                 AM_G = m[["AM_G"]], AM_E = m[["AM_E"]], k = base$k, j = base$j)
    if (regime == "Eq"){
        q <- do.call(uniIterativeMath_Cascade_Binary, args)
        checkCascadeEquilibriumConsistency(q)
        o <- solve_VYo_Cascade_Equilibrium_Binary(base$delta_o, base$a_o, q)
    } else {
        q <- do.call(uniIterativeMath_Cascade_DisEq_Binary, args)
        checkCascadeDisEqConsistency(q)
        o <- solve_VYo_Cascade_DisEq_Binary(base$delta_o, base$a_o, q)
    }
    list(regime = regime, mech = mech, q = q, o = o)
}

## RDR heritability input, computed with the same formula as the scripts' RDR constraint
h2RDR <- function(a, delta, VE, VY, k = .5, j = .5) (2*a^2*j + 2*delta^2*k) * (2*a^2*j + 2*delta^2*k + VE) / VY

## How the mating multipliers are passed to each fit. gamma~ has no scale of its own, so one multiplier
## is fixed: AM_E = 1 (the MVN script's convention) unless the truth has AM_E = 0 (genetic homogamy), in
## which case AM_G = 1 is fixed and AM_E estimated. With latent parental phenotypes the mechanism is not
## identified, so both multipliers are fixed at the hypothesized (here: true) mechanism.
multiplierArgs <- function(mech, latent){
    m <- MECHANISMS[[mech]]
    if (latent) return(list(AM_G_value = m[["AM_G"]], AM_G_free = FALSE, AM_E_value = m[["AM_E"]], AM_E_free = FALSE))
    if (m[["AM_E"]] == 0) return(list(AM_G_value = 1, AM_G_free = FALSE, AM_E_value = .5, AM_E_free = TRUE))
    list(AM_G_value = .5, AM_G_free = TRUE, AM_E_value = 1, AM_E_free = FALSE)
}

DESIGNS <- c("SameTrait_ObservedYPYM", "SameTrait_LatentYPYM", "DiffTrait_LatentYPYM",
             "DiffTrait_ObservedYPYM_EstimatedAParent", "DiffTrait_ObservedYPYM_FixedAParent")
## which simulated dataset each design is fit to
DESIGN_DATA <- c(SameTrait_ObservedYPYM = "SameTrait_Observed", SameTrait_LatentYPYM = "SameTrait_Latent",
                 DiffTrait_LatentYPYM = "DiffTrait_Latent",
                 DiffTrait_ObservedYPYM_EstimatedAParent = "DiffTrait_Observed",
                 DiffTrait_ObservedYPYM_FixedAParent = "DiffTrait_Observed")

## Fitting function name for a regime x design x (binary?)
fitFunctionName <- function(regime, design, binary){
    fn <- paste0("fitUniSEMPGS_Cascade_", if (regime == "DisEq") "DisEq_" else "", if (binary) "Binary_" else "", design)
    sub("SameTrait_LatentYPYM$", "SameTrait_LatentParents", fn)
}
scriptPath <- function(cascadeDir, regime, design, binary){
    folder <- paste0(if (regime == "Eq") "Equilibrium" else "Disequilibrium", if (binary) "Binary" else "")
    file.path(cascadeDir, folder, paste0("UniSEMPGS_Cascade_", if (regime == "DisEq") "DisEq_" else "",
                                         if (binary) "Binary_" else "", design, ".R"))
}

## Truth (by OpenMx label) and extra fit arguments for one design, plus the model-implied covariance
## matrix (exactly the entries of that script's expCov) for the exact-covariance coding check.
designSpec <- function(tr, design, binary){
    q <- tr$q; o <- tr$o
    latent <- grepl("Latent", design)
    diff   <- grepl("DiffTrait", design)
    mult   <- multiplierArgs(tr$mech, latent)
    common <- c(AMGenMulti = q$AM_G, AMEnvMulti = q$AM_E, f11 = q$f, mu11 = q$mu)
    if (tr$regime == "Eq") common <- c(common, VF11 = q$VF, w11 = q$w, v11 = q$v, ic11 = q$ic, gc11 = q$gc, hc11 = q$hc)
    # variance of Yo in the model: the offspring trait's (DiffTrait), the offspring's own (DisEq), or VY (Eq)
    VY_o_model <- if (diff) o$VY_o else if (tr$regime == "DisEq") q$VY_off else q$VY
    if (!diff){
        truth <- c(common, VY11 = q$VY, VE11 = q$VE, delta11 = q$delta, a11 = q$a, Omega11 = q$Omega, Gamma11 = q$Gamma)
        if (tr$regime == "Eq") truth <- c(truth, gt11 = q$gt, ht11 = q$ht)
        thT <- q$thetaT; thNT <- q$thetaNT; YoYp <- q$Yo_Yp
    } else {
        truth <- c(common, VE_parent = q$VE, VE_offspring = o$VE_o, VE_p11 = q$VE, VE_o11 = o$VE_o,
                   VY11 = q$VY, VY_parent = q$VY, VY_offspring = o$VY_o, VY_p11 = q$VY, VY_o11 = o$VY_o,
                   delta_p11 = q$delta, a_p11 = q$a, delta_o11 = o$delta_o, a_o11 = o$a_o,
                   Omega_p11 = q$Omega, Gamma_p11 = q$Gamma)
        thT <- o$thetaT_o; thNT <- o$thetaNT_o; YoYp <- o$Yo_Yp_o
    }
    if (binary) truth <- c(truth, thresh_Yp11 = THRESH[["parent"]], thresh_Ym11 = THRESH[["parent"]],
                           thresh_Yo11 = if (diff) THRESH[["offspring_diff"]] else THRESH[["parent"]])

    args <- mult
    if (design == "SameTrait_LatentYPYM") args$h2_RDR <- h2RDR(q$a, q$delta, q$VE, q$VY)
    if (design == "DiffTrait_LatentYPYM"){
        args$h2_RDR_parent    <- h2RDR(q$a, q$delta, q$VE, q$VY)
        args$h2_RDR_offspring <- h2RDR(o$a_o, o$delta_o, o$VE_o, q$VY)   # VY_o == VY in this design
        args$mu_fixed         <- q$mu
    }
    if (design == "DiffTrait_ObservedYPYM_EstimatedAParent")
        args$h2_RDR_offspring <- h2RDR(o$a_o, o$delta_o, o$VE_o, o$VY_o)
    if (design == "DiffTrait_ObservedYPYM_FixedAParent"){
        args$h2_RDR_parent    <- h2RDR(q$a, q$delta, q$VE, q$VY)
        args$h2_RDR_offspring <- h2RDR(o$a_o, o$delta_o, o$VE_o, o$VY_o)
    }

    gc <- if (tr$regime == "Eq") q$gc else 0
    if (latent){
        cols <- c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")
        CM <- rbind(c(VY_o_model, thT, thNT, thT, thNT),
                    c(thT,  q$k + gc, gc, q$gt, q$gt),
                    c(thNT, gc, q$k + gc, q$gt, q$gt),
                    c(thT,  q$gt, q$gt, q$k + gc, gc),
                    c(thNT, q$gt, q$gt, gc, q$k + gc))
    } else {
        cols <- c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1")
        CM <- rbind(c(q$VY, q$Yp_Ym, YoYp, q$Omega, q$Omega, q$Yp_PGSm, q$Yp_PGSm),
                    c(q$Yp_Ym, q$VY, YoYp, q$Yp_PGSm, q$Yp_PGSm, q$Omega, q$Omega),
                    c(YoYp, YoYp, VY_o_model, thT, thNT, thT, thNT),
                    c(q$Omega, q$Yp_PGSm, thT, q$k + gc, gc, q$gt, q$gt),
                    c(q$Omega, q$Yp_PGSm, thNT, gc, q$k + gc, q$gt, q$gt),
                    c(q$Yp_PGSm, q$Omega, thT, q$gt, q$gt, q$k + gc, gc),
                    c(q$Yp_PGSm, q$Omega, thNT, q$gt, q$gt, gc, q$k + gc))
    }
    dimnames(CM) <- list(cols, cols)
    thresholds <- if (latent) c(Yo1 = if (diff) THRESH[["offspring_diff"]] else THRESH[["parent"]])
                  else c(Yp1 = THRESH[["parent"]], Ym1 = THRESH[["parent"]],
                         Yo1 = if (diff) THRESH[["offspring_diff"]] else THRESH[["parent"]])
    list(truth = truth, args = args, CM = CM, thresholds = thresholds, dataset = DESIGN_DATA[[design]])
}
