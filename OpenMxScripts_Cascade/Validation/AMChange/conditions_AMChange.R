## Shared conditions for validating the Cascade AMChange scripts (two-regime disequilibrium: a change in assortment
## strength), sourced by 00_truth.R and 02_fit_models.R after they set `baseDir` (this folder) and `cascadeDir`.
##
## Generating values: BASE and the three mating mechanisms of ../conditions.R. The grandparents (_pp) mate with
## spousal gamma~ correlation AM_PP = .4 (the Eq validation's am) long enough to reach equilibrium; the parents (_p)
## then mate once with spousal correlation am_p, one of BRANCHES. Within every mechanism (delta, a, VE) are rescaled
## so the parents' phenotypic variance is exactly 1, which lets the continuous and the binary models be fit to the
## same simulated populations (as in ../conditions.R). The rescaling depends only on the grandparental equilibrium,
## so all branches of a mechanism share delta, a, VE; only the DiffTrait offspring trait's VE_o (solved so its
## variance is 1) differs between branches.

source(file.path(dirname(baseDir), "conditions.R"))     # BASE, MECHANISMS, PREV, THRESH, h2RDR, multiplierArgs, DESIGNS

AM_PP       <- .4
BRANCHES    <- c(weaker = .2, stopped = 0, stronger = .6)    # am_p: the parents' own spousal gamma~ correlation
MAIN_BRANCH <- "weaker"
BASE_MODEL1 <- modifyList(BASE, list(delta = sqrt(.5), a = 0))   # BASE's .5 genetic variance, all in the PGS

AMC_DESIGNS     <- c(DESIGNS, "Model1")
AMC_DESIGN_DATA <- c(DESIGN_DATA, Model1 = "SameTrait_Latent")

amchangeSetup <- function(cascadeDir){
    source(file.path(cascadeDir, "DisequilibriumAMChangeBinary", "00-UniIterativeMath_Cascade_AMChange_Binary.R"),
           local = FALSE)
    invisible(TRUE)
}

## Iterative-math truth for one mechanism x branch: the parent-trait point `q` (VY_p == 1) and the offspring trait `o`
amchangeTruth <- function(mech, branch, base = BASE){
    m <- MECHANISMS[[mech]]
    q <- uniIterativeMath_Cascade_AMChange_Binary(delta = base$delta, a = base$a, f = base$f, VE = base$VE,
                                                  am_pp = AM_PP, am_p = BRANCHES[[branch]],
                                                  AM_G = m[["AM_G"]], AM_E = m[["AM_E"]], k = base$k, j = base$j)
    checkCascadeAMChangeConsistency(q)
    o <- solve_VYo_Cascade_AMChange_Binary(base$delta_o, base$a_o, q)
    list(mech = mech, branch = branch, q = q, o = o)
}

amchangeFitName <- function(design, binary){
    fn <- paste0("fitUniSEMPGS_Cascade_AMChange_", if (binary) "Binary_" else "", design)
    sub("SameTrait_LatentYPYM$", "SameTrait_LatentParents", fn)
}
amchangeScriptPath <- function(cascadeDir, design, binary){
    file.path(cascadeDir, paste0("DisequilibriumAMChange", if (binary) "Binary" else ""),
              paste0("UniSEMPGS_Cascade_AMChange_", if (binary) "Binary_" else "", design, ".R"))
}

## Truth (by OpenMx label), extra fit arguments, and the model-implied covariance matrix (exactly the entries of the
## script's expCov) for one design.
designSpecAMChange <- function(tr, design, binary){
    q <- tr$q; o <- tr$o
    latent <- grepl("Latent|Model1", design)
    diff   <- grepl("DiffTrait", design)
    model1 <- design == "Model1"
    # With mu_p = 0 (AM stopped) the parents' mating carries no information on what spouses assort on, so the
    # mechanism is fixed at the hypothesized (true) one even with observed parents: in the DiffTrait_Observed designs a
    # free AM_G then trades off against mu_pp at no cost in fit.
    mult   <- multiplierArgs(tr$mech, latent || q$mu_p == 0)

    truth <- c(AMGenMulti = q$AM_G, AMEnvMulti = q$AM_E, f11 = q$f, mu_pp11 = q$mu_pp, mu_p11 = q$mu_p,
               VY_p11 = q$VY_p, VF_p11 = q$VF_p, w_p11 = q$w_p, gc_p11 = q$gc_p, Omega_p11 = q$Omega_p)
    if (!model1) truth <- c(truth, v_p11 = q$v_p, hc_p11 = q$hc_p, ic_p11 = q$ic_p, Gamma_p11 = q$Gamma_p)
    if (!diff){
        truth <- c(truth, VE11 = q$VE, delta11 = q$delta, a11 = q$a)
        VYo <- q$VY_o; thT <- q$thetaT; thNT <- q$thetaNT; YoYp <- q$Yo_Yp
    } else {
        truth <- c(truth, VE_p11 = q$VE, VE_o11 = o$VE_o, VY_o11 = o$VY_o, delta_p11 = q$delta, a_p11 = q$a,
                   delta_o11 = o$delta_o, a_o11 = o$a_o)
        VYo <- o$VY_o; thT <- o$thetaT_o; thNT <- o$thetaNT_o; YoYp <- o$Yo_Yp_o
    }
    if (model1) truth <- truth[names(truth) != "a11"]
    if (binary){
        truth <- truth[!names(truth) %in% c("VY_p11", "VY_o11")]       # fixed at 1 on the liability scale
        truth <- c(truth, thresh_Yp11 = THRESH[["parent"]], thresh_Ym11 = THRESH[["parent"]],
                   thresh_Yo11 = if (diff) THRESH[["offspring_diff"]] else THRESH[["parent"]])
    }

    args <- mult
    if (design %in% c("SameTrait_LatentYPYM")) args$h2_RDR <- h2RDR(q$a, q$delta, q$VE, q$VY_p)
    if (design == "DiffTrait_LatentYPYM"){
        args$h2_RDR_parent    <- h2RDR(q$a, q$delta, q$VE, q$VY_p)
        args$h2_RDR_offspring <- h2RDR(o$a_o, o$delta_o, o$VE_o, q$VY_p)      # VY_o == VY_p in this design
        # one of the two copaths must be supplied (see the script header). The parents' own mu_p is the natural
        # one, but with mu_p = 0 (AM stopped) gt_p = 0 carries no information on Omega~, so mu_pp is fixed instead.
        if (q$mu_p != 0) args$mu_p_fixed <- q$mu_p else args$mu_pp_fixed <- q$mu_pp
    }
    if (design == "DiffTrait_ObservedYPYM_EstimatedAParent")
        args$h2_RDR_offspring <- h2RDR(o$a_o, o$delta_o, o$VE_o, o$VY_o)
    if (design == "DiffTrait_ObservedYPYM_FixedAParent"){
        args$h2_RDR_parent    <- h2RDR(q$a, q$delta, q$VE, q$VY_p)
        args$h2_RDR_offspring <- h2RDR(o$a_o, o$delta_o, o$VE_o, o$VY_o)
    }

    gc <- q$gc_p; gt <- q$gt_p; k <- q$k
    if (latent){
        cols <- c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")
        CM <- rbind(c(VYo,  thT,    thNT,   thT,    thNT),
                    c(thT,  k + gc, gc,     gt,     gt),
                    c(thNT, gc,     k + gc, gt,     gt),
                    c(thT,  gt,     gt,     k + gc, gc),
                    c(thNT, gt,     gt,     gc,     k + gc))
    } else {
        cols <- c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1")
        CM <- rbind(c(q$VY_p,    q$Yp_Ym,   YoYp, q$Omega_p, q$Omega_p, q$Yp_PGSm, q$Yp_PGSm),
                    c(q$Yp_Ym,   q$VY_p,    YoYp, q$Yp_PGSm, q$Yp_PGSm, q$Omega_p, q$Omega_p),
                    c(YoYp,      YoYp,      VYo,  thT,       thNT,      thT,       thNT),
                    c(q$Omega_p, q$Yp_PGSm, thT,  k + gc,    gc,        gt,        gt),
                    c(q$Omega_p, q$Yp_PGSm, thNT, gc,        k + gc,    gt,        gt),
                    c(q$Yp_PGSm, q$Omega_p, thT,  gt,        gt,        k + gc,    gc),
                    c(q$Yp_PGSm, q$Omega_p, thNT, gt,        gt,        gc,        k + gc))
    }
    dimnames(CM) <- list(cols, cols)
    thresholds <- if (latent) c(Yo1 = if (diff) THRESH[["offspring_diff"]] else THRESH[["parent"]])
                  else c(Yp1 = THRESH[["parent"]], Ym1 = THRESH[["parent"]],
                         Yo1 = if (diff) THRESH[["offspring_diff"]] else THRESH[["parent"]])
    list(truth = truth, args = args, CM = CM, thresholds = thresholds, dataset = AMC_DESIGN_DATA[[design]])
}
