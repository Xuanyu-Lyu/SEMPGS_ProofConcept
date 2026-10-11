## Step 0: iterative-math truth for the AMChange validation, plus algebra-level checks.
##   A. Reductions: the AMChange math with am_pp = 0 must equal the Disequilibrium/ math, and with am_pp = am_p the
##      Equilibrium/ math, for every quantity the fitting scripts use (all three mechanisms and a mixed one).
##   B. Truth for every mechanism x branch (parents' VY_p == 1), written to output/truth/, and the generating values
##      01_simulate.py reads (output/truth/generating_params.csv).
##   C. Every design's model-implied covariance matrix is positive definite.
##   D. The offspring quantities next to the PDF Part VII Section 38 forms they correct.
## Run with: Rscript 00_truth.R

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
cascadeDir <- dirname(dirname(baseDir))
suppressPackageStartupMessages(library(data.table))
source(file.path(baseDir, "conditions_AMChange.R"))
cascadeValidationSetup(cascadeDir)      # the Eq and DisEq math, for the reductions
amchangeSetup(cascadeDir)
outDir <- file.path(baseDir, "output", "truth"); dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

## ---------------- A. Reductions ----------------
cat("A. Reductions of the AMChange math (max absolute difference over all compared quantities)\n")
red <- rbindlist(lapply(list(primary = c(1, 1), genetic = c(1, 0), social = c(0, 1), mixed = c(.6, 1.3)), function(m){
    args <- list(delta = BASE$delta, a = BASE$a, f = BASE$f, VE = BASE$VE, AM_G = m[1], AM_E = m[2])
    out <- list()
    for (am in c(.2, .4)){
        q <- do.call(uniIterativeMath_Cascade_AMChange, c(args, am_pp = 0, am_p = am)); checkCascadeAMChangeConsistency(q)
        d <- do.call(uniIterativeMath_Cascade_DisEq, c(args, am = am))
        dd <- c(q$VY_p - d$VY, q$Omega_p - d$Omega, q$Gamma_p - d$Gamma, q$w_p - d$w, q$VF_p - d$VF, q$gc_p,
                q$Omega_td_p - d$Omega_td, q$tau_p - d$tau, q$mu_p - d$mu, q$gt_p - d$gt, q$it_p - d$itol,
                q$w_o - d$w_o, q$v_o - d$v_o, q$VF_o - d$VF_o, q$thetaNT - d$thetaNT, q$thetaT - d$thetaT,
                q$VY_o - d$VY_off, q$Yo_Yp - d$Yo_Yp, q$Yp_Ym - d$Yp_Ym, q$Yp_PGSm - d$Yp_PGSm)
        q2 <- do.call(uniIterativeMath_Cascade_AMChange, c(args, am_pp = am, am_p = am)); checkCascadeAMChangeConsistency(q2)
        e  <- do.call(uniIterativeMath_Cascade, c(args, am = am))
        de <- c(q2$VY_p - e$VY, q2$VY_o - e$VY, q2$Omega_p - e$Omega, q2$gc_p - e$gc, q2$gt_p - e$gt, q2$ic_p - e$ic,
                q2$it_p - e$ic, q2$w_p - e$w, q2$w_o - e$w, q2$VF_o - e$VF, q2$thetaNT - e$thetaNT,
                q2$thetaT - e$thetaT, q2$Yo_Yp - e$Yo_Yp, q2$Yp_Ym - e$Yp_Ym, q2$mu_pp - e$mu, q2$mu_p - e$mu)
        out[[length(out) + 1]] <- data.table(AM_G = m[1], AM_E = m[2], am = am,
                                             vs_DisEq_at_am_pp_0 = max(abs(dd)), vs_Eq_at_am_pp_eq_am_p = max(abs(de)))
    }
    rbindlist(out)
}), idcol = "mechanism")
print(red, digits = 2)
stopifnot(max(red$vs_DisEq_at_am_pp_0) < 1e-10, max(red$vs_Eq_at_am_pp_eq_am_p) < 1e-10)

## ---------------- B. Truth ----------------
grid <- expand.grid(mech = names(MECHANISMS), branch = names(BRANCHES), stringsAsFactors = FALSE)
TR  <- lapply(seq_len(nrow(grid)), function(i) amchangeTruth(grid$mech[i], grid$branch[i]))
TR1 <- lapply(seq_len(nrow(grid)), function(i) amchangeTruth(grid$mech[i], grid$branch[i], base = BASE_MODEL1))
scalarPart <- function(x) unlist(x[vapply(x, function(v) is.numeric(v) && length(v) == 1, TRUE)])
truthLong <- rbindlist(lapply(c(TR, TR1), function(tr){
    v <- c(scalarPart(tr$q), setNames(scalarPart(tr$o), paste0("diff_", names(scalarPart(tr$o)))),
           setNames(unlist(tr$q$pdf_partVII), paste0("pdf38_", names(tr$q$pdf_partVII))))
    data.table(base = if (tr$q$a == 0) "BASE_MODEL1" else "BASE", mechanism = tr$mech, branch = tr$branch,
               quantity = names(v), value = unname(v))
}))
fwrite(truthLong, file.path(outDir, "truth_long.csv"))
gp <- rbindlist(lapply(TR, function(tr) with(tr$q, data.table(
    mechanism = tr$mech, branch = tr$branch, delta = delta, a = a, VE = VE, f = f, am_pp = am_pp, am_p = am_p,
    AM_G = AM_G, AM_E = AM_E, k = k, j = j, delta_o = tr$o$delta_o, a_o = tr$o$a_o, VE_o = tr$o$VE_o,
    scale = scale, iter_generations = generations))))
fwrite(gp, file.path(outDir, "generating_params.csv"))
cat("\nB. Generating values (parents' VY_p = 1), written to output/truth/generating_params.csv\n")
print(gp, digits = 4)

## ---------------- C. Positive definiteness ----------------
pd <- rbindlist(lapply(seq_len(nrow(grid)), function(i) rbindlist(lapply(AMC_DESIGNS, function(ds){
    tr <- if (ds == "Model1") TR1[[i]] else TR[[i]]
    CM <- designSpecAMChange(tr, ds, FALSE)$CM
    data.table(mechanism = grid$mech[i], branch = grid$branch[i], design = ds,
               min_eigenvalue = min(eigen(CM, only.values = TRUE)$values), symmetric = max(abs(CM - t(CM))) < 1e-12)
}))))
cat("\nC. Smallest eigenvalue of each design's implied covariance matrix:", sprintf("%.4f", min(pd$min_eigenvalue)),
    if (all(pd$min_eigenvalue > 0) && all(pd$symmetric)) "(all positive definite)\n" else "(NOT all PD)\n")
stopifnot(all(pd$min_eigenvalue > 0), all(pd$symmetric))

## ---------------- D. Corrected offspring quantities vs PDF Part VII Section 38 ----------------
cat("\nD. Offspring quantities (per data column) vs PDF Part VII Section 38 as printed\n")
cmp <- rbindlist(lapply(TR, function(tr) with(tr$q, data.table(
    mechanism = tr$mech, branch = tr$branch, mu_pp = mu_pp, mu_p = mu_p, gc_p = gc_p, gt_p = gt_p,
    thetaNT = thetaNT, thetaNT_pdf38 = pdf_partVII$thetaNT2 / 2,
    thetaT = thetaT, thetaT_pdf38 = pdf_partVII$thetaT2 / 2,
    VY_o = VY_o, VY_o_pdf38 = pdf_partVII$VY2))))
print(cmp, digits = 4)
fwrite(cmp, file.path(outDir, "offspring_vs_pdf38.csv"))
