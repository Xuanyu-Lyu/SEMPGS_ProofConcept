## Step 0: iterative-math truth for every regime x mating mechanism.
## Writes (under output/truth/):
##   generating_params.csv    the (VY-standardized) generating values that 01_simulate.py feeds to GeneEvolve
##   iterative_math.csv       every scalar the Cascade iterative math returns, per regime x mechanism
##   iterative_trajectory.csv the iterative math's generation-by-generation trajectory
## Run with: Rscript 00_truth.R   (from anywhere)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
cascadeDir <- dirname(baseDir)
source(file.path(baseDir, "conditions.R"))
cascadeValidationSetup(cascadeDir)

outDir <- file.path(baseDir, "output", "truth"); dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

gen_rows <- list(); val_rows <- list(); traj_rows <- list()
for (regime in REGIMES) for (mech in names(MECHANISMS)){
    tr <- conditionTruth(regime, mech); q <- tr$q; o <- tr$o
    gen_rows[[length(gen_rows) + 1]] <- data.frame(
        regime = regime, mechanism = mech, delta = q$delta, a = q$a, VE = q$VE, f = q$f, am = q$am,
        AM_G = q$AM_G, AM_E = q$AM_E, k = q$k, j = q$j,
        delta_o = o$delta_o, a_o = o$a_o, VE_o = o$VE_o, scale = q$scale,
        iter_generations = q$generations)
    scal <- q[sapply(q, function(x) is.numeric(x) && length(x) == 1)]
    if (regime == "Eq") scal$VY_off <- q$VY    # at equilibrium the offspring's VY is the parents'
    vals <- c(unlist(scal), hapvar_obs = q$k + q$gc, hapvar_lat = q$j + q$hc, r_gamma = q$am,
              it = q$itol, VY_o = o$VY_o, VE_o = o$VE_o, Yo_Yp_o = o$Yo_Yp_o,
              thetaT_o = o$thetaT_o, thetaNT_o = o$thetaNT_o)
    if (regime == "DisEq") vals <- c(vals, setNames(unlist(q$pdf_partIV), paste0("PDF_", names(q$pdf_partIV))))
    val_rows[[length(val_rows) + 1]] <- data.frame(regime = regime, mechanism = mech,
                                                   quantity = names(vals), value = unname(vals))
    traj_rows[[length(traj_rows) + 1]] <- cbind(regime = regime, mechanism = mech, q$trajectory)
    cat(sprintf("%-5s %-8s VY=%.4f mu=%.4f Vgamma=%.4f tau=%.4f gt=%.5f VE_o=%.4f (s=%.4f, %d iterations)\n",
                regime, mech, q$VY, q$mu, q$Vgamma, q$tau, q$gt, o$VE_o, q$scale, q$generations))
}
write.csv(do.call(rbind, gen_rows), file.path(outDir, "generating_params.csv"), row.names = FALSE)
write.csv(do.call(rbind, val_rows), file.path(outDir, "iterative_math.csv"), row.names = FALSE)
traj <- do.call(rbind, lapply(traj_rows, function(d){ d[setdiff(Reduce(union, lapply(traj_rows, names)), names(d))] <- NA; d }))
write.csv(traj, file.path(outDir, "iterative_trajectory.csv"), row.names = FALSE)
cat("wrote", outDir, "\n")
