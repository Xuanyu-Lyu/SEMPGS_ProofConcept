## Step 5: do RDR-anchored models reduce the G x Age bias of f? Four Cascade models that use RDR heritability are fit to
## the same G x Age trio data as Model 2 (03_fit_sim.R):
##   EstA     UniSEMPGS_Cascade_DiffTrait_ObservedYPYM_EstimatedAParent.R  observed parents; RDR fixes the offspring's a_o,
##                                                                        a_p is estimated
##   FixA     UniSEMPGS_Cascade_DiffTrait_ObservedYPYM_FixedAParent.R      observed parents; RDR fixes a_p and a_o
##   SameLat  UniSEMPGS_Cascade_SameTrait_LatentYPYM.R                     latent parents (Yo and the 4 PGS haplotypes only);
##                                                                        one delta and a; RDR of the child fixes a
##   DiffLat  UniSEMPGS_Cascade_DiffTrait_LatentYPYM.R                     latent parents; own delta and a for parents and
##                                                                        offspring, RDR fixes both a; mu fixed at its true
##                                                                        value (the script needs it); one VY shared by
##                                                                        the parents' and the child's trait
## Observed-parent models use the multipliers of Model 2 (AM_G or AM_E estimated); latent-parent models fix both at the
## true mechanism, as their scripts require.
## RDR inputs: the true heritability RDR would estimate, in the scripts' RDR form h2 = VG0 * (VG0 + VE) / VY, where VG0
## is the base-population genetic variance and VE the residual (h2RDR() in Validation/conditions.R):
##   parents    the adult trait: VG0 = 2 delta^2 k + 2 a^2 j, over the parents' VY.
##   offspring  the child phenotype as RDR in children would see it: the child's whole genetic variance, shared part
##              r^2 VG0 plus child-specific part (1 - r^2) VG0 = VG0, over the child's own VY.
## The models treat the child-specific part as residual (it is uncorrelated with the PGS and the latent score), so this
## child RDR is larger than the value that would make them correctly specified (r^2 VG0 in place of VG0).
## Each fit: the free model from the true parental values and again from f = -.2, 0, .1, .2, .4 (as testF); best kept.
## One CSV per fit in output/rdr_fits/ (resumable), combined into output/rdr_fits.csv.
## Run with: Rscript 05_fit_rdr.R [--workers 3] [--out output]   (--out: the output root holding data/; output_snpsets or
##           output_overlap for 06_simulate_overlap.py)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
getArg  <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else as.numeric(a[i + 1])
}
workers <- getArg("--workers", 3)
outArg  <- match("--out", commandArgs(trailingOnly = TRUE))
outRoot <- if (is.na(outArg)) "output" else commandArgs(trailingOnly = TRUE)[outArg + 1]   # output root (default output/)

# one thread per process, and one OpenMx thread per fit (the fitting scripts ask for parallel::detectCores())
Sys.setenv(OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(MASS); library(parallel) })
assignInNamespace("detectCores", function(...) 1L, ns = "parallel")
source(file.path(baseDir, "conditions.R"))

RDR_MODELS <- c(EstA    = "OpenMxScripts_Cascade/Equilibrium/UniSEMPGS_Cascade_DiffTrait_ObservedYPYM_EstimatedAParent.R",
                FixA    = "OpenMxScripts_Cascade/Equilibrium/UniSEMPGS_Cascade_DiffTrait_ObservedYPYM_FixedAParent.R",
                SameLat = "OpenMxScripts_Cascade/Equilibrium/UniSEMPGS_Cascade_SameTrait_LatentYPYM.R",
                DiffLat = "OpenMxScripts_Cascade/Equilibrium/UniSEMPGS_Cascade_DiffTrait_LatentYPYM.R")
DESIGN_OF  <- c(EstA = "DiffTrait_ObservedYPYM_EstimatedAParent", FixA = "DiffTrait_ObservedYPYM_FixedAParent",
                SameLat = "SameTrait_LatentYPYM", DiffLat = "DiffTrait_LatentYPYM")
LATENT     <- c(EstA = FALSE, FixA = FALSE, SameLat = TRUE, DiffLat = TRUE)
scriptPath <- c(scriptPath, RDR_MODELS)              # loadBuilder() (FTestValidity) reads scriptPath
RDR_BUILDERS <- lapply(setNames(names(RDR_MODELS), names(RDR_MODELS)), loadBuilder)

## truth by model (each design's labels; parents' values, f = .15), RDR inputs and mu for one mechanism x u
rdrSpec <- function(mech, u){
    s <- gxageSpec(mech, F_TRUE, u); q <- s$q
    list(truth = lapply(DESIGN_OF, function(d) designSpec(s$tr, d, binary = FALSE)$truth),
         h2_parent = h2RDR(q$a, q$delta, q$VE, q$VY),
         h2_child  = h2RDR(q$a, q$delta, q$VE, s$o$VY_o),
         mu = q$mu)
}
## arguments of each model's fitting function (besides data_path)
rdrArgs <- function(model, mech, sp){
    mult <- multiplierArgs(mech, latent = LATENT[[model]])
    switch(model,
           EstA    = c(mult, list(h2_RDR_offspring = sp$h2_child)),
           FixA    = c(mult, list(h2_RDR_parent = sp$h2_parent, h2_RDR_offspring = sp$h2_child)),
           SameLat = c(mult, list(h2_RDR = sp$h2_child)),
           DiffLat = c(mult, list(h2_RDR_parent = sp$h2_parent, h2_RDR_offspring = sp$h2_child, mu_fixed = sp$mu)))
}
USES_PARENT_RDR <- c(EstA = FALSE, FixA = TRUE, SameLat = FALSE, DiffLat = TRUE)

## the free-model search of testF (FTestValidity), without the f = 0 fit
fitFree <- function(builder, data_path, args, truth, extraTries = 5){
    m <- do.call(builder, c(list(data_path = data_path), args))
    lab <- intersect(names(truth), names(omxGetParameters(m)))
    lb  <- omxGetParameters(m, fetch = "lbound")[lab]
    m   <- omxSetParameters(m, labels = lab, values = ifelse(is.na(lb), truth[lab], pmax(truth[lab], lb + .01)))
    m0  <- mxOption(mxOption(m, "Calculate Hessian", "No"), "Standard Errors", "No")
    fit <- safeFit(m0, extraTries)
    for (fs in F_STARTS) fit <- better(fit, safeFit(omxSetParameters(m0, labels = "f11", values = fs), extraTries))
    if (!is.finite(m2LL(fit))) return(NULL)
    se <- tryCatch(suppressWarnings(mxRun(mxOption(mxOption(fit, "Calculate Hessian", "Yes"), "Standard Errors", "Yes"),
                                          silent = TRUE, suppressWarnings = TRUE)), error = function(e) NULL)
    if (!is.null(se) && m2LL(se) <= m2LL(fit) + 1e-6) fit <- se
    fit
}

dataDir <- file.path(baseDir, outRoot, "data")
fitDir  <- file.path(baseDir, outRoot, "rdr_fits"); dir.create(fitDir, recursive = TRUE, showWarnings = FALSE)
tmpDir  <- file.path(baseDir, outRoot, "tmp");      dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(dataDir, "^Eq_(primary|genetic|social)_(rep[0-9]{3}|pooled)_u[0-9]{3}\\.tsv$")
if (length(files) == 0) stop("No datasets in ", dataDir, " - run 01_simulate.py first")
tasks <- CJ(model = names(RDR_MODELS), file = files)
tasks[, out := file.path(fitDir, paste0(model, "_", sub("\\.tsv$", ".csv", file)))]
## one spec per mechanism x u level in the data (U_VALUES from 01_simulate.py; u = .33 or .5 from 06_simulate_overlap.py)
conds <- unique(data.table(mech = sub("^Eq_([a-z]+)_.*$", "\\1", files),
                           u = as.integer(sub("^.*_u([0-9]{3})\\.tsv$", "\\1", files)) / 100))
SPECS <- list()
for (i in seq_len(nrow(conds))) SPECS[[sprintf("%s_%03d", conds$mech[i], round(100 * conds$u[i]))]] <- rdrSpec(conds$mech[i], conds$u[i])

fitOne <- function(k){
    t <- tasks[k]
    if (file.exists(t$out)) return(invisible(t$out))
    m    <- regmatches(t$file, regexec("^Eq_([a-z]+)_(rep([0-9]{3})|pooled)_u([0-9]{3})\\.tsv$", t$file))[[1]]
    mech <- m[2]; u <- as.integer(m[5]) / 100
    repl <- if (m[3] == "pooled") NA_integer_ else as.integer(m[4])
    sp   <- SPECS[[sprintf("%s_%03d", mech, round(100 * u))]]
    args <- rdrArgs(t$model, mech, sp)
    path <- file.path(dataDir, t$file)
    if (LATENT[[t$model]]){                          # latent-parent scripts take Yo1 and the 4 PGS haplotypes only
        p5 <- file.path(tmpDir, paste0(t$model, "_", t$file))
        fwrite(fread(path, select = COLS5), p5, sep = "\t")
        path <- p5
    }
    truth <- sp$truth[[t$model]]
    t0 <- Sys.time()
    fit <- fitFree(RDR_BUILDERS[[t$model]], path, args, truth)
    id  <- data.frame(model = t$model, mechanism = mech, rep = repl, pooled = is.na(repl), u = u, f = F_TRUE,
                      n = nrow(fread(path, select = 1L)),
                      h2_RDR_parent = if (USES_PARENT_RDR[[t$model]]) sp$h2_parent else NA, h2_RDR_child = sp$h2_child)
    if (LATENT[[t$model]]) unlink(path)
    if (is.null(fit)){
        res <- cbind(id, parameter = NA, estimate = NA, SE = NA, true = NA, status = NA, m2LL = NA)
    } else {
        p <- summary(fit)$parameters
        res <- cbind(id[rep(1, nrow(p)), ], parameter = p$name, estimate = p$Estimate, SE = p$Std.Error,
                     true = unname(truth[p$name]), status = fit$output$status$code, m2LL = m2LL(fit))
    }
    res$time_sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    fwrite(res, t$out)
    cat(sprintf("%-7s %-26s f=%.3f status=%s (%s s)\n", t$model, sub("\\.tsv$", "", t$file),
                if (is.null(fit)) NA else res$estimate[res$parameter == "f11"], res$status[1], res$time_sec[1]))
    invisible(t$out)
}

todo <- which(!file.exists(tasks$out))
cat(nrow(tasks), "fits (", length(files), "datasets x", length(RDR_MODELS), "models ),", length(todo), "to do;", workers, "workers\n")
t0 <- Sys.time()
invisible(mclapply(todo, function(k) tryCatch(fitOne(k), error = function(e) cat(tasks$out[k], "ERROR:", conditionMessage(e), "\n")),
                   mc.cores = workers, mc.preschedule = FALSE))
cat("Done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

all <- rbindlist(lapply(list.files(fitDir, "\\.csv$", full.names = TRUE), fread), fill = TRUE)
fwrite(all, file.path(baseDir, outRoot, "rdr_fits.csv"))
cat("Wrote", nrow(all), "rows from", nrow(unique(all[, .(model, mechanism, rep, u)])), "fits to", file.path(outRoot, "rdr_fits.csv"), "\n\n")
fr <- all[parameter == "f11" & pooled == FALSE]
print(dcast(fr[, .(mean_f = round(mean(estimate), 4)), by = .(model, mechanism, u)], mechanism + u ~ model, value.var = "mean_f"))
