## Tests for the four Cascade Model 1 scripts (Balbona, Kim & Keller 2021, Model 1: the PGS explains all heritability,
## a = 0; only the offspring phenotype and the parental haplotypic PGS are observed):
##   Equilibrium/UniSEMPGS_Cascade_Model1.R                 Disequilibrium/UniSEMPGS_Cascade_DisEq_Model1.R
##   EquilibriumBinary/UniSEMPGS_Cascade_Binary_Model1.R    DisequilibriumBinary/UniSEMPGS_Cascade_DisEq_Binary_Model1.R
##
## A. Coding check, as in 02_fit_models.R --mode mvn: for each regime x mating mechanism, the Cascade iterative math
##    is run at BASE with all genetic variance moved into the PGS (delta^2 = .5, a = 0) and rescaled to parental
##    VY = 1; one exact-covariance dataset is drawn (continuous n = 32,000; binary n = 50,000 with Yo1 dichotomized
##    at the parent threshold) and fit with the script's DEFAULT search (extraTries = 30, exhaustive = F). With
##    latent parents both multipliers are fixed at the true mechanism.
## B. With AM_G = AM_E = 1 the Cascade scripts must reproduce OpenMxScripts/{Equilibrium,Disequilibrium}/*Model1.R:
##    both are fit to the primary-AM continuous datasets of A and their -2LL and estimates compared.
## C. Convergence of the default search on noisy data: the GeneEvolve SameTrait_Latent trios (output/data/, from
##    01_simulate.py) are fit with the default search and with exhaustive = T. These populations have a latent
##    genetic score (a > 0), so Model 1 is misspecified there and only the -2LL comparison is meaningful.
##    Skipped when output/data/ is absent.
## Run with: Rscript 08_model1.R [--workers 3]   (outputs in output/model1/)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
cascadeDir <- dirname(baseDir)
origDir    <- file.path(dirname(cascadeDir), "OpenMxScripts")
getArg <- function(flag, default){
    a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
    if (is.na(i) || i == length(a)) default else a[i + 1]
}
workers <- as.integer(getArg("--workers", 3))

suppressPackageStartupMessages({ library(OpenMx); library(data.table); library(MASS); library(parallel) })
# cap the threads the fitting scripts request, as in 02_fit_models.R
threads <- max(1, parallel::detectCores() %/% workers)
assignInNamespace("detectCores", function(...) threads, ns = "parallel")
source(file.path(baseDir, "conditions.R"))
cascadeValidationSetup(cascadeDir)
source(file.path(cascadeDir, "Equilibrium", "UniSEMPGS_Cascade_Model1.R"))
source(file.path(cascadeDir, "Disequilibrium", "UniSEMPGS_Cascade_DisEq_Model1.R"))
source(file.path(cascadeDir, "EquilibriumBinary", "UniSEMPGS_Cascade_Binary_Model1.R"))
source(file.path(cascadeDir, "DisequilibriumBinary", "UniSEMPGS_Cascade_DisEq_Binary_Model1.R"))
source(file.path(origDir, "Equilibrium", "UniSEMPGS_Model1.R"))
source(file.path(origDir, "Disequilibrium", "UniSEMPGS_DisEq_Model1.R"))

outDir  <- file.path(baseDir, "output", "model1"); dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
simDir  <- file.path(baseDir, "output", "data")
tol_cont <- .02; tol_bin <- .05

BASE_MODEL1 <- modifyList(BASE, list(delta = sqrt(.5), a = 0))   # BASE's .5 genetic variance, all in the PGS
fitFn <- function(regime, binary) get(paste0("fitUniSEMPGS_Cascade_", if (regime == "DisEq") "DisEq_" else "",
                                             if (binary) "Binary_" else "", "Model1"))
summaryTable <- function(s, truth){
    p <- s$parameters
    data.frame(parameter = p$name, estimate = p$Estimate, SE = p$Std.Error, true = unname(truth[p$name]),
               statusCode = as.character(s$statusCode)[1], minus2LL = s$Minus2LogLikelihood, stringsAsFactors = FALSE)
}

## ================= A. Coding check =================
specModel1 <- function(regime, mech, binary){
    q <- conditionTruth(regime, mech, base = BASE_MODEL1)$q
    VYo <- if (regime == "Eq") q$VY else q$VY_off
    gc  <- if (regime == "Eq") q$gc else 0
    CM <- rbind(c(VYo, q$thetaT, q$thetaNT, q$thetaT, q$thetaNT),
                c(q$thetaT,  q$k + gc, gc, q$gt, q$gt),
                c(q$thetaNT, gc, q$k + gc, q$gt, q$gt),
                c(q$thetaT,  q$gt, q$gt, q$k + gc, gc),
                c(q$thetaNT, q$gt, q$gt, gc, q$k + gc))
    dimnames(CM) <- list(c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1"), c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1"))
    truth <- c(VY11 = q$VY, VE11 = q$VE, delta11 = q$delta, Omega11 = q$Omega, mu11 = q$mu, f11 = q$f)
    if (regime == "Eq") truth <- c(truth, gt11 = q$gt, gc11 = q$gc, w11 = q$w, VF11 = q$VF)
    if (binary) truth <- c(truth[names(truth) != "VY11"], thresh_Yo11 = THRESH[["parent"]])
    truth <- c(truth, setNames(rep(0, 5), c("meanYo1", "meanTp1", "meanNTp1", "meanTm1", "meanNTm1")))
    list(CM = CM, truth = truth, mult = multiplierArgs(mech, latent = TRUE))
}

tasks <- expand.grid(regime = REGIMES, mech = names(MECHANISMS), binary = c(TRUE, FALSE), stringsAsFactors = FALSE)
fitA <- function(i){
    t <- tasks[i, ]
    tag  <- sprintf("%s_%s_%s", t$regime, t$mech, if (t$binary) "binary" else "continuous")
    spec <- specModel1(t$regime, t$mech, t$binary)
    stopifnot(max(abs(spec$CM - t(spec$CM))) < 1e-10, min(eigen(spec$CM, only.values = TRUE)$values) > 0)
    set.seed(2026)
    dat <- as.data.frame(mvrnorm(n = if (t$binary) 5e4 else 32000, mu = rep(0, 5), Sigma = spec$CM, empirical = TRUE))
    if (t$binary) dat$Yo1 <- as.integer(dat$Yo1 > THRESH[["parent"]])
    tsv <- file.path(outDir, paste0("mvn_", tag, ".tsv")); fwrite(dat, tsv, sep = "\t")
    t0 <- Sys.time()
    res <- tryCatch(summaryTable(do.call(fitFn(t$regime, t$binary), c(list(data_path = tsv), spec$mult)), spec$truth),
                    error = function(e) data.frame(parameter = NA, estimate = NA, SE = NA, true = NA,
                                                   statusCode = paste("ERROR:", conditionMessage(e)), minus2LL = NA))
    cbind(regime = t$regime, mechanism = t$mech, trait = if (t$binary) "binary" else "continuous", res,
          time_sec = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
A <- rbindlist(mclapply(seq_len(nrow(tasks)), fitA, mc.cores = workers, mc.preschedule = FALSE), fill = TRUE)
fwrite(A, file.path(outDir, "A_coding_check.csv"))

cat("\n================ A. Coding check (exact-covariance data, default search) ================\n")
for (key in unique(A[, paste(regime, mechanism, trait)])){
    d <- A[paste(regime, mechanism, trait) == key]
    struct <- !grepl("^mean", d$parameter)
    tol <- if (d$trait[1] == "binary") tol_bin else tol_cont
    d[, z := (estimate - true) / SE]
    bad <- struct & (is.na(d$true) | abs(d$estimate - d$true) > tol)
    cat(sprintf("%-32s status=%-4s max|est-true|=%.1e  max|z|=%5.2f  %4.0fs  %s\n", key, d$statusCode[1],
                max(abs(d$estimate - d$true)[struct]), max(abs(d$z)[struct]), d$time_sec[1],
                if (any(bad) || is.na(d$minus2LL[1])) paste("FLAG:", paste(d$parameter[bad], collapse = " ")) else "PASS"))
}

## ================= B. Reduction to the original Model 1 scripts (AM_G = AM_E = 1) =================
cat("\n================ B. Cascade (AM_G = AM_E = 1) vs original Model 1, primary-AM continuous data ================\n")
B <- rbindlist(lapply(REGIMES, function(rg){
    tsv  <- file.path(outDir, sprintf("mvn_%s_primary_continuous.tsv", rg))
    orig <- if (rg == "Eq") fitUniSEMPGS_Model1(tsv) else fitUniSEMPGS_DisEq_Model1(tsv)
    casc <- A[A$regime == rg & A$mechanism == "primary" & A$trait == "continuous"]
    p <- orig$parameters
    m <- match(p$name, casc$parameter)
    data.frame(regime = rg, parameter = p$name, original = p$Estimate, cascade = casc$estimate[m],
               m2LL_original = orig$Minus2LogLikelihood, m2LL_cascade = casc$minus2LL[1])
}))
fwrite(B, file.path(outDir, "B_reduction.csv"))
for (rg in REGIMES){
    d <- B[B$regime == rg]
    cat(sprintf("%-6s -2LL original %.4f  cascade %.4f  | max |estimate difference| over %d shared parameters = %.1e\n",
                rg, d$m2LL_original[1], d$m2LL_cascade[1], nrow(d), max(abs(d$original - d$cascade))))
}

## ================= C. Default vs exhaustive search on GeneEvolve trios =================
cat("\n================ C. Default vs exhaustive search on GeneEvolve trios (Model 1 misspecified: a > 0) ================\n")
simTasks <- expand.grid(regime = REGIMES, mech = names(MECHANISMS), stringsAsFactors = FALSE)
simTasks$path <- file.path(simDir, sprintf("%s_%s_SameTrait_Latent.tsv", simTasks$regime, simTasks$mech))
simTasks <- simTasks[file.exists(simTasks$path), ]
if (nrow(simTasks) == 0){
    cat("output/data/ not found: run 01_simulate.py first. Skipped.\n")
} else {
    fitC <- function(i){
        t <- simTasks[i, ]
        args <- c(list(data_path = t$path), multiplierArgs(t$mech, latent = TRUE))
        d <- do.call(fitFn(t$regime, FALSE), args)
        x <- do.call(fitFn(t$regime, FALSE), c(args, list(exhaustive = TRUE)))
        data.frame(regime = t$regime, mechanism = t$mech, n = d$numObs,
                   m2LL_default = d$Minus2LogLikelihood, m2LL_exhaustive = x$Minus2LogLikelihood,
                   status_default = as.character(d$statusCode)[1])
    }
    C <- rbindlist(mclapply(seq_len(nrow(simTasks)), fitC, mc.cores = workers, mc.preschedule = FALSE))
    C[, diff := m2LL_default - m2LL_exhaustive]
    fwrite(C, file.path(outDir, "C_default_vs_exhaustive.csv"))
    print(C, digits = 9)
}
