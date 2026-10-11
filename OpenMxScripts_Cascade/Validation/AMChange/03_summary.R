## Step 3: summarize the AMChange validation (writes output/results/summary_*.csv and prints the tables).
##   A. Moments: GeneEvolve vs the iterative math (the scripts' formulas) vs PDF Part VII Section 38 as printed.
##      Data-level moments are computed per population from the exported trios (exactly the entries of the scripts'
##      covariance matrices, one random offspring per family); the offspring's F moments come from 01_simulate.py's
##      population records. mean and SE are over populations; z = (sim - math)/SE.
##   B. Coding check (02_fit_models.R --mode mvn): max |estimate - truth| per fit.
##   C. GeneEvolve fits (--mode sim): z = (estimate - truth)/SE per parameter, max |z| per fit.
##   D. Reduction fits (--mode reduction): -2LL vs the original Eq / DisEq script, and the estimated mu_pp.
## Run with: Rscript 03_summary.R

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- sub("^--file=", "", grep("^--file=", cmdArgs, value = TRUE))
baseDir <- if (length(fileArg) > 0) dirname(normalizePath(fileArg)) else getwd()
cascadeDir <- dirname(dirname(baseDir))
suppressPackageStartupMessages(library(data.table))
source(file.path(baseDir, "conditions_AMChange.R"))
amchangeSetup(cascadeDir)
resDir <- file.path(baseDir, "output", "results"); dir.create(resDir, recursive = TRUE, showWarnings = FALSE)
partsDir <- file.path(baseDir, "output", "sim", "parts")
options(width = 200)

## ---------------- A. Moments ----------------
if (dir.exists(partsDir)){
    grid <- CJ(branch = names(BRANCHES), mechanism = names(MECHANISMS))
    simMoments <- rbindlist(lapply(seq_len(nrow(grid)), function(i){
        b <- grid$branch[i]; m <- grid$mechanism[i]
        files <- list.files(partsDir, sprintf("^%s_%s_rep[0-9]+_SameTrait_Observed\\.tsv$", b, m), full.names = TRUE)
        rbindlist(lapply(files, function(fs){
            d  <- fread(fs); dd <- fread(sub("SameTrait_Observed", "DiffTrait_Observed", fs))
            S  <- cov(d); Sd <- cov(dd)
            hap <- function(X, a, b) mean(c(X[paste0(a, "p1"), paste0(b, "p1")], X[paste0(a, "m1"), paste0(b, "m1")]))
            data.table(branch = b, mechanism = m, rep = as.integer(sub(".*_rep([0-9]+)_.*", "\\1", basename(fs))),
                       VY_p = mean(c(S["Yp1", "Yp1"], S["Ym1", "Ym1"])),
                       Omega_p = mean(c(S["Yp1", "Tp1"], S["Yp1", "NTp1"], S["Ym1", "Tm1"], S["Ym1", "NTm1"])),
                       hapvar = mean(diag(S)[c("Tp1", "NTp1", "Tm1", "NTm1")]),
                       gc_p = hap(S, "T", "NT"),
                       gt_p = mean(S[c("Tp1", "NTp1"), c("Tm1", "NTm1")]),
                       Yp_Ym = S["Yp1", "Ym1"],
                       Yp_PGSm = mean(c(S["Yp1", c("Tm1", "NTm1")], S["Ym1", c("Tp1", "NTp1")])),
                       Yo_Yp = mean(S["Yo1", c("Yp1", "Ym1")]), VY_o = S["Yo1", "Yo1"],
                       thetaT = mean(S["Yo1", c("Tp1", "Tm1")]), thetaNT = mean(S["Yo1", c("NTp1", "NTm1")]),
                       VY_o_diff = Sd["Yo1", "Yo1"], Yo_Yp_diff = mean(Sd["Yo1", c("Yp1", "Ym1")]),
                       thetaT_diff = mean(Sd["Yo1", c("Tp1", "Tm1")]), thetaNT_diff = mean(Sd["Yo1", c("NTp1", "NTm1")]))
        }))
    }))
    simMoments <- melt(simMoments, id.vars = c("branch", "mechanism", "rep"), variable.name = "quantity", value.name = "value")
    # population-level moments of the parents' mating and of the offspring's F (from 01_simulate.py's records)
    mf <- fread(file.path(baseDir, "output", "sim", "moments_final.csv"))
    pop <- mf[quantity %in% c("mu", "r_gamma", "w_off", "VF_off", "Y_Fo")]
    pop[quantity == "mu", quantity := "mu_p"]; pop[quantity == "w_off", quantity := "w_o"]
    pop[quantity == "VF_off", quantity := "VF_o"]
    simMoments <- rbind(simMoments, pop[, .(branch, mechanism, rep, quantity, value)])
    simSum <- simMoments[, .(sim = mean(value), SE = sd(value) / sqrt(.N), reps = .N), by = .(branch, mechanism, quantity)]

    mathOf <- function(b, m){
        tr <- amchangeTruth(m, b); q <- tr$q; o <- tr$o; P <- q$pdf_partVII
        data.table(quantity = c("VY_p", "Omega_p", "hapvar", "gc_p", "gt_p", "Yp_Ym", "Yp_PGSm", "Yo_Yp", "VY_o",
                                "thetaT", "thetaNT", "VY_o_diff", "Yo_Yp_diff", "thetaT_diff", "thetaNT_diff",
                                "mu_p", "r_gamma", "w_o", "VF_o", "Y_Fo"),
                   math = c(q$VY_p, q$Omega_p, q$k + q$gc_p, q$gc_p, q$gt_p, q$Yp_Ym, q$Yp_PGSm, q$Yo_Yp, q$VY_o,
                            q$thetaT, q$thetaNT, o$VY_o, o$Yo_Yp_o, o$thetaT_o, o$thetaNT_o,
                            q$mu_p, q$am_p, q$w_o, q$VF_o, q$Y_Fo),
                   # PDF Part VII Section 38 as printed (per data column: half the papers' two-haplotype theta; the
                   # PDF has no gc term in the gamete variance, i.e. var(T) = k)
                   pdf38 = c(NA, NA, q$k, 0, NA, NA, NA, NA, P$VY2, P$thetaT2 / 2, P$thetaNT2 / 2,
                             rep(NA, 4), NA, NA, P$w2, P$VF2, NA))
    }
    math <- rbindlist(lapply(seq_len(nrow(grid)), function(i) cbind(grid[i], mathOf(grid$branch[i], grid$mechanism[i]))))
    A <- merge(simSum, math, by = c("branch", "mechanism", "quantity"))
    A[, `:=`(z = (sim - math) / SE, rel = (sim - math) / abs(math), z_pdf38 = (sim - pdf38) / SE)]
    A[, branch := factor(branch, names(BRANCHES))]; A[, mechanism := factor(mechanism, names(MECHANISMS))]
    setorder(A, branch, mechanism, quantity)
    fwrite(A, file.path(resDir, "summary_moments.csv"))
    cat("A. Moments: GeneEvolve (mean over", max(A$reps), "populations) vs the scripts' math\n")
    print(A[, .(n = .N, max_abs_z = round(max(abs(z)), 2), worst = quantity[which.max(abs(z))],
                max_rel_pct = round(100 * max(abs(rel[abs(math) > .01])), 1)), by = .(branch, mechanism)])
    cat("\n   The quantities PDF Part VII Section 38 also predicts (z = (sim - prediction)/SE):\n")
    print(dcast(A[!is.na(pdf38) & quantity %in% c("thetaNT", "thetaT", "VY_o", "hapvar", "gc_p")],
                branch + mechanism ~ quantity, value.var = c("z", "z_pdf38"), fun.aggregate = function(x) round(x, 1)))
}

## ---------------- B. Coding check ----------------
summarizeFits <- function(mode, E = NULL){
    f <- file.path(resDir, paste0("estimates_", mode, ".csv"))
    if (is.null(E)){
        if (!file.exists(f) || file.size(f) == 0) return(NULL)
        E <- fread(f)
    }
    E[, struct := !is.na(true) & !grepl("^mean", parameter)]
    E[, .(status = statusCode[1], n = n[1], max_abs_err = max(abs(estimate - true)[struct]),
          worst = parameter[struct][which.max(abs(estimate - true)[struct])],
          max_abs_z = max(abs((estimate - true) / SE)[struct], na.rm = TRUE),
          worst_z = parameter[struct][which.max(abs((estimate - true) / SE)[struct])],
          big_miss_no_SE = any(struct & abs(estimate - true) > .05 & is.na(SE)),
          m2LL = minus2LL[1], m2LL_original = m2LL_original[1], mu_pp = estimate[parameter == "mu_pp11"][1],
          mu_pp_true = true[parameter == "mu_pp11"][1], time_sec = time_sec[1]),
      by = .(branch, mechanism, design, trait)]
}
B <- summarizeFits("mvn")
if (!is.null(B)){
    Est <- fread(file.path(resDir, "estimates_mvn.csv"))
    if (!"seed" %in% names(Est)) Est[, seed := 2026L]
    Est[is.na(seed), seed := 2026L]
    Est[, struct := !is.na(true) & !grepl("^mean", parameter)]
    # fits with replicate datasets (--seed), together with their seed-2026 original
    repKeys <- Est[seed != 2026L, unique(paste(branch, mechanism, design, trait))]
    reps <- Est[struct & paste(branch, mechanism, design, trait) %in% repKeys]
    B <- summarizeFits("mvn", Est[seed == 2026L])
    # continuous data reproduce the implied matrix exactly, so the estimates must hit the truth; dichotomized data
    # carry sampling noise, so binary fits are judged by z (or a small absolute miss)
    B[, result := ifelse(trait == "continuous", ifelse(max_abs_err <= .02, "PASS", "FLAG"),
                         ifelse(max_abs_err <= .05 | (max_abs_z <= 3 & !big_miss_no_SE), "PASS", "FLAG"))]
    fwrite(B, file.path(resDir, "summary_mvn.csv"))
    cat("\nB. Coding check (exact-covariance data): PASS = max |estimate - truth| <= .02 (continuous); <= .05, or max |z| <= 3\n",
        "   with no parameter missing by > .05 without an SE (binary)\n")
    print(B[, .(branch, mechanism, design, trait, status, max_abs_err = signif(max_abs_err, 2), worst,
                max_abs_z = round(max_abs_z, 2), result, time_sec)])
    if (nrow(reps) > 0){
        reps[, z := (estimate - true) / SE]
        R <- reps[, .(datasets = uniqueN(seed), mean_z = round(mean(z), 2), sd_z = round(sd(z), 2),
                      max_abs_z = round(max(abs(z)), 2)), by = .(branch, mechanism, design, trait, parameter)]
        fwrite(R, file.path(resDir, "summary_mvn_replicates.csv"))
        cat("\n   Replicate datasets (--seed) of a binary coding check: z per parameter across datasets\n")
        print(R)
    }
}

## ---------------- C. GeneEvolve fits ----------------
## Likelihood-ratio test of the truth (continuous fits): -2LL with every structural parameter at its true value (the
## means at their MLE, the sample means) minus the fit's -2LL, against chi-square with df = the free structural
## dimensions (free parameters minus constraints). It does not rely on the Hessian SEs, which are unreliable for
## some of these fits (NaN, or too small; see the README).
m2LLatTruth <- function(CM, path){
    X <- as.matrix(fread(path)[, colnames(CM), with = FALSE]); n <- nrow(X); S <- cov(X) * (n - 1) / n
    n * (ncol(X) * log(2 * pi) + as.numeric(determinant(CM)$modulus) + sum(diag(solve(CM, S))))
}
structuralDims <- function(design, mechFixed){
    dims <- c(SameTrait_ObservedYPYM = 7, SameTrait_LatentYPYM = 5, DiffTrait_LatentYPYM = 5,
              DiffTrait_ObservedYPYM_EstimatedAParent = 9, DiffTrait_ObservedYPYM_FixedAParent = 8, Model1 = 5)[[design]]
    dims - (mechFixed && grepl("Observed", design))          # observed designs estimate one multiplier
}
C <- summarizeFits("sim")
if (!is.null(C)){
    C[, c("d_m2LL_truth", "df", "p_truth") := {
        out <- lapply(seq_len(.N), function(i){
            if (trait[i] != "continuous") return(list(NA_real_, NA_real_, NA_real_))
            tr <- amchangeTruth(mechanism[i], branch[i]); sp <- designSpecAMChange(tr, design[i], FALSE)
            path <- file.path(baseDir, "output", "data", sprintf("%s_%s_%s.tsv", branch[i], mechanism[i], sp$dataset))
            d  <- m2LLatTruth(sp$CM, path) - m2LL[i]
            df <- structuralDims(design[i], !isTRUE(sp$args$AM_G_free) && !isTRUE(sp$args$AM_E_free))
            list(d, df, pchisq(d, df, lower.tail = FALSE))
        })
        list(sapply(out, `[[`, 1), sapply(out, `[[`, 2), sapply(out, `[[`, 3))
    }]
    fwrite(C, file.path(resDir, "summary_sim.csv"))
    cat("\nC. Fits to the pooled GeneEvolve trios: max |z| = |estimate - truth| / SE per fit, and the LR test of the truth\n")
    print(C[, .(branch, mechanism, design, trait, status, n, max_abs_z = round(max_abs_z, 2), worst_z,
                max_abs_err = signif(max_abs_err, 2), d_m2LL_truth = round(d_m2LL_truth, 2), df,
                p_truth = signif(p_truth, 2), time_sec)])
    E <- fread(file.path(resDir, "estimates_sim.csv"))[!is.na(true) & !grepl("^mean", parameter)]
    E[, z := (estimate - true) / SE]
    cat("\n   mu_pp and mu_p across the GeneEvolve fits:\n")
    print(E[parameter %in% c("mu_pp11", "mu_p11"), .(branch, mechanism, design, trait, parameter,
                                                     estimate = round(estimate, 3), true = round(true, 3), z = round(z, 2))])
}

## ---------------- D. Reduction ----------------
D <- summarizeFits("reduction")
if (!is.null(D)){
    fwrite(D, file.path(resDir, "summary_reduction.csv"))
    cat("\nD. AMChange scripts on the Eq / DisEq exact-covariance data (mu_pp must be mu_p for Eq and 0 for DisEq)\n")
    print(D[, .(regime = branch, mechanism, design, status, d_m2LL = signif(m2LL - m2LL_original, 2),
                mu_pp = signif(mu_pp, 4), mu_pp_true = signif(mu_pp_true, 4), max_abs_err = signif(max_abs_err, 2), worst)])
}
