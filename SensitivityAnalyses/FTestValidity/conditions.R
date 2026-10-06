## Shared setup for the f-test validity check (sourced by 01_population_lrt.R and 02_monte_carlo.R).
##
## Claim under test: the likelihood-ratio test of f = 0 (vertical transmission) gives the same answer whether one
## fits Model 2 or Model 1, whatever the type of assortative mating (AM) and whether or not the trait has a latent
## genetic component that the PGS does not capture (which violates Model 1's assumption).
##
## Equilibrium, continuous trait. The data are generated from the SEM-PGS Cascade model, which is the true model for
## every condition, and three models are fit to each dataset:
##   M2_Cascade  Model 2, Cascade version (UniSEMPGS_Cascade_SameTrait_ObservedYPYM.R): parental phenotypes observed,
##               latent genetic score modelled, mating mechanism estimated. Correctly specified everywhere.
##   M1_Cascade  Model 1, Cascade version (UniSEMPGS_Cascade_Model1.R): offspring phenotype and PGS only, no latent
##               genetic score, mating mechanism fixed at the true one.
##   M1          Model 1 (UniSEMPGS_Model1.R): offspring phenotype and PGS only, no latent genetic score, primary
##               phenotypic AM assumed.
## For each model f is tested with mxCompare(free-f fit, f = 0 fit): a 1-df likelihood-ratio test.

findRepo <- function(dir){
    while (!file.exists(file.path(dir, "OpenMxScripts_Cascade"))) dir <- dirname(dir)
    dir
}
REPO <- findRepo(if (exists("baseDir")) baseDir else getwd())
source(file.path(REPO, "OpenMxScripts_Cascade", "Validation", "conditions.R"))   # BASE, MECHANISMS, conditionTruth, designSpec
cascadeValidationSetup(file.path(REPO, "OpenMxScripts_Cascade"))

## ---------------- conditions ----------------
## latent = FALSE: all base-population genetic variance (.5) is in the PGS (Model 1's assumption holds);
## latent = TRUE:  the BASE split, delta^2 = .2 in the PGS and a^2 = .3 in the latent genetic score.
## mating "none" is random mating (am = 0). f = 0 gives the type I error; f > 0 the power.
F_VALUES <- c(0, .05, .15)
MATING   <- c("none", "primary", "genetic", "social")
GRID <- expand.grid(mating = MATING, latent = c(FALSE, TRUE), f = F_VALUES, stringsAsFactors = FALSE)
GRID$cond <- with(GRID, sprintf("%s_%s_f%03d", mating, ifelse(latent, "latent", "nolatent"), round(100 * f)))
MODELS <- c("M2_Cascade", "M1_Cascade", "M1")

conditionBase <- function(latent, f, mating){
    b <- modifyList(BASE, list(f = f, am = if (mating == "none") 0 else BASE$am))
    if (!latent) b <- modifyList(b, list(delta = sqrt(.5), a = 0))
    b
}
## Truth (parental VY rescaled to 1) and the 7x7 population covariance of Yp1, Ym1, Yo1, Tp1, NTp1, Tm1, NTm1
conditionSpec <- function(row){
    mech <- if (row$mating == "none") "primary" else row$mating
    tr   <- conditionTruth("Eq", mech, base = conditionBase(row$latent, row$f, row$mating))
    spec <- designSpec(tr, "SameTrait_ObservedYPYM", binary = FALSE)
    stopifnot(min(eigen(spec$CM, only.values = TRUE)$values) > 0)
    spec$q <- tr$q
    spec
}

## Mating multipliers passed to the Cascade models. Model 2 estimates the mechanism (conditions.R's convention:
## AM_E fixed at 1 and AM_G free, or AM_G fixed at 1 and AM_E free for genetic homogamy); Model 1 Cascade has latent
## parents, so its mechanism is fixed at the truth. Under random mating mu = 0 and the multipliers are not
## identified, so both are fixed at 1.
multipliers <- function(model, mating){
    if (mating == "none") return(list(AM_G_value = 1, AM_G_free = FALSE, AM_E_value = 1, AM_E_free = FALSE))
    multiplierArgs(mating, latent = model == "M1_Cascade")
}

## ---------------- model builders ----------------
## Each fitting script ends with mxTryHard(Model1, ...) and return(summary(...)). The builders below are the same
## functions with those two lines replaced, so they return the UNFITTED MxModel and the harness controls the fit.
scriptPath <- c(M2_Cascade = "OpenMxScripts_Cascade/Equilibrium/UniSEMPGS_Cascade_SameTrait_ObservedYPYM.R",
                M1_Cascade = "OpenMxScripts_Cascade/Equilibrium/UniSEMPGS_Cascade_Model1.R",
                M1         = "OpenMxScripts/Equilibrium/UniSEMPGS_Model1.R")
loadBuilder <- function(model){
    src <- paste(readLines(file.path(REPO, scriptPath[[model]])), collapse = "\n")
    fit <- regmatches(src, regexpr("fitModel1 <- mxTryHard\\(Model1,[^)]*c\\(0,1\\)[^)]*\\)", src))
    stopifnot(length(fit) == 1, grepl("return\\(summary\\(fitModel1, verbose = TRUE\\)\\)", src))
    src <- sub(fit, "fitModel1 <- Model1", src, fixed = TRUE)
    src <- sub("return(summary(fitModel1, verbose = TRUE))", "return(fitModel1)", src, fixed = TRUE)
    env <- new.env()
    eval(parse(text = src), envir = env)
    get(ls(env)[1], envir = env)
}
BUILDERS <- lapply(setNames(MODELS, MODELS), loadBuilder)

## ---------------- data files ----------------
COLS7 <- c("Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1")
COLS5 <- c("Yo1", "Tp1", "NTp1", "Tm1", "NTm1")
writeData <- function(dat, stem){
    p7 <- paste0(stem, "_7.tsv"); p5 <- paste0(stem, "_5.tsv")
    data.table::fwrite(as.data.frame(dat)[, COLS7], p7, sep = "\t")
    data.table::fwrite(as.data.frame(dat)[, COLS5], p5, sep = "\t")
    c(M2_Cascade = p7, M1_Cascade = p5, M1 = p5)
}

## ---------------- the likelihood-ratio test of f = 0 ----------------
## 1. free fit: mxTryHard from the model's starting values (M2_Cascade: the generating values, as is usual in Monte
##    Carlo studies; Model 1: the scripts' own sample-moment starts), and again from f = -.2, 0, .1, .2 and .4 (the
##    Model 1 Cascade likelihood under social homogamy has several optima, so one start is not enough); the best
##    fit is kept;
## 2. f = 0 fit: mxTryHard from the best free solution with f fixed at 0, and from the original starting values;
##    the better one is kept;
## 3. if the f = 0 fit beats the free fit, or the free fit ends with a status other than 0 or 1, the free model is
##    refit from the f = 0 solution (with f starting at exactly 0, so it cannot end worse than the f = 0 fit).
## These steps guarantee -2LL(free) <= -2LL(f = 0) and guard both fits against local optima, which would otherwise
## show up as spurious rejections, missed rejections or negative statistics.
## Returns the free-model estimate and SE of f and the mxCompare() result.
F_STARTS <- c(-.2, 0, .1, .2, .4)
tryHard <- function(m, extraTries){
    out <- NULL
    invisible(capture.output(out <- suppressWarnings(suppressMessages(mxTryHard(m, extraTries = extraTries,
        OKstatuscodes = c(0, 1), silent = TRUE, exhaustive = FALSE, jitterDistrib = "rnorm", loc = .5, scale = .1)))))
    out
}
m2LL <- function(fit){            # a fit that ended without a -2LL counts as failed
    v <- if (is.null(fit)) NULL else fit$output$Minus2LogLikelihood
    if (length(v) == 0 || !is.finite(v)) Inf else v
}
better <- function(a, b) if (m2LL(b) < m2LL(a)) b else a
safeFit <- function(m, extraTries) tryCatch(tryHard(m, extraTries), error = function(e) NULL)
okStatus <- function(fit) is.finite(m2LL(fit)) && isTRUE(fit$output$status$code %in% c(0, 1))

testF <- function(model, data_path, mult, truth, extraTries = 5){
    args <- c(list(data_path = data_path), if (model != "M1") mult)
    m <- do.call(BUILDERS[[model]], args)
    if (model == "M2_Cascade"){
        lab <- intersect(names(truth), names(omxGetParameters(m)))
        lb  <- omxGetParameters(m, fetch = "lbound")[lab]
        start <- ifelse(is.na(lb), truth[lab], pmax(truth[lab], lb + .01))   # e.g. a = 0 sits below a's lbound
        m <- omxSetParameters(m, labels = lab, values = start)
    }
    noSE   <- function(x) mxOption(mxOption(x, "Calculate Hessian", "No"), "Standard Errors", "No")
    withSE <- function(x) mxOption(mxOption(x, "Calculate Hessian", "Yes"), "Standard Errors", "Yes")
    fix0   <- function(x) omxSetParameters(x, labels = "f11", free = FALSE, values = 0)
    m0   <- noSE(m)
    free <- safeFit(m0, extraTries)
    for (fs in F_STARTS) free <- better(free, safeFit(omxSetParameters(m0, labels = "f11", values = fs), extraTries))
    f0   <- if (!is.null(free)) safeFit(fix0(free), extraTries) else NULL
    f0   <- better(f0, safeFit(fix0(m0), extraTries))
    reseeded <- FALSE
    if (!is.null(f0) && (m2LL(f0) < m2LL(free) - 1e-6 || !okStatus(free))){
        free <- better(free, safeFit(omxSetParameters(f0, labels = "f11", free = TRUE, values = 0), extraTries))
        reseeded <- TRUE
    }
    if (!is.finite(m2LL(free)) || !is.finite(m2LL(f0))) return(data.frame(model = model, f_est = NA, f_SE = NA, m2LL_free = NA,
        m2LL_f0 = NA, LRT = NA, p = NA, status_free = NA, status_f0 = NA, reseeded = reseeded))
    # standard errors from the chosen free solution (a run that starts at the optimum)
    se <- tryCatch(suppressWarnings(mxRun(withSE(free), silent = TRUE, suppressWarnings = TRUE)), error = function(e) NULL)
    if (!is.null(se) && m2LL(se) <= m2LL(free) + 1e-6) free <- se
    cmp <- mxCompare(free, f0)
    s <- summary(free)
    data.frame(model = model,
               f_est = s$parameters$Estimate[s$parameters$name == "f11"],
               f_SE  = s$parameters$Std.Error[s$parameters$name == "f11"],
               m2LL_free = m2LL(free), m2LL_f0 = m2LL(f0),
               LRT = cmp$diffLL[2], p = cmp$p[2],
               status_free = free$output$status$code, status_f0 = f0$output$status$code, reseeded = reseeded)
}
