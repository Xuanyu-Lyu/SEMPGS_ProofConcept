## Shared setup for the G x Age sensitivity check (sourced by 02_population_fits.R and 03_fit_sim.R).
##
## Question: parents are measured as adults and offspring as children. If the two ages share only part of their
## genetic component, does Model 2 (UniSEMPGS_Cascade_SameTrait_ObservedYPYM.R, Equilibrium), which assumes one delta
## and one a for both generations, still estimate f (vertical transmission) without much bias?
##
## G x Age model (common factor, symmetric). With VG0 = 2*delta^2*k + 2*a^2*j the base-population genetic variance,
##   G_adult = G_shared + U_adult,   G_child = G_shared + U_child,   var(U_adult) = var(U_child) = u * VG0,
## so the genetic correlation between adult and child expression is r = 1 - u, and both ages keep variance VG0.
## The PGS carries the adult effects (adult GWAS), spouses assort on the adult mating phenotype, and F = f*(Yp + Ym)
## uses the parents' adult phenotypes. Only the offspring's phenotype, measured in childhood, uses the child value:
##   Yo = r*[delta*(Tp + Tm) + a*(LTp + LTm)] + sqrt(1 - r^2)*G_U + F + E,
## where G_U, the child-specific genetic value (variance VG0), is uncorrelated with the PGS, the LGS, F and the parents.
## Yo's moments are therefore those of the DiffTrait offspring trait with delta_o = r*delta, a_o = r*a and
## VE_o = VE + (1 - r^2)*VG0. 01_simulate.py implements the same model in GeneEvolve, with child effects
## beta_child = r*beta_adult + sqrt(1 - r^2)*beta_U at the same causal variants.

## BASE, conditionTruth, designSpec (Validation) and conditionBase, multipliers, BUILDERS, testF (FTestValidity)
source(file.path(if (exists("baseDir")) baseDir else getwd(), "..", "FTestValidity", "conditions.R"))
## Model 2 as shipped (fitUniSEMPGS_Cascade_SameTrait_ObservedYPYM), for the local-optimum check
source(file.path(REPO, scriptPath[["M2_Cascade"]]))

U_VALUES <- c(0, .1, .2, .3, .4, .5)       # unique share of each generation's genetic variance (r = 1 - u)
F_TRUE   <- .15                            # f in the GeneEvolve populations (BASE)

## Equilibrium truth for one mating type x f x u: the 7x7 population covariance of Yp1, Ym1, Yo1, Tp1, NTp1, Tm1, NTm1
## with Yo the child phenotype, the parents' values by Model 2 label (the truth Model 2 should recover), and the
## offspring moments. mating "none" is random mating (am = 0); latent = TRUE keeps the BASE PGS / LGS split.
gxageSpec <- function(mating, f, u){
    mech <- if (mating == "none") "primary" else mating
    tr   <- conditionTruth("Eq", mech, base = conditionBase(latent = TRUE, f = f, mating = mating))
    q    <- tr$q
    r    <- 1 - u
    same <- designSpec(tr, "SameTrait_ObservedYPYM", binary = FALSE)
    VG0  <- 2*q$k*q$delta^2 + 2*q$j*q$a^2
    tr$o <- diffTraitOffspring_Cascade(q, delta_o = r*q$delta, a_o = r*q$a, VE_o = q$VE + (1 - r^2)*VG0)
    CM   <- designSpec(tr, "DiffTrait_ObservedYPYM_EstimatedAParent", binary = FALSE)$CM
    if (u == 0) stopifnot(max(abs(CM - same$CM)) < 1e-10)      # no G x Age: the SameTrait covariance
    stopifnot(min(eigen(CM, only.values = TRUE)$values) > 0)
    list(CM = CM, truth = same$truth, q = q, o = tr$o, VG0 = VG0, r = r, tr = tr)
}
