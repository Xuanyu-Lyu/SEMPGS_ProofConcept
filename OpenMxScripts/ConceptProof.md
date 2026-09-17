# Proof of Concept: Univariate SEM-PGS Model Variants

This document describes the twenty univariate SEM-PGS fitting functions in this directory — five
model designs, each fit under two regimes (**equilibrium**, `Equilibrium/`, vs **disequilibrium**,
`Disequilibrium/`) and two trait types (**continuous**, fit as an ordinary phenotype, vs **binary**,
fit as a liability-threshold trait in `EquilibriumBinary/`/`DisequilibriumBinary/`) — how each model
is designed and identified, how the parameter-recovery simulations are constructed, and the test
results.

All models are scalar (univariate) versions of the bivariate SEM-PGS framework: every matrix in the
bivariate math is treated as a scalar. Fitting is done in OpenMx (NPSOL, raw-data ML, nonlinear
equality constraints via `mxConstraint`, `mxTryHard` with jittered restarts).

## 1. Notation

| Symbol | Meaning |
|---|---|
| `delta`, `a` | effects of the observed PGS and the latent (unmeasured) genetic score on the phenotype |
| `k`, `j` | haplotypic variances of the observed / latent PGS (fixed at .5) |
| `f` | vertical transmission (VT) path from each parental phenotype to the offspring's familial environment |
| `mu` | assortative-mating (AM) co-path; spousal covariance = `VY*mu*VY` |
| `Omega`, `Gamma` | within-person covariance between the phenotype and the observed / latent haplotypic PGS |
| `gt`, `ht` | between-spouse haplotype covariances induced by AM (observed–observed, latent–latent) |
| `gc`, `hc`, `ic` | within-person between-haplotype covariances (inherited from AM in the previous generation) |
| `w`, `v` | covariance between the familial environment and the observed / latent PGS |
| `VF`, `VE`, `VY` | familial-environment, residual-environment, phenotypic (liability, for binary traits) variance |
| `thetaT`, `thetaNT` | covariance of the offspring phenotype with a transmitted / non-transmitted parental haplotype |
| `Th`, threshold | (binary models only) the liability cutpoint(s) separating the two observed categories of a dichotomous phenotype |

Observed variables: parental phenotypes `Yp1`, `Ym1` (when observed), offspring phenotype `Yo1`, and
the transmitted / non-transmitted haplotypic PGS from each parent `Tp1`, `NTp1`, `Tm1`, `NTm1`. In
the binary variants, whichever of `Yp1`/`Ym1`/`Yo1` are observed are dichotomous (0/1, fit via
`mxFactor()` + a liability threshold); the PGS columns are always continuous.

`itlo`/`itol` (the transmitted-vs-non-transmitted haplotype covariances used in the bivariate model)
never appear in these univariate scripts: for scalars, `Gamma*mu*Omega == Omega*mu*Gamma`, so both
collapse into the single quantity `ic` (equilibrium) or are computed directly as part of `thetaNT`
(disequilibrium) — see §2.

## 2. The two regimes

### Equilibrium (AM + VT at their joint fixed point)

The population has experienced AM and VT for many generations, so all moments satisfy the
recursion's fixed point:

```
Omega = 2*delta*gc + 2*a*ic + delta*k + .5*w
Gamma = 2*a*hc + 2*delta*ic + a*j + .5*v
VY    = 2*delta*Omega + 2*a*Gamma + delta*w + a*v + VF + VE
VF    = 2*f^2*VY + 2*f^2*mu*VY^2
w     = 2*f*Omega + 2*f*VY*mu*Omega ;  v = 2*f*Gamma + 2*f*VY*mu*Gamma
gt    = Omega^2*mu ; ht = Gamma^2*mu ; ic = Omega*Gamma*mu
gc    = gt ; hc = ht
thetaNT = 2*delta*gc + 2*a*ic + .5*w ;  thetaT = delta*k + thetaNT
Yo_Yp = (delta*Omega + a*Gamma)*(1 + mu*VY) + f*VY*(1 + mu*VY)
```

(For binary traits, `VY` is fixed to 1 instead of estimated — see §3's Identification note.)

### Disequilibrium (VT equilibrated under random mating; AM for the first time in the parents' generation)

`gc = hc = 0` (no AM in any earlier generation), all `mu` terms drop out of the within-generation
quantities, and the equations become self-referential VT fixed points:

```
Omega = delta*k + .5*w   with  w = 2*f*Omega     (=> Omega = delta*k/(1-f))
Gamma = a*j + .5*v       with  v = 2*f*Gamma
VY    = 2*delta*Omega + 2*a*Gamma + delta*w + a*v + VF + VE ,  VF = 2*f^2*VY
gt    = Omega^2*mu ; ht = Gamma^2*mu ; ic = Omega*Gamma*mu   (one generation of AM)
thetaNT = 2*a*ic + 2*delta*gt + w ;  thetaT = 2*delta*k + thetaNT
Yo_Yp = delta*Omega + a*Gamma + f*VY + (delta*Omega + a*Gamma)*mu*VY
```

In the disequilibrium models `gt/ht/ic/w/v` are pure `mxAlgebra`s (functions of the free
parameters), whereas in the equilibrium models they are free parameters tied down by
`mxConstraint`s.

## 3. The twenty models

Every model shares the same PGS block of the expected covariance matrix (`k+gc` on the diagonal,
`gc` within person, `gt` across spouses) and the `thetaT`/`thetaNT` column for the offspring
phenotype. They differ in (a) whether parental phenotypes are observed, (b) whether parent and
offspring express the same trait, (c) how `a` is identified, and (d) whether the phenotype(s) are
continuous or binary.

### Equilibrium (`Equilibrium/`, `EquilibriumBinary/`)

| # | Continuous script | Binary script | Data | Identification of `a` |
|---|---|---|---|---|
| 1 | `UniSEMPGS_SameTrait_ObservedYPYM.R` | `UniSEMPGS_Binary_SameTrait_ObservedYPYM.R` | 7 vars: Yp, Ym, Yo, Tp, NTp, Tm, NTm | from the observed parental phenotypes (no RDR needed) |
| 2 | `UniSEMPGS_SameTrait_LatentYPYM.R` | `UniSEMPGS_Binary_SameTrait_LatentYPYM.R` | 5 vars: Yo, Tp, NTp, Tm, NTm | RDR constraint: `(2a²j+2δ²k)(2a²j+2δ²k+VE) = h2_RDR·VY` |
| 3 | `UniSEMPGS_DiffTrait_ LatentYPYM.R` | `UniSEMPGS_Binary_DiffTrait_LatentYPYM.R` | 5 vars | RDR for parent **and** offspring `a`; parental and offspring VY constrained equal; needs `mu_fixed` (see Identification below) |
| 4 | `UniSEMPGS_DiffTrait_ ObservedYPYM_EstiatedAParent.R` | `UniSEMPGS_Binary_DiffTrait_ObservedYPYM_EstimatedAParent.R` | 7 vars | RDR for offspring `a` only; parental `a` estimated from the observed parental phenotypes |
| 5 | `UniSEMPGS_DiffTrait_ ObservedYPYM_FixedAParent.R` | `UniSEMPGS_Binary_DiffTrait_ObservedYPYM_FixedAParent.R` | 7 vars | RDR anchors both parental and offspring `a` |

The DiffTrait models carry generation-specific `delta_p, a_p, delta_o, a_o, VE_p, VE_o` (and, in
#4/#5, separate `VY_p`, `VY_o`). The offspring-trait variance under equilibrium is

```
VY_o = 2*delta_o^2*k + 4*delta_o^2*gc + 2*a_o^2*j + 4*a_o^2*hc + 8*delta_o*a_o*ic
     + 2*delta_o*w + 2*a_o*v + 2*f^2*VY_p*(1 + VY_p*mu) + VE_o
```

with `gc, hc, ic, w, v, mu` taken from the parental-trait equilibrium process; `thetaT`/`thetaNT`/`Yo_Yp` use `delta_o, a_o` in place of `delta, a`.

### Disequilibrium (`Disequilibrium/`, `DisequilibriumBinary/`)

| # | Continuous script | Binary script | Data | Identification of `a` |
|---|---|---|---|---|
| 6 | `UniSEMPGS_DisEq_SameTrait_ObservedYPYM.R` | `UniSEMPGS_DisEq_Binary_SameTrait_ObservedYPYM.R` | 7 vars | observed parental phenotypes |
| 7 | `UniSEMPGS_DisEq_SameTrait_LatentYPYM.R` | `UniSEMPGS_DisEq_Binary_SameTrait_LatentYPYM.R` | 5 vars | RDR |
| 8 | `UniSEMPGS_DisEq_DiffTrait_LatentYPYM.R` | `UniSEMPGS_DisEq_Binary_DiffTrait_LatentYPYM.R` | 5 vars | RDR for both generations; shared VY; needs `mu_fixed` (see Identification below) |
| 9 | `UniSEMPGS_DisEq_DiffTrait_ObservedYPYM_EstimatedAParent.R` | `UniSEMPGS_DisEq_Binary_DiffTrait_ObservedYPYM_EstimatedAParent.R` | 7 vars | RDR for offspring only |
| 10 | `UniSEMPGS_DisEq_DiffTrait_ObservedYPYM_FixedAParent.R` | `UniSEMPGS_DisEq_Binary_DiffTrait_ObservedYPYM_FixedAParent.R` | 7 vars | RDR for both generations |

Same layout as #1–#5 but with the disequilibrium equations of §2: `gc=hc=0` fixed, `gt/ht/ic/w/v` as
algebras, and the offspring-trait variance simplifying to

```
VY_o = 2*delta_o^2*k + 2*a_o^2*j + 2*delta_o*w + 2*a_o*v + 2*f^2*VY_p + VE_o
```

### Binary variants: liability-threshold identification

`EquilibriumBinary/` and `DisequilibriumBinary/` mirror the ten continuous models exactly (same
algebra, constraints, and file/function naming with `_Binary_` inserted); only the phenotype
handling changes. Whichever of `Yp1`/`Ym1`/`Yo1` are observed are dichotomized and fit as an OpenMx
**liability-threshold model**: `mxFactor()` on the phenotype columns, an `mxMatrix` of thresholds
passed to `mxExpectationNormal(thresholds=, threshnames=)`, and phenotype means fixed at 0 (the PGS
columns stay continuous with freely estimated means, exactly as in the continuous scripts).

A liability-threshold model can't separately estimate the phenotype's total variance and its
threshold, so `VY` (and `VY_p`/`VY_o` in the DiffTrait models) is **fixed to 1** instead of freely
estimated; the pre-existing `VY_Constraint`-style `mxConstraint`s then force the SEM-implied
variance algebra to equal 1, exactly as they previously forced it to equal a free `VY`. In the
SameTrait models `Yp1`/`Ym1`/`Yo1` share a single free threshold (same trait ⇒ same assumed
prevalence); in the DiffTrait models `Yp1`/`Ym1` (parent trait) share one threshold and `Yo1` (a
different trait) gets its own.

**Identifying the DiffTrait_LatentYPYM models (#3, #8).** With both traits fully latent on the
parent side, the 5-variable data (`Yo1, Tp1, NTp1, Tm1, NTm1`) supply only 4 independent moments
(`VY, thetaT, thetaNT, gt`), one fewer than the model's free dimensions even after both RDR
constraints. Concretely, `gc_eq`/`gt_eq` only pin down the **product** `Omega_p^2*mu` (via
`gc == Omega_p*mu*Omega_p`), not `Omega_p` and `mu` individually — fixing an unrelated parameter
(we first tried `f`) does **not** resolve this: it still converges cleanly, just to a confidently
**wrong** point, not merely a noisy one (verified empirically). Passing an externally-known
`mu_fixed` (e.g. the known/assumed assortative-mating co-path for this trait) breaks exactly that
degeneracy — `Omega_p` becomes uniquely determined by `gc == Omega_p^2*mu_fixed` — and makes the
model point-identified (verified to recover every parameter to ~1e-4 in the continuous case).
Leaving `mu_fixed = NULL` reproduces the original under-identified behavior (a warning is emitted).

Fixing `mu` alone leaves one more gap: nothing in the original scripts constrained the *sign* of
`Omega_p`, `delta_p`, or `f` (only `a_p`/`a_o` had a `lbound = .001`), so once `mu` is fixed the
optimizer can land on a mirror-image solution with `Omega_p`, `delta_p`, `f` (and everything that
depends on them) all wrong. All four `DiffTrait_LatentYPYM` scripts now give `delta_p`, `Omega_p`,
and `f` the same `lbound = .001` that `a_p`/`a_o` already had, which eliminates that spurious
solution (verified: `Omega_p` recovers to ~1e-6 with the bound in place, vs. landing on `-Omega_p`
without it).

## 4. Simulation design (parameter-recovery tests)

Test scripts: `Equilibrium/test_EquilibriumModels.R`, `Disequilibrium/test_DisequilibriumModels.R`,
`EquilibriumBinary/test_BinaryEquilibriumModels.R`, and
`DisequilibriumBinary/test_BinaryDisequilibriumModels.R` (run each with `Rscript`; outputs go to
`test_data/` and `test_results/` inside the folder).

**Step 1 — scalar iterative math.** `Equilibrium/00-UniIterativeMath.R` iterates the full AM+VT
recursion (the scalar translation of `PaperScripts/09-fitWithMVN.R`, mate correlation held fixed at
`am`, so `mu = am/VY` each generation) from a random-mating base until every quantity changes by
< 1e-12 (~37 generations). `Disequilibrium/00-UniIterativeMath_DisEq.R` iterates the VT-only
recursion (`mu = 0`) to its fixed point — checked against the closed forms `Omega = delta*k/(1-f)`
etc. — and then applies **one** generation of AM (`mu = am/VY`, `gt = Omega²mu`, …). Both helpers
assert, after convergence, that the results satisfy the corresponding fitting scripts' constraint
algebra to 1e-8.

For the binary suites, `EquilibriumBinary/00-UniIterativeMath_Binary.R` and
`DisequilibriumBinary/00-UniIterativeMath_DisEq_Binary.R` wrap the same continuous recursions: each
runs the recursion once with arbitrary trial `delta,a,VE` to find the natural (non-1) `VY`, then
re-runs it with `(delta,a,VE) -> (s·delta, s·a, s²·VE)` where `s = 1/sqrt(VY)`, holding `f` and `am`
fixed. This exploits an exact homogeneity of the fixed point under liability rescaling (verified
numerically to machine precision) and lands the fixed point at `VY == 1` exactly, with `mu` becoming
exactly `am`, `Omega/Gamma/w/v/thetaT/thetaNT/Yp_PGSm` scaling by `s`, `Yo_Yp/Yp_Ym` scaling by
`s²`, and the PGS-PGS covariances (`gt/ht/gc/hc/ic`) unchanged. For the DiffTrait models,
`solve_VYo_*_Binary()` then picks offspring `delta_o/a_o` freely and solves `VE_o` in closed form so
the offspring's implied liability variance is *also* exactly 1, using the already-standardized
parent-side quantities — avoiding any cross-generation rescaling of the shared `f` parameter.

**Step 2 — true parameters.** Parental / same-trait values follow the uniModelBias study (trait 1):
`delta = .2` (vg = .64 with 60/64 latent), `a = √.60`, `VE = .36`, `f = .15`, mate correlation
`am = .4`, `k = j = .5`. Offspring trait for DiffTrait conditions: continuous suites use
`delta_o = .3`, `a_o = .65`, `VE_o = .5` for the observed-parent models; binary suites use
`delta_o = .2`, `a_o = .3` chosen directly on the (already-standardized) liability scale. For the
Latent DiffTrait models (#3, #8), which constrain `VY_o = VY_p`, `VE_o` is instead **solved** so the
constraint holds exactly in the generating process (`VY_p − (all other VY_o terms)` continuous;
`solve_VYo_*_Binary(target = 1)` binary).

**Step 3 — RDR heritability inputs.** Each `h2_RDR*` argument is computed from the truth through the
same formula the script's constraint uses, e.g. `h2_RDR = (2a²j+2δ²k)(2a²j+2δ²k+VE) / VY` (with the
generation-specific `a, δ, VE` and whichever `VY` the script's `rdr_right` references), then passed
into the fitting function. The DiffTrait_Latent models additionally receive `mu_fixed` set to the
truth's `mu`.

**Step 4 — expected covariance/thresholds and data.** For each model the 7×7 (or 5×5, latent
conditions) expected covariance matrix is assembled with exactly the entries of that script's
`expCov` algebra, checked for symmetry and positive definiteness (`chol`), and one liability dataset
is drawn with `MASS::mvrnorm(..., empirical = TRUE)` — so the sample covariance equals the expected
matrix exactly. The continuous suites use n = 32,000 and fit the liability directly. The binary
suites use n = 100,000 and dichotomize whichever phenotype columns are observed at fixed thresholds
on the standardized liability scale (`qnorm(0.75)` for SameTrait/parent traits, `qnorm(0.60)` for
the DiffTrait offspring trait — deliberately different from the parent's, to exercise the
shared-vs-separate threshold structure). Data are written as TSV with the column names/order each
script expects and fitted via the script's own function.

**Step 5 — comparison.** Every free parameter's estimate is joined to its true value by OpenMx
label; a model is flagged when any structural parameter misses truth by more than a tolerance
(`.02` for the continuous suites, `.05` for the binary suites — see §5). Per-model tables are saved
as `test_results/*_recovery.csv`.

Because the simulated covariance is generated from the same equations the models fit, these are
**coding + identification** tests; they do not independently validate the equilibrium/disequilibrium
derivations themselves. The continuous suites additionally get an *exact* check: because
`mvrnorm(empirical = TRUE)` makes the sample covariance equal the population target bit-for-bit, a
correctly-coded, identified model must recover the generating parameters to near-machine precision.
The binary suites don't have that luxury — dichotomizing destroys the exact-covariance property —
so recovery there is only as good as any finite-sample ordinal MLE, and is expected to carry real
sampling noise (see §5).

## 5. Results

Continuous suites (`extraTries = 8`, n = 32,000, flag tolerance `.02`):

| Model | Result | Max abs. error (structural) |
|---|---|---|
| Eq SameTrait Observed | **PASS** | 6e-05 |
| Eq SameTrait Latent (RDR) | **PASS** | 5e-04 |
| Eq DiffTrait Latent (2×RDR, `mu_fixed`) | **PASS** | ~1e-4 |
| Eq DiffTrait Observed, estimated a_p | **PASS** | 1e-04 |
| Eq DiffTrait Observed, RDR-fixed a_p | **PASS** | 6e-05 |
| DisEq SameTrait Observed | **PASS** | 4e-05 |
| DisEq SameTrait Latent (RDR) | **PASS** | 5e-05 |
| DisEq DiffTrait Latent (2×RDR, `mu_fixed`) | **PASS** | ~1e-4 |
| DisEq DiffTrait Observed, estimated a_p | **PASS** | 5e-05 |
| DisEq DiffTrait Observed, RDR-fixed a_p | **PASS** | 3e-05 |

All ten continuous models recover every free parameter essentially exactly, which verifies the
OpenMx implementations against the iterative math. (The two DiffTrait_Latent models originally
flagged as under-identified are now identified and passing — see §3's Identification note and §6.)

Binary suites (`extraTries = 15`, `exhaustive = TRUE`, n = 100,000, flag tolerance `.05`; ordinal ML
is noisier than exact-covariance continuous ML, so a tighter tolerance would just flag sampling
noise):

| Model | Result |
|---|---|
| Eq Binary SameTrait Observed | **PASS** |
| Eq Binary SameTrait Latent (RDR) | **PASS** |
| Eq Binary DiffTrait Latent (2×RDR, `mu_fixed`) | FLAG — `VE_parent`/`a_p`/`Gamma_p` within 1.6 SE |
| Eq Binary DiffTrait Observed, estimated a_p | FLAG — `VE_parent`/`a_p` within 1.6 SE |
| Eq Binary DiffTrait Observed, RDR-fixed a_p | **PASS** |
| DisEq Binary SameTrait Observed | **PASS** |
| DisEq Binary SameTrait Latent (RDR) | **PASS** |
| DisEq Binary DiffTrait Latent (2×RDR, `mu_fixed`) | **PASS** |
| DisEq Binary DiffTrait Observed, estimated a_p | FLAG — `VE_parent`/`a_p` within 0.9 SE |
| DisEq Binary DiffTrait Observed, RDR-fixed a_p | **PASS** |

**On the three remaining binary flags: verified as sampling noise, not a coding or identification
bug.** Every flagged parameter misses truth by *less than 1.6 standard errors* — unremarkable for a
single finite (n = 100,000) ordinal-ML draw. We did not take that on faith: a multi-seed check (5
independent datasets) for the `DiffTrait_ObservedYPYM_EstimatedAParent` design showed estimates
scattering on both sides of truth, and a from-scratch check for `DiffTrait_LatentYPYM` showed the
same. That check *also* surfaced two real problems the single-seed run had masked — see §6 items 6
and 7 — both now fixed; the three flags above are what's left over after fixing them, and are
consistent in size with ordinary sampling variation (they shrink as `n` grows, and don't recur with
a fixed sign across seeds). The `EstimatedAParent` variant is consistently noisier than
`FixedAParent` in both the continuous and binary suites (the RDR anchor on the parental `a` removes
a source of estimation noise), so it is the more likely of the two designs to brush up against
whatever flag tolerance is chosen.

## 6. Fixes discovered by the tests

The first test runs surfaced several latent bugs in the fitting scripts, all now fixed:

1. `Equilibrium/UniSEMPGS_DiffTrait_ ObservedYPYM_EstiatedAParent.R` was missing its closing `}` (never parsed).
2. All RDR scripts referenced the `h2_RDR*` function argument directly inside `mxAlgebra`. That resolves under a direct `mxRun` but fails inside `mxTryHard` ("Unknown reference"); the value is now stored in a fixed `mxMatrix` (`h2mat*`).
3. The DiffTrait scripts used identical strings for parameter labels and entity names (e.g. `label="delta_p", name="delta_p"`), which OpenMx rejects; labels now carry an `11` suffix (`delta_p11`, …).
4. The disequilibrium scripts' original starting values implied a non-positive-definite expected covariance matrix (min eigenvalue −0.013), making every NPSOL attempt start infeasible; starts were changed to a PD-verified set.
5. `Equilibrium/UniSEMPGS_SameTrait_LatentYPYM.R` defined but never included `Omega_Constraint`, `gt_constraint`, `gc_constraint`; and the three equilibrium DiffTrait scripts left `w`, `v`, `ic` entirely unconstrained. Both omissions let the optimizer reach alternative solutions; the missing constraints were added.
6. **The DiffTrait_LatentYPYM design (#3, #8) was genuinely under-identified**, not just imprecisely estimated: the 5-variable latent data supply one fewer independent moment than the model's free dimensions. Fixing an externally-known `mu_fixed` (not `f` — fixing `f` alone converges cleanly to a *wrong* point, since the remaining `Omega_p`/`mu` degeneracy is unrelated to `f`) resolves it; see §3.
7. Once `mu_fixed` was added, all four `DiffTrait_LatentYPYM` scripts (continuous + binary, both regimes) turned out to have a second, previously-invisible problem: `Omega_p`, `delta_p`, and `f` had no sign constraint (only `a_p`/`a_o` did), so the optimizer could land on a mirror-image solution with all three (and everything depending on them) wrong. Added `lbound = .001` to `Omega_p`, `delta_p`, and `f` in all four scripts, matching the pre-existing `a_p`/`a_o` convention; verified `Omega_p` now recovers to ~1e-6 instead of landing on `-Omega_p`.
8. The `itlo`/`itol` pair in the disequilibrium scripts (both regimes' `SameTrait`/`DiffTrait` variants, continuous and binary — 10 scripts total) was dead weight: for scalars `Gamma*mu*Omega == Omega*mu*Gamma`, so `itlo` and `itol` always held the same value and were only ever used as `a*itlo + a*itol`. Replaced with a single `ic` quantity (`ic_Algebra = Omega*mu*Gamma`) and `thetaNT = 2*a*ic + 2*delta*gt + w`, matching the equilibrium scripts' existing `ic`-based convention. `itlo`/`itol` only need to be distinguished once `Omega`/`Gamma` become asymmetric (bivariate) matrices.

## 7. Reproducing

```sh
cd Scripts4Dinka
Rscript Equilibrium/test_EquilibriumModels.R                    # ~5 fits, writes Equilibrium/test_results/
Rscript Disequilibrium/test_DisequilibriumModels.R               # ~5 fits, writes Disequilibrium/test_results/
Rscript EquilibriumBinary/test_BinaryEquilibriumModels.R          # ~5 fits, writes EquilibriumBinary/test_results/
Rscript DisequilibriumBinary/test_BinaryDisequilibriumModels.R    # ~5 fits, writes DisequilibriumBinary/test_results/
```

The binary suites take noticeably longer per fit (~3 minutes total each) because `extraTries = 15`
with `exhaustive = TRUE` is needed for reliable convergence on ordinal likelihoods (see §5/§6);
that setting is also recommended for real-data use of the `DiffTrait_ObservedYPYM_EstimatedAParent`
and `DiffTrait_LatentYPYM` binary scripts specifically, which were the ones observed to need it.

Full console output for the latest run is kept in each folder's `test_results/test_run.log`.
