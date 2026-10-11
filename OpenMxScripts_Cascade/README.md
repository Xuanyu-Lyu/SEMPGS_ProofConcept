# SEM-PGS Cascade model: OpenMx scripts and iterative math

This folder mirrors `OpenMxScripts/`, but in every script spouses assort on a **latent mating phenotype** γ̃
instead of on the phenotype Y (`SEM_PGS_Cascade_model.pdf`, Parts III and IV):

```
gamma~ = delta~*(T + NT) + a~*(LT + LNT) + 1~*F + 1~*E
```

Primary phenotypic AM, genetic homogamy and social homogamy are special cases. The rest of each model
(VT through Y, the RDR constraints, the binary liability-threshold machinery, identification via `mu_fixed`)
is unchanged from the corresponding original script.

## Layout

| Folder | Iterative math | Fitting scripts (`fitUniSEMPGS_Cascade_*`) |
|---|---|---|
| `Equilibrium/` | `00-UniIterativeMath_Cascade.R` | `UniSEMPGS_Cascade_{SameTrait_ObservedYPYM, SameTrait_LatentYPYM, DiffTrait_LatentYPYM, DiffTrait_ObservedYPYM_EstimatedAParent, DiffTrait_ObservedYPYM_FixedAParent}.R` |
| `Disequilibrium/` | `00-UniIterativeMath_Cascade_DisEq.R` | `UniSEMPGS_Cascade_DisEq_*.R` (same five designs) |
| `EquilibriumBinary/` | `00-UniIterativeMath_Cascade_Binary.R` | `UniSEMPGS_Cascade_Binary_*.R` |
| `DisequilibriumBinary/` | `00-UniIterativeMath_Cascade_DisEq_Binary.R` | `UniSEMPGS_Cascade_DisEq_Binary_*.R` |
| `DisequilibriumAMChange/` | `00-UniIterativeMath_Cascade_AMChange.R` | `UniSEMPGS_Cascade_AMChange_*.R` (the five designs and Model 1) |
| `DisequilibriumAMChangeBinary/` | `00-UniIterativeMath_Cascade_AMChange_Binary.R` | `UniSEMPGS_Cascade_AMChange_Binary_*.R` |
| `Validation/` | the tests: GeneEvolve simulation, fits, and `03_results.ipynb` (`Validation/AMChange/`: the AMChange tests) | |

Function names are the originals with `Cascade_` inserted. For example,
`fitUniSEMPGS_SameTrait_ObservedYPYM` becomes `fitUniSEMPGS_Cascade_SameTrait_ObservedYPYM`. The arguments are
unchanged, plus four new ones for the mating mechanism (next section).

Each of the four folders also holds a Cascade version of Model 1 from Balbona et al. (2021):
`UniSEMPGS_Cascade_Model1.R`, `UniSEMPGS_Cascade_DisEq_Model1.R`, `UniSEMPGS_Cascade_Binary_Model1.R` and
`UniSEMPGS_Cascade_DisEq_Binary_Model1.R`. See [Model 1](#model-1-no-latent-genetic-score).

## The mating multipliers

Following the ETFD Cascade MVN scripts, every path into γ̃ is the corresponding path into Y times a multiplier:

* **`AM_G`** (label `AMGenMulti`) multiplies the genetic paths: δ̃ = AM_G·δ and ã = AM_G·a.
* **`AM_E`** (label `AMEnvMulti`) is the PDF's single 1̃ on the F and E paths. Both paths share one multiplier,
  as in the PDF.

| Mechanism | AM_G | AM_E | γ̃ |
|---|---|---|---|
| primary phenotypic AM | 1 | 1 | Y (every script then reduces exactly to the original) |
| genetic homogamy | 1 | 0 | δ(T+NT) + a(LT+LNT) |
| social homogamy | 0 | 1 | F + E |

The new arguments are `AM_G_value, AM_G_free, AM_E_value, AM_E_free`. γ̃ has no scale of its own: multiplying γ̃
by c leaves every observable unchanged if µ is divided by c². **At least one multiplier must therefore be fixed.**

* **Observed-parent designs** default to `AM_E` fixed at 1, the analogue of the MVN script's fixed `AM_U`, with
  `AM_G` estimated. Genetic homogamy (AM_E = 0) cannot be represented on that scale. To fit it, fix `AM_G = 1`
  and free `AM_E` instead.
* **Latent-parent designs** default to both multipliers fixed. With latent parental phenotypes the mechanism is
  not identified, so it must be supplied as the hypothesized mechanism under test.

µ is the copath between spouses' γ̃: µ = cov(γ̃p, γ̃m) / Vγ̃².

## What changes relative to the original scripts

These are the Cascade quantities from the PDF, as named in the scripts:

| Quantity | Equilibrium (Part III) | Disequilibrium (Part IV) |
|---|---|---|
| `Omega_td` = cov(γ̃, [N]T) | 2δ̃gc + 2ã·ic + δ̃k + ½·1̃·w | δ̃k + ½·1̃·w |
| `Gamma_td` = cov(γ̃, L[N]T) | 2ãhc + 2δ̃ic + ãj + ½·1̃·v | ãj + ½·1̃·v |
| `zeta` = cov(γ̃, F) | δ̃w + ãv + 1̃·VF | same, with VF = 2f²VY |
| `tau` = cov(Y, γ̃) | 2aΓ̃ + 2δΩ̃ + ζ + 1̃·VE | same |
| `Vgamma` = var(γ̃) | 2ãΓ̃ + 2δ̃Ω̃ + 1̃ζ + 1̃²VE | same |
| gt, ht, ic | Ω̃µΩ̃, Γ̃µΓ̃, Ω̃µΓ̃ | same |
| w, v | 2fΩ + 2fτµΩ̃, 2fΓ + 2fτµΓ̃ | parents: unchanged (2fΩ, 2fΓ); offspring: `w_o`, `v_o` = 2fΩ + 2fτµΩ̃, 2fΓ + 2fτµΓ̃ (PDF §18) |
| VF | 2f²VY + 2f²τ²µ (now a free parameter `VF11` held to this by a constraint, because VF and τ depend on each other) | parents: unchanged (2f²VY); offspring: `VF_o` = 2f²VY + 2f²τ²µ (PDF §18) |
| cov(Yp, Ym), cov(Yp, [N]Tm) | τµτ, τµΩ̃ | same |
| cov(Yo, Yp) | δΩ + aΓ + fVY + τµ(δΩ̃ + aΓ̃ + fτ) | same (PDF §19) |
| θNT, θT (per data column) | unchanged: 2δgc + 2a·ic + ½w, δk + θNT | δ·gt + a·ic + ½w_o, δk + θNT (derived; see below) |
| variance of Yo | VY (unchanged) | the offspring's own `VY_o_Algebra` = 2δ²(k+gt) + 2a²(j+ht) + 4δa·ic + 2δw_o + 2av_o + VF_o + VE (derived) |

Everything else is copied from the original scripts: Ω, Γ and VY of the parents, the equilibrium offspring-trait
variance, the RDR constraints, the constraint set each design includes, starting values, and the binary thresholds
and covariates.

### Disequilibrium offspring side, and where the PDF is wrong

The Disequilibrium offspring side is **derived directly** from Yo = δ(Tp + Tm) + a(LTp + LTm) + Fo + Eo. Neither the
original DisEq scripts nor PDF Part IV §18 matched the simulation here.

* **θ per data column.** In the papers θ is the sum over the father's and the mother's haplotype,
  θNT = cov(Yo, NTp + NTm). The data columns (Tp1, NTp1, …) hold one haplotype each, so the scripts' θ is half the
  papers' θ. The scripts say so in a comment next to θ. The first-AM parents' own two haplotypes are uncorrelated,
  so cov(Yo, NTp) = δ·gt + a·ic + ½w_o. The original DisEq scripts used the papers' two-haplotype expression
  instead (twice this), which also made the implied covariance matrix non-positive-definite at δ² = .2.
* **The offspring's F.** The offspring's F comes from AM-mated parents, so w_o, v_o and VF_o carry the τµ terms (PDF
  §18's w₂, v₂, VF₂). The parents' own w and VF have no µ feedback.
* **The offspring's own VY_o.** In both SameTrait scripts (Yo's variance) and the DiffTrait scripts (the offspring
  trait's variance), it replaces the parental VY and the original DiffTrait formula.
* **cov(Yo, Yp)** follows PDF §19, which includes f·cov(Ym, Yp) = fτ²µ. The original DisEq scripts omit this term.

The PDF's own generation-2 formulas in §18 (Ω₂, Γ₂, VY₂, θNT₂) count the AM-induced haplotype covariances twice.
**`PDF_Errata.md`** lists these and the other errors found, with derivations and simulation values.

## Model 1 (no latent genetic score)

The four `*Model1.R` scripts are the Cascade versions of `OpenMxScripts/*/UniSEMPGS_*Model1.R`: Model 1 of
Balbona et al. (2021), in which the PGS explains all heritability. They are the `SameTrait_LatentYPYM` scripts with
`a = 0`, so `a, j, Γ, Γ̃, v, h, i` and the RDR constraint drop out. The data are `Yo1, Tp1, NTp1, Tm1, NTm1`, and γ̃
reduces to δ̃(T + NT) + 1̃·F + 1̃·E:

| Quantity | Equilibrium | Disequilibrium |
|---|---|---|
| `Omega_td` | 2δ̃gc + δ̃k + ½·1̃·w | δ̃k + ½·1̃·w |
| `zeta`, `tau` | δ̃w + 1̃·VF, 2δΩ̃ + ζ + 1̃·VE | same, with VF = 2f²VY |
| gt | Ω̃µΩ̃ (= gc) | Ω̃µΩ̃ (gc = 0) |
| w, VF | 2fΩ + 2fτµΩ̃, 2f²VY + 2f²τ²µ | offspring: `w_o`, `VF_o` of the same form |
| θNT, θT | 2δgc + ½w, δk + θNT | δ·gt + ½w_o, δk + θNT |

As in the latent-parent designs, both multipliers are fixed by default (`AM_G_value = AM_E_value = 1`): with
latent parents the mating mechanism is not identified and must be supplied. The free structural parameters are
then δ, f, µ and VE, which the four moments VY_o, θT, θNT and gt just-identify. Under social homogamy gt is tiny
(Ω̃ = ½w), so µ is weakly identified there, as in the other latent-parent designs (Validation Section 5).

The continuous scripts take their starting values from the sample moments (δ from θT − θNT, Ω from θT, VY from
var(Yo)), as the original Model 1 scripts do; with fixed starts the default search can stop at a far-off point
(`OpenMxScripts/ConceptProof.md` §5).

## Change in assortment strength (AMChange)

`DisequilibriumAMChange/` and `DisequilibriumAMChangeBinary/` implement PDF Part VII: assortative mating was going on
before the parents, but its strength changed in the parents' own mating. Three generations are involved:

* **`_pp`, the grandparents** and every generation before them, mated on γ̃ with copath `mu_pp` long enough to
  reach the Part III equilibrium.
* **`_p`, the parents**, are an ordinary member of that equilibrium population. They mate once, on the same γ̃
  (same AM_G and AM_E), with their own copath `mu_p`.
* **`_o`, the offspring**, are the first generation formed under `mu_p`.

`mu_pp = 0` is the Disequilibrium model and `mu_pp = mu_p` the Equilibrium model. Both are nested, and the scripts
reproduce them exactly (below). The six designs, their data columns and their arguments are those of the other
folders. Function names insert `AMChange_` (and `Binary_`), e.g. `fitUniSEMPGS_Cascade_AMChange_SameTrait_ObservedYPYM`.
The SameTrait_LatentYPYM file defines `fitUniSEMPGS_Cascade_AMChange_SameTrait_LatentParents`, as in the other folders.

| Quantity | Labels / names in the scripts | Model |
|---|---|---|
| Parents (the Eq fixed point at `mu_pp`) | free, held by constraints: `VY_p11`, `VF_p11`, `Omega_p11`, `Gamma_p11`, `w_p11`, `v_p11`, `gc_p11`, `hc_p11`, `ic_p11` | Part III §8–11 at µ_pp. `Omega_td_p`, `Gamma_td_p`, `zeta_p`, `tau_p`, `Vgamma_p` are the parents' shortcuts to γ̃ |
| The grandparents' mating | `mu_pp11`; algebras `gt_pp`, `ht_pp`, `it_pp` | gt_pp = Ω̃_p²µ_pp, held equal to the parents' cis covariance gc_p (and hc_p, ic_p) |
| The parents' mating | `mu_p11`; algebras `gt_p`, `ht_p`, `it_p`, `Yp_Ym`, `Yp_PGSm` | gt_p = Ω̃_p²µ_p, cov(Yp, Ym) = τ_p²µ_p, cov(Yp, [N]Tm) = τ_pµ_pΩ̃_p |
| The offspring's F | `w_o`, `v_o`, `VF_o` | 2f(Ω_p + τ_pµ_pΩ̃_p), 2f(Γ_p + τ_pµ_pΓ̃_p), 2f²(VY_p + τ_p²µ_p) |
| θ per data column | `thetaNT`, `thetaT` | δ(gc_p + gt_p) + a(ic_p + it_p) + ½w_o, and δk + θNT |
| The offspring's variance | `VY_o_Algebra` | 2δ²(k + gc_p + gt_p) + 2a²(j + hc_p + ht_p) + 4δa(ic_p + it_p) + 2δw_o + 2av_o + VF_o + VE |
| cov(Yo, Yp) | `Yo_Yp` | δΩ_p + aΓ_p + fVY_p + τ_pµ_p(δΩ̃_p + aΓ̃_p + fτ_p) (PDF §38) |
| PGS block | | var(Tp) = k + gc_p, cov(Tp, NTp) = gc_p, cov(Tp, Tm) = gt_p |

The offspring rows are derived from Yo = δ(Tp + Tm) + a(LTp + LTm) + Fo + Eo. A parent's transmitted gamete takes each
locus from either of the parent's haplotypes, so it keeps the cis covariance gc_p the parent inherited. PDF §38 drops
this term and uses the new trans term gt_p twice instead, which is exact only when µ_pp = µ_p
(`PDF_Errata.md`, error 5). In the DiffTrait designs, θ, VY_o and `Yo_Yp` use the offspring trait's δ_o, a_o and VE_o.

**What identifies the two copaths.** µ_pp is identified mainly by the parents' cis covariance gc_p = cov(Tp, NTp),
and µ_p by gt_p = cov(Tp, Tm) (and, with observed parents, by cov(Yp, Ym)). Compared with the DisEq scripts the model
adds one parameter (µ_pp) and one moment (gc_p), so the latent designs stay just-identified. Design notes:

* **`DiffTrait_LatentYPYM`** needs one of `mu_p_fixed` and `mu_pp_fixed`. gc_p and gt_p fix only Ω̃_p²µ_pp and
  Ω̃_p²µ_p, so one copath must pin Ω̃_p, as `mu_fixed` does in the DisEq script. Supply `mu_p_fixed` unless the
  parents mate at random (µ_p = 0), where gt_p = 0 says nothing about Ω̃_p; then supply `mu_pp_fixed`.
* **If the parents mate at random (µ_p ≈ 0)**, nothing in the data shows what spouses assorted on, so fix the
  mechanism (AM_G and AM_E) even with observed parents. In the DiffTrait_ObservedYPYM designs a free AM_G then
  trades off against µ_pp at no cost in fit: AM_G = −0.04 with µ_pp = 75 fit the exact data as well as the truth
  (AM_G = 1, µ_pp = .4).
* **Under social homogamy µ_pp is close to unidentified.** The grandparents' assortment leaves gc_p ≈ .001 in the
  genotypes, its only direct trace, so µ_pp has SEs near 1 and wanders by ±.1–.4 between fits that match the truth.
  With a free AM_G (the observed-parent designs) Ω̃_p ≈ ½w can also be cancelled. In one binary
  DiffTrait_ObservedYPYM_EstimatedAParent coding check, AM_G = −0.17 brought Ω̃_p to ≈ .01 and µ_pp ran off to 16
  (NaN SEs). The latent-parent designs are weakly identified under social homogamy, as in the Eq and DisEq scripts.

**Differences from the conventions of the other folders.**

* SameTrait_ObservedYPYM constrains Omega_p and gc_p, which the Eq script, copied from the original, leaves free.
  Here gc_p carries µ_pp.
* The continuous scripts take their starting values from the sample moments (VY_p, Omega_p, gc_p, δ), like the
  Model 1 scripts.
* Data columns other than the design's are ignored, so their order does not matter.

**Iterative math.** `uniIterativeMath_Cascade_AMChange(delta, a, f, VE, am_pp, am_p, AM_G, AM_E)` runs the Part III
recursion at the grandparents' spousal γ̃ correlation `am_pp` to convergence, then the parents' one mating at `am_p`.
µ_pp = am_pp/Vγ̃_p and µ_p = am_p/Vγ̃_p. It returns every `_pp`, `_p` and `_o` quantity, and `$pdf_partVII` (§38 as
printed, for comparison only). `checkCascadeAMChangeConsistency()` asserts the scripts' constraint algebra, and
`diffTraitOffspring_Cascade_AMChange()` gives the DiffTrait offspring trait. The binary file adds
`uniIterativeMath_Cascade_AMChange_Binary()` (VY_p = 1) and `solve_VYo_Cascade_AMChange_Binary()`.

**Validation (`Validation/AMChange/`).** Generating values: BASE and the three mechanisms of `Validation/conditions.R`,
rescaled to VY_p = 1. The grandparents' spousal γ̃ correlation is .4, and the parents' is .2 (`weaker`, the main
condition), 0 (`stopped`) or .6 (`stronger`). GeneEvolve: 12 populations per mechanism, each 30 generations at .4 and
then branched into the three parental matings (~48,700 trios per cell). Results are in
`Validation/AMChange/output/results/summary_*.csv`:

* **Algebra.** At am_pp = 0 the math equals the DisEq math, and at am_pp = am_p the Eq math (max difference 4e-13,
  four mechanisms). Every implied covariance matrix is positive definite.
* **Moments (GeneEvolve vs the math).** All 180 comparisons agree: 20 moments × 3 branches × 3 mechanisms, max
  |z| 3.2. These include the PGS block (var(Tp) = k + gc_p, cov(Tp, NTp), cov(Tp, Tm)), θ, VY_o, cov(Yo, Yp) and
  the offspring's w_o and VF_o. The PDF §38 forms miss θNT, θT and VY_o by |z| 4–20 under primary and genetic AM
  (`PDF_Errata.md`, error 5).
* **Coding check (exact-covariance data, exhaustive search).** 53 of the 54 continuous fits recover the truth within
  .006, most within 4e-4. The exception is AM stopped + social homogamy + DiffTrait_Latent (VE_p off by .024, SE ≈ 1).
  With the scripts' default search, 45 of 54 were within .02: the misses were local optima under social homogamy, and
  the stopped branch before the mechanism was fixed.
  Binary (main condition, n = 50,000, default search): 16 of 18 fits pass (max |z| ≤ 2.6, or within .05). Both
  flags are DiffTrait_ObservedYPYM_EstimatedAParent, which was also the noisiest binary design before
  (`OpenMxScripts/ConceptProof.md` §5). Under genetic homogamy AM_E is .10 off, with NaN SEs for 9 parameters. Under
  social homogamy µ_pp runs off (above). The primary-AM SameTrait_ObservedYPYM fit missed AM_G and both copaths by
  2.1–2.6 SE. On three more datasets (`--seed 11`–`13`) every parameter's mean z over the four fits lies between
  −0.6 and 0.7. The primary-AM DiffTrait_ObservedYPYM_FixedAParent fit was rerun with `--extraTries 2` (max
  error .005, OpenMx status 6) after the default search ran for two hours.
* **Nesting.** On the Eq and DisEq coding-check data, all 30 continuous AMChange fits reach the −2LL of the
  corresponding Eq or DisEq script (within 3e-4). They return µ_pp = µ (Eq) or 0 (DisEq) and recover every other
  parameter within .007. The exception is DiffTrait_Latent under social homogamy on the DisEq data: the −2LL matches,
  but the parental-trait parameters sit .09 away on a flat ridge (NaN SEs).
* **GeneEvolve fits (continuous, 45 fits, exhaustive search).** A likelihood-ratio test never rejects the truth
  (all parameters at their true values vs the fit): p ≥ .09 in every fit. Both copaths are recovered: under primary
  AM, µ_pp = .39–.47 (truth .40) and µ_p = .18–.22 (truth .20) in the main condition, and µ_p = .60–.63 (truth .60)
  in the `stronger` one. The Wald z of single parameters is less trustworthy. SameTrait_ObservedYPYM reports NaN or
  too-small SEs for several parameters (|z| up to 28 for AM_G = 1.017, truth 1). Under social homogamy the estimates
  wander along the weakly identified directions (with p_truth ≥ .14). Use likelihood-based intervals (`mxCI`) or a
  bootstrap with these scripts, as recommended for the DisEq SameTrait_Observed script
  (`OpenMxScripts/ConceptProof.md` §5).
* **GeneEvolve fits (binary, main condition).** Six fits: SameTrait_Latent, DiffTrait_Latent and
  SameTrait_Observed under primary and genetic AM, all with max |z| ≤ 2.4. The binary DiffTrait_ObservedYPYM
  GeneEvolve fits were not run. With ~48,700 trios one optimizer attempt takes 10–25 min, and the default search
  (`extraTries = 30`) had not stopped after an hour. Their algebra is that of the continuous scripts, checked on the
  same data, and their binary coding is checked by the coding check above.

## Iterative math

* `uniIterativeMath_Cascade()` (Eq) iterates PDF Part III from a random-mating, no-VT generation 0 (the same
  starting state as GeneEvolve). It returns the fixed point, all Cascade shortcuts, and the generation-by-
  generation `trajectory`. `checkCascadeEquilibriumConsistency()` asserts the fitting scripts' constraint algebra,
  including that τ's two PDF expressions agree.
* `uniIterativeMath_Cascade_DisEq()` runs the original VT-only recursion, then one mating event on γ̃. It returns
  the offspring generation's `w_o`, `v_o`, `VF_o`, `VY_off` and the per-haplotype θ. It also returns the PDF Part IV
  §18 formulas as printed in `$pdf_partIV`, for comparison only.
* `diffTraitOffspring_Cascade*()` return the DiffTrait offspring-trait quantities. They solve VE_o when
  `VY_o_target` is given.
* The binary files add `uniIterativeMath_Cascade*_Binary()`. These use the homogeneity (δ, a, VE) → (s·δ, s·a, s²·VE),
  under which γ̃ scales with Y and µ by 1/s², to land exactly on VY = 1. They also add `solve_VYo_Cascade_*_Binary()`.
* With AM_G = AM_E = 1, the equilibrium recursion reproduces `OpenMxScripts/Equilibrium/00-UniIterativeMath.R` to
  1e-13. The DisEq recursion reproduces the parental quantities of the original DisEq math. Its offspring
  quantities differ as described above.

## Validation

`Validation/` needs GeneEvolve-Python on the `cascade-mating` branch (the `mating_paths` argument):

```sh
cd OpenMxScripts_Cascade/Validation
Rscript 00_truth.R                                 # iterative-math truth + generating values
python 01_simulate.py --reps 12 --workers 6        # GeneEvolve populations -> output/data, output/sim
python 01_simulate.py --reps 36 --only DisEq       # 24 more DisEq populations (36 in total; see the notebook)
Rscript 02_fit_models.R --mode mvn --workers 4     # exact-covariance coding check (all 60 fits)
Rscript 02_fit_models.R --mode sim --workers 4     # fits to the GeneEvolve trios (all 60 fits)
# the social-homogamy investigation
Rscript 04_social_homogamy.R --workers 10          # per-population fits (3 designs) + profile likelihood of mu
Rscript 05_moment_sensitivity.R                    # which data moments move the SameTrait_Obs estimates
python 06_replicate_eq_social.py                   # 24 more Eq social-homogamy populations (output/replication)
Rscript 06_replicate_eq_social.R                   # ... and their fits
Rscript 07_original_diseq_check.R                  # the corrected ORIGINAL DisEq script on GeneEvolve data
Rscript 08_model1.R --workers 3                    # the Cascade Model 1 scripts (output/model1)
jupyter nbconvert --to notebook --execute --inplace 03_results.ipynb
# the AMChange (change in assortment strength) scripts
cd AMChange
Rscript 00_truth.R                                 # reductions to the Eq/DisEq math, truth, generating values
python 01_simulate.py --reps 12 --workers 3        # GeneEvolve: am_pp = .4 history, then am_p = .2 / 0 / .6
Rscript 02_fit_models.R --mode mvn                 # coding check: 54 continuous (all branches) + 18 binary (main branch)
Rscript 02_fit_models.R --mode mvn --binary 1 --mech primary --design SameTrait_ObservedYPYM --seed 11   # also 12, 13
Rscript 02_fit_models.R --mode reduction           # the AMChange scripts on the Eq/DisEq exact data (needs ../ mvn run)
Rscript 02_fit_models.R --mode sim --binary 0 --branch all     # 45 continuous GeneEvolve fits
Rscript 02_fit_models.R --mode sim --binary 1 --mech primary,genetic \
    --design SameTrait_LatentYPYM,DiffTrait_LatentYPYM,SameTrait_ObservedYPYM   # 6 binary GeneEvolve fits
Rscript 03_summary.R                               # all tables -> output/results/summary_*.csv
```

`03_results.ipynb` holds all test results. Every table and figure has a note on how to read it:

1. Coding check: every script fit to data drawn from its own implied matrix.
2. Whether GeneEvolve reaches equilibrium.
3. Simulated moments vs the iterative math, including the PDF §18 comparison.
4. OpenMx estimates from GeneEvolve data vs the truth.
5. Why some social-homogamy fits are off.
6. The corrected original DisEq scripts on GeneEvolve data.

The notebook ends with a summary of findings. In brief:

* **Equilibrium.**
  * The Cascade math agrees with GeneEvolve for primary AM, genetic homogamy and social homogamy. All moments are
    within |z| ≤ 2.5, mostly 1–2%.
  * The equilibrium models, continuous and binary, are unbiased on GeneEvolve data, except for the social-homogamy
    cases below.
* **Disequilibrium, after the offspring-side fix.**
  * All 66 moment comparisons agree with GeneEvolve over 36 populations: max |z| 3.2, largest relative
    difference 4.4%.
  * All 30 DisEq models pass the coding check.
  * On GeneEvolve data, 14 of 15 continuous and 14 of 15 binary fits have max |z| ≤ 3.1. The exceptions are
    SameTrait-Observed under social homogamy (below), and binary DiffTrait-EstimatedAParent under genetic
    homogamy, where the offspring-trait split is ~2% off (z = 6).
* **Social homogamy (notebook Section 5).** Under social homogamy the mating leaves almost no trace in the
  genotypes: Ω̃ ≈ ½w, so gt ≈ 0.001. Only three designs are affected, each for a different reason:
  * **Latent-parent designs.** µ (SameTrait) and the parent-trait split (DiffTrait) are identified only through gt,
    so in practice they can't be estimated. Exact data still recover the truth, so this is a limit of the design,
    not a coding error.
  * **DisEq SameTrait-Observed.** The estimates are unbiased across 36 populations. OpenMx's reported SEs in this
    script are 2–14× too small, worst for µ under social homogamy. The apparent µ bias (z = 9.5) is about 1.5 real
    SE, which the profile likelihood confirms.
  * **Eq SameTrait-Observed.** A shift along a flat a↓/VE↑/f↑/µ↓ direction, driven by sampling noise in θNT, gc and
    var(T). It shrinks to within noise in 24 replication populations, and the simulated equilibrium matches the
    math within 1%.
* **Original DisEq scripts.** After the same correction, the SameTrait-Observed script is unbiased on the GeneEvolve
  primary-AM populations. Its SEs for a, f, VE and Γ are 2–2.8× too small (`07_original_diseq_check.R`).
* **Model 1 (`08_model1.R`; not in the notebook).** Truth: BASE with all genetic variance in the PGS (δ² = .5,
  a = 0), every fit run with the scripts' default search.
  * All 12 coding-check fits (2 regimes × 3 mechanisms × continuous/binary) recover the truth. The continuous fits
    are within 1e-3, and the largest binary miss is |z| = 1.9 (Eq genetic homogamy). Across 12 more datasets for
    that fit, the mean z of every parameter lies between −0.11 and 0.16. Under social homogamy µ has a large SE
    (.26–.47 continuous), so the binary µ misses of .06 are |z| ≤ 0.3. These two fits exceed the .05 tolerance
    on µ only.
  * With AM_G = AM_E = 1 the Cascade scripts reproduce the original Model 1 scripts: same −2LL, estimates within
    2e-5.
  * On all six GeneEvolve SameTrait_Latent populations, the default search reaches the exhaustive-search optimum
    (−2LL within 1e-4). Model 1 is misspecified there (a > 0), so only convergence is checked.
* **The PDF's errors.** See `PDF_Errata.md`.
* **Run settings.** Continuous fits use `extraTries = 30, exhaustive = TRUE`; without it they can stop at local
  optima. Binary fits (single-threaded ordinal FIML) use the scripts' default search.
