# Errata for `SEM_PGS_Cascade_model.pdf`

This list collects the errors found while implementing the Cascade model and checking it against GeneEvolve.

The equilibrium Cascade algebra is **correct**. This covers Part III §8–12: Ω̃, Γ̃, ζ, τ and Vγ̃, gt = Ω̃²µ, w/v/VF
with the τµ terms, the cross-mate covariances S_A µ S_C, and cov(Yo, Y\*). The one-AM-generation parental and
cross-mate algebra is also correct: Part IV §15–17 and §19. The simulation matches all of them within a few percent
(`Validation/03_results.ipynb`).

The errors concern the **offspring of the first assortative mating** (disequilibrium), the offspring of a change
in assortment strength (Part VII, error 5), and three entries of the Supp. Table II transcription. Part VII's
parental and cross-mate algebra (§34–37) is correct (`Validation/AMChange/`).

## Notation used below

* Y is the phenotype, T and NT are one parent's transmitted and non-transmitted PGS haplotypes, and LT and LNT the
  corresponding latent-score haplotypes. F is the family environment; k and j are the haplotype variances.
* Subscript 1 is the parental generation and subscript 2 their offspring. The parents come from a random-mating
  population, so their own two haplotypes are uncorrelated (g₁ = h₁ = i₁ = 0).
* One assortative mating (copath µ) makes the spouses' haplotypes covary: g₂ = Ω̃₁²µ, h₂ = Γ̃₁²µ and i₂ = Ω̃₁µΓ̃₁.
* The offspring phenotype is Yo = δ(Tp + Tm) + a(LTp + LTm) + Fo + Eo, with Fo = f(Yp + Ym).
* **θ in the papers is a sum over both parents:** θNT = cov(Yo, NTp + NTm). A data column holds one parental
  haplotype, so its covariance with Yo is θ/2.

The offspring's structure follows from direct expansion:

| Quantity | Value in generation 2 |
|---|---|
| var(Tp) | k (no within-haplotype AM covariance yet) |
| cov(Tp, Tm) | g₂ (the offspring's two haplotypes covary only *across* parents) |
| cov(Tp, NTp) | 0 |
| cov(LTm, Tp) | i₂ |
| w₂ = cov(Fo, Tp + Tm) | 2f(Ω₁ + τ₁µΩ̃₁) |
| VF₂ | 2f²(VY₁ + τ₁²µ) |

## Errors

### 1. Part IV §18: generation-2 same-person quantities

The PDF states these are "unchanged in form from Part III Section 11". That form is exact only at equilibrium, where
the within-haplotype AM covariance and the cis covariance both equal gc (so var(T) = k + gc and cov(T, NT) = gc).
After a single assortative mating from a random-mating base, only the cross-parent covariance g₂ exists. The
equilibrium form therefore counts g₂, h₂ and i₂ twice.

| Quantity | PDF §18 | Correct |
|---|---|---|
| Ω₂ = cov(Y₂, own haplotype) | δk + 2δg₂ + 2a·i₂ + ½w₂ | **δk + δg₂ + a·i₂ + ½w₂** |
| Γ₂ | aj + 2ah₂ + 2δi₂ + ½v₂ | **aj + ah₂ + δi₂ + ½v₂** |
| θNT,₂ = cov(Yo, NTp + NTm) | 4δg₂ + 4a·i₂ + w₂ | **2δg₂ + 2a·i₂ + w₂** (per haplotype: δg₂ + a·i₂ + ½w₂) |
| θT,₂ | 2δk + θNT,₂ | 2δk + θNT,₂, with the corrected θNT,₂ (per haplotype: δk + δg₂ + a·i₂ + ½w₂) |
| VY₂ | 2aΓ₂ + 2δΩ₂ + av₂ + δw₂ + VF₂ + Vε with the §18 Ω₂, Γ₂ = 2δ²k + **4**δ²g₂ + 2a²j + **4**a²h₂ + **8**aδi₂ + 2δw₂ + 2av₂ + VF₂ + Vε | **2δ²(k + g₂) + 2a²(j + h₂) + 4aδ·i₂ + 2δw₂ + 2av₂ + VF₂ + Vε** |

The formula for VY₂ is right as written; it is wrong only because it uses the wrong Ω₂ and Γ₂. The §18 expressions
for w₂, v₂ and VF₂ are **correct**.

Simulation check (12 × 10,000-person populations; per-haplotype θ; sim = mean ± SE):

| Mechanism | Quantity | PDF §18 | Corrected | GeneEvolve |
|---|---|---|---|---|
| primary AM | Ω₂ = θT per haplotype | 0.297 | 0.274 | 0.270 ± 0.002 |
| primary AM | θNT per haplotype | 0.095 | 0.072 | 0.067 ± 0.001 |
| primary AM | VY₂ | 1.167 | 1.121 | 1.111 ± 0.007 |
| genetic homogamy | Ω₂ = θT per haplotype | 0.332 | 0.292 | 0.286 ± 0.003 |
| genetic homogamy | θNT per haplotype | 0.130 | 0.090 | 0.088 ± 0.002 |
| genetic homogamy | VY₂ | 1.230 | 1.149 | 1.138 ± 0.006 |
| social homogamy | VY₂ | 1.025 | 1.023 | 1.017 ± 0.003 |

Social homogamy barely moves because g₂ and i₂ are tiny (γ̃ carries almost no PGS).

### 2. Part IV §20: the reduction check is incomplete

§20 says every Part IV equation "collapses exactly onto Part I, Section 3 (Supp. Table II) … term for term". It
checks only g₂, w₂, VF₂, cov(Yp, Ym), cov(Yp, [N]Tm), cov(Y\*, Fo) and cov(Yo, Y\*).

§18's θNT,₂ = 4δg₂ + 4a·i₂ + w₂ does **not** reduce to Supp. Table II's θNT = w₂ + 2δg₂ + 2a·i₂. The Table II
form is the correct one. Ω₂ and Γ₂ do not reduce either.

### 3. Part I §3 (Supp. Table II transcription), AM columns

| Row | Printed | Correct | Reason |
|---|---|---|---|
| Ω₂ (AM, latent) | δk + **2a·i₂** + ½w₂ + δg₂ | δk + δg₂ + **a·i₂** + ½w₂ | cov(LTp, Tp) = 0 and cov(LTm, Tp) = i₂. GeneEvolve, primary AM: printed 0.288, corrected 0.274, simulated 0.270 ± 0.002. |
| Γ₂ (AM, latent) | aj + **2δ·i₂** + ½v₂ + ah₂ | aj + ah₂ + **δ·i₂** + ½v₂ | Same reason. |
| w₁, v₁ (both AM columns) | 2Ω₁µ(VF₁ + δw₁), 2Γ₁µ(VF₁ + δw₁ + av₁) | w₁ = 2fΩ₁, v₁ = 2fΓ₁ | Generation 1 comes from a random-mating base, as in the No-AM columns and PDF Part IV §15. The printed expressions have no f and are shaped like a cross-mate covariance. |
| θNT (AM, no latent) | **2**w₂ + 2δg₂ | w₂ + 2δg₂ | The w₂ coefficient should be 1, as in the AM/latent column. The CORRECTION NOTE in Part I says Table II uses coefficient 1, but this column prints 2. |

The ⁽¹⁾ footnote on VF₂ is right to restore the µ term: the simulation gives VF₂ = 2f²(VY₁ + τ₁²µ).

### 4. Parts II §6 and VI §28–31 (multivariate disequilibrium): the offspring's F carries AM feedback

Both parts use the **parents'** w_p/m = f_pΩ_p + f_mΩ_m and VF = f_pVY_pf_pᵀ + f_mVY_mf_mᵀ in the offspring-side
quantities: θNT, θLNT, VF and VY_o. Part VI's Setup says this lack of feedback is "a property of the source material".
Part II's note explains it by the parents being "the first generation of AM".

That explanation holds for the **parents'** own F. The offspring's F, however, is formed from AM-mated parents, so
its covariance with their haplotypes and its variance carry the cross-mate terms. These are the matrix analogues of
§18's w₂ and VF₂:

```
w_o = f_pΩ_p + f_mΩ_m + f_p τ_p µ Ω̃_m + f_m τ_m µᵀ Ω̃_p
VF_o = f_pVY_pf_pᵀ + f_mVY_mf_mᵀ + f_p τ_p µ τ_mᵀ f_mᵀ + f_m τ_m µᵀ τ_pᵀ f_pᵀ
```

(Part II, without the Cascade: replace τ, Ω̃ by VY, Ω.)

Univariate simulation, primary AM:

| Quantity | Parents' value (as used in Parts II/VI) | Offspring value (corrected) | GeneEvolve |
|---|---|---|---|
| w | 0.071 | 0.100 | 0.098 |
| VF | 0.045 | 0.063 | 0.062 |

Part VI §31's `Ω_p/m,2 = 2a·i_c,2 + 2δ·g_c,2 + δk + ½w` (and hence `VY_p/m,2`) also has the §18 double count. There, gc,2
= ½(gt + gtᵀ) and ic,2 are cross-parent covariances, so their coefficients should be 1, and w should be w_o.

### 5. Parts VII §38 and VIII §46 (two-regime disequilibrium): the inherited cis covariance is dropped

These parts reuse §18's Ω₂, Γ₂, θNT,₂ and VY₂ (VIII §46: Part VI's), with generation-1 parents at a µ₀ equilibrium.
§35 correctly gives those parents the cis covariance gc,₁ = Ω̃₁²µ₀ (and hc,₁, ic,₁). A parent's gamete keeps that
covariance: each locus comes from either of the parent's haplotypes, so var(Tp) = k + gc,₁ and cov(Tp, NTp) = gc,₁.
§38 instead uses the new trans covariance g₂ = Ω̃₁²µ₁ twice and drops gc,₁. Per data column, with g₂, h₂, i₂ from §36:

| Quantity | PDF §38 | Correct |
|---|---|---|
| var(Tp), cov(Tp, NTp) | k, 0 (implicit) | k + gc,₁, gc,₁ |
| θNT per haplotype | 2δg₂ + 2a·i₂ + ½w₂ | **δ(gc,₁ + g₂) + a(ic,₁ + i₂) + ½w₂** |
| θT per haplotype | δk + θNT | δk + θNT, with the corrected θNT |
| VY₂ | 2aΓ₂ + 2δΩ₂ + av₂ + δw₂ + VF₂ + Vε with §18's Ω₂, Γ₂ | **2δ²(k + gc,₁ + g₂) + 2a²(j + hc,₁ + h₂) + 4aδ(ic,₁ + i₂) + 2δw₂ + 2av₂ + VF₂ + Vε** |

The corrected forms reduce to error 1's corrections at µ₀ = 0 and to Part III at µ₀ = µ₁. The printed forms are
exact only at µ₀ = µ₁, where gc,₁ = g₂. §36–37 (g₂, h₂, i₂, the cross-mate covariances), §38's w₂, v₂, VF₂ and
cov(Yo, Y\*) are **correct**. Part VIII §46 has the same omission in matrix form (the offspring's cis input should be
g_c,1 plus the new trans term); it was not implemented or simulated.

Simulation check (`Validation/AMChange/`): 12 GeneEvolve populations per cell; grandparents at a γ̃ correlation of .4
until equilibrium, then the parents mate once at .2, 0 or .6. Per haplotype; sim = mean ± SE:

| Parents' AM | Mechanism | Quantity | PDF §38 | Corrected | GeneEvolve |
|---|---|---|---|---|---|
| weaker (.2) | primary | θNT | 0.085 | 0.101 | 0.098 ± 0.002 |
| weaker (.2) | primary | VY₂ | 0.905 | 0.932 | 0.943 ± 0.008 |
| weaker (.2) | genetic | θNT | 0.116 | 0.143 | 0.141 ± 0.002 |
| weaker (.2) | genetic | VY₂ | 0.872 | 0.916 | 0.920 ± 0.005 |
| stopped (0) | primary | θNT | 0.045 | 0.076 | 0.077 ± 0.003 |
| stopped (0) | genetic | θNT | 0.052 | 0.106 | 0.106 ± 0.003 |
| stopped (0) | genetic | VY₂ | 0.743 | 0.832 | 0.827 ± 0.006 |
| stronger (.6) | primary | θNT | 0.166 | 0.151 | 0.153 ± 0.003 |
| stronger (.6) | genetic | VY₂ | 1.128 | 1.084 | 1.076 ± 0.008 |
| all | primary, genetic | var(Tp) | 0.500 | 0.536, 0.567 | 0.536–0.537, 0.565–0.566 |

All 180 moment comparisons agree with the corrected forms (max |z| 3.2). Under social homogamy gc,₁ ≈ .001, so the
two forms barely differ there.

## Related, but not PDF errors (the original OpenMx DisEq scripts)

The original `OpenMxScripts/Disequilibrium*/` scripts use the papers' two-haplotype θ (θNT = 2a·i + 2δ·gt + w,
θT = 2δk + θNT) for the single-haplotype data columns. That doubles θ, and at δ² = .2 it makes the implied
covariance matrix non-positive-definite.

The same scripts also:
* set the offspring's variance to the parental VY;
* use the parental w and VF;
* omit f·cov(Ym, Yp) from cov(Yo, Yp).

The Cascade DisEq scripts (`OpenMxScripts_Cascade/Disequilibrium*/`) now use the corrected per-haplotype θ with
w_o, the offspring's own VY_o, and the full cov(Yo, Yp).
