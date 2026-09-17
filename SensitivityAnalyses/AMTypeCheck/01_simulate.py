#!/usr/bin/env python3
"""Simulate replicate populations (trait 1: genetic AM, trait 2: social AM) and export one
SEM-PGS trio dataset per trait and replicate, plus the generating / empirical truth values.

Run from anywhere:
    python 01_simulate.py --dry-run     # print the conditions and implied quantities, no simulation
    python 01_simulate.py               # all replicates in sim_conditions.py (skips ones already on disk)
    python 01_simulate.py --reps 3      # first 3 replicates only

Outputs (all under output/):
    data/rep###_geneticAM.tsv, data/rep###_socialAM.tsv   7-column trio data for the OpenMx script
    truth/rep###_truth.csv                                 true value of every model parameter, per condition
    truth/rep###_summary.csv                               diagnostics (n trios, VY, h2, spousal correlations)
"""
import argparse
import contextlib
import io
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "GeneEvolve-Python"))
from core_simulation import AssortativeMatingSimulation  # noqa: E402
import sim_conditions as C  # noqa: E402

OUT_DATA = os.path.join(HERE, "output", "data")
OUT_TRUTH = os.path.join(HERE, "output", "truth")
CONDITION_NAME = {"genotypic": "geneticAM", "social": "socialAM", "phenotypic": "phenotypicAM"}


class TrioExportSimulation(AssortativeMatingSimulation):
    """Keeps the last parental generation so trios can be built without storing the full history."""
    parent_phen_df = None

    def _post_offspring_hook(self, offspring_data, gen_abs_idx):
        if gen_abs_idx == self.n_burn_in + self.num_generations - 1:
            self.parent_phen_df = self.phen_df.copy()
        return offspring_data


def model_matrices():
    d = np.sqrt(C.PROP_OBS)
    a = np.sqrt(1.0 - C.PROP_OBS)
    return dict(
        cove_mat=np.diag([C.VE, C.VE]),
        f_mat=np.diag([C.F, C.F]),
        s_mat=np.zeros((2, 2)),
        a_mat=np.diag([a, a]),
        d_mat=np.diag([d, d]),
        covy_mat=np.eye(2),
        k2_matrix=np.zeros((2, 2)),
    )


def base_variances(cv_info):
    """Generation-0 (HWE, no LD) variance of the observed and latent genetic value, per trait."""
    p = cv_info["maf"].values
    w = 2 * p * (1 - p)
    out = {}
    for i in (1, 2):
        alpha = cv_info[f"alpha{i}"].values
        mask = cv_info[f"mask_obs{i}"].values
        out[i] = (float(np.sum(w * (alpha * mask) ** 2)), float(np.sum(w * (alpha * (1 - mask)) ** 2)))
    return out


def fcol(df, name):
    return df[name].values.astype(float)


def cov(x, y):
    return float(np.cov(x, y)[0, 1])


def build_trios_and_truth(off, par, i, vao, val, rng):
    """One offspring per family; PGS haplotypes rescaled so that k = j = .5 at generation 0."""
    s_o, s_l = np.sqrt(vao), np.sqrt(val)
    par = par.copy()
    par["ID"] = par["ID"].astype(np.int64)
    par = par.set_index("ID")
    off = off.copy()
    off["Father.ID"] = off["Father.ID"].astype(np.int64)
    off["Mother.ID"] = off["Mother.ID"].astype(np.int64)
    fam = off.sample(frac=1.0, random_state=int(rng.integers(2**31))).drop_duplicates(["Father.ID", "Mother.ID"])
    fa = par.loc[fam["Father.ID"].values]
    mo = par.loc[fam["Mother.ID"].values]

    trio = pd.DataFrame({
        "Yp1": fcol(fa, f"Y{i}"), "Ym1": fcol(mo, f"Y{i}"), "Yo1": fcol(fam, f"Y{i}"),
        "Tp1": fcol(fam, f"TPO{i}") / s_o, "NTp1": fcol(fam, f"NTPO{i}") / s_o,
        "Tm1": fcol(fam, f"TMO{i}") / s_o, "NTm1": fcol(fam, f"NTMO{i}") / s_o,
    })

    # Empirical population truth on the parental generation, in the model's units.
    # Within-person quantities use each parent's own two haplotypes (TPO/TMO = observed, TPL/TML = latent);
    # spouse-level quantities use the (father, mother) pairs behind the exported trios.
    Y, Fv = fcol(par, f"Y{i}"), fcol(par, f"F{i}")
    To = [fcol(par, f"TPO{i}") / s_o, fcol(par, f"TMO{i}") / s_o]
    Tl = [fcol(par, f"TPL{i}") / s_l, fcol(par, f"TML{i}") / s_l]
    Yp, Ym = fcol(fa, f"Y{i}"), fcol(mo, f"Y{i}")
    Tp_o = [fcol(fa, f"TPO{i}") / s_o, fcol(fa, f"TMO{i}") / s_o]
    Tm_o = [fcol(mo, f"TPO{i}") / s_o, fcol(mo, f"TMO{i}") / s_o]
    Tp_l = [fcol(fa, f"TPL{i}") / s_l, fcol(fa, f"TML{i}") / s_l]
    Tm_l = [fcol(mo, f"TPL{i}") / s_l, fcol(mo, f"TML{i}") / s_l]
    VY = float(np.var(Y, ddof=1))

    truth = {
        "delta11": s_o, "a11": s_l, "f11": C.F, "VE11": C.VE,
        "VY11": VY,
        "mu11": cov(Yp, Ym) / VY**2,
        "Omega11": np.mean([cov(Y, h) for h in To]),
        "Gamma11": np.mean([cov(Y, h) for h in Tl]),
        "gc11": cov(To[0], To[1]),
        "hc11": cov(Tl[0], Tl[1]),
        "ic11": np.mean([cov(To[0], Tl[1]), cov(To[1], Tl[0])]),
        "gt11": np.mean([cov(a, b) for a in Tp_o for b in Tm_o]),
        "ht11": np.mean([cov(a, b) for a in Tp_l for b in Tm_l]),
        "itlo11": np.mean([cov(a, b) for a in Tp_l for b in Tm_o]),
        "itol11": np.mean([cov(a, b) for a in Tp_o for b in Tm_l]),
        "w11": 2 * np.mean([cov(Fv, h) for h in To]),
        "v11": 2 * np.mean([cov(Fv, h) for h in Tl]),
        "VF_sim": float(np.var(Fv, ddof=1)),
    }

    A_p = fcol(fa, f"AO{i}") + fcol(fa, f"AL{i}")
    A_m = fcol(mo, f"AO{i}") + fcol(mo, f"AL{i}")
    S_p = fcol(fa, f"F{i}") + fcol(fa, f"E{i}")
    S_m = fcol(mo, f"F{i}") + fcol(mo, f"E{i}")
    A_all = fcol(par, f"AO{i}") + fcol(par, f"AL{i}")
    summary = {
        "n_trios": len(trio),
        "VY": VY,
        "h2": float(np.var(A_all, ddof=1)) / VY,
        "r_spouse_Y": float(np.corrcoef(Yp, Ym)[0, 1]),
        "r_spouse_A": float(np.corrcoef(A_p, A_m)[0, 1]),
        "r_spouse_FE": float(np.corrcoef(S_p, S_m)[0, 1]),
        "VAO_base": vao, "VAL_base": val,
    }
    return trio, truth, summary


def run_rep(r):
    seed = C.SEED_BASE + r
    np.random.seed(seed)
    cv = AssortativeMatingSimulation.prepare_CV_random_selection(
        C.N_CV, 0.0, C.MAF_MIN, C.MAF_MAX, C.PROP_OBS, C.PROP_OBS)
    cv[["alpha1", "alpha2"]] *= np.sqrt(C.VG)   # engine scales effects to VA ~ 1; rescale to VG
    am_mat = np.diag([C.AM[C.TRAITS[1]], C.AM[C.TRAITS[2]]])
    with contextlib.redirect_stdout(io.StringIO()):
        sim = TrioExportSimulation(
            cv_info=cv,
            num_generations=C.N_GENERATIONS, pop_size=C.POP_SIZE,
            mating_type=(C.TRAITS[1], C.TRAITS[2]), mate_on_trait=None,
            am_list=[am_mat] * C.N_GENERATIONS,
            seed=seed, save_each_gen=False, save_covs=False,
            **model_matrices(),
        )
        res = sim.run_simulation()

    bases = base_variances(cv)
    rng = np.random.default_rng(seed)
    truth_rows, summary_rows = [], []
    for i in (1, 2):
        cond = CONDITION_NAME[C.TRAITS[i]]
        trio, truth, summary = build_trios_and_truth(res["PHEN"], sim.parent_phen_df, i, *bases[i], rng)
        trio.to_csv(os.path.join(OUT_DATA, f"rep{r:03d}_{cond}.tsv"), sep="\t", index=False)
        truth_rows += [{"rep": r, "condition": cond, "parameter": k, "true": v} for k, v in truth.items()]
        summary_rows.append({"rep": r, "condition": cond, "seed": seed, **summary})
    pd.DataFrame(truth_rows).to_csv(os.path.join(OUT_TRUTH, f"rep{r:03d}_truth.csv"), index=False)
    pd.DataFrame(summary_rows).to_csv(os.path.join(OUT_TRUTH, f"rep{r:03d}_summary.csv"), index=False)
    return summary_rows


def print_conditions():
    expected_trios = C.POP_SIZE / 2 * (1 - np.exp(-2))
    print("=== Simulation conditions (sim_conditions.py) ===")
    for k in ["N_REPS", "SEED_BASE", "POP_SIZE", "N_GENERATIONS", "N_CV", "MAF_MIN", "MAF_MAX",
              "VG", "PROP_OBS", "VE", "F", "R_SPOUSE_Y", "TRAITS"]:
        print(f"  {k:<14} = {getattr(C, k)}")
    print("\n=== Implied generating values (model scale, k = j = .5) ===")
    print(f"  delta            = sqrt(VG * PROP_OBS)       = {np.sqrt(C.VG * C.PROP_OBS):.4f}")
    print(f"  a                = sqrt(VG * (1 - PROP_OBS)) = {np.sqrt(C.VG * (1 - C.PROP_OBS)):.4f}")
    print(f"  VE               =                              {C.VE:.4f}")
    print(f"  f                =                              {C.F:.4f}")
    print(f"  gen-0 VY         = VG + VE                    = {C.VG + C.VE:.4f}   (grows with AM + VT)")
    print(f"  gen-0 h2         = VG / (VG + VE)             = {C.VG / (C.VG + C.VE):.4f}")
    print(f"  observed CVs     = int(N_CV * PROP_OBS)       = {int(C.N_CV * C.PROP_OBS)} of {C.N_CV}")
    print(f"  trios / dataset  ~ POP_SIZE/2 * (1 - e^-2)    = {expected_trios:.0f}")
    print(f"\n=== Assortment: mating-variable correlations calibrated (00_calibrate_am.py) to "
          f"r(Y_p, Y_m) ~ {C.R_SPOUSE_Y:.2f} at generation {C.N_GENERATIONS} ===")
    mating_var = {"genotypic": "A = AO + AL", "social": "F + E", "phenotypic": "Y"}
    for i in (1, 2):
        t = C.TRAITS[i]
        print(f"  trait {i} ({t:>10}): AM = {C.AM[t]:.3f} on {mating_var[t]}")
    print("  Realised spousal correlations, VY and h2 are recorded per replicate in output/truth/rep###_summary.csv.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=C.N_REPS)
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true", help="regenerate replicates that already exist on disk")
    args = ap.parse_args()

    print_conditions()
    if args.dry_run:
        return
    os.makedirs(OUT_DATA, exist_ok=True)
    os.makedirs(OUT_TRUTH, exist_ok=True)
    conds = [CONDITION_NAME[C.TRAITS[i]] for i in (1, 2)]
    print()
    for r in range(args.start, args.reps + 1):
        if not args.force and all(os.path.exists(os.path.join(OUT_DATA, f"rep{r:03d}_{c}.tsv")) for c in conds):
            print(f"rep {r:03d}: exists, skipped")
            continue
        t0 = time.time()
        rows = run_rep(r)
        msg = "  ".join(f"{s['condition']}: n={s['n_trios']} VY={s['VY']:.3f} h2={s['h2']:.3f} rY={s['r_spouse_Y']:.3f}"
                        for s in rows)
        print(f"rep {r:03d}: {msg}  ({time.time() - t0:.1f}s)")


if __name__ == "__main__":
    main()
