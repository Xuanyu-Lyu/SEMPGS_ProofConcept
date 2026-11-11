#!/usr/bin/env python3
"""Step 1: GeneEvolve populations with G x Age in the final generation.

The populations are those of the Cascade validation (OpenMxScripts_Cascade/Validation/01_simulate.py, Eq regime):
same engine, generating values (Validation/output/truth/generating_params.csv), seeds and trio sampling. The u = 0
trio files are therefore byte-identical to Validation's parts Eq_{mech}_rep###_SameTrait_Observed.tsv.

G x Age (the model is described in conditions.R). The parents' adult phenotypes, the assortment and F = f*(Yp + Ym)
are unchanged. In the final generation the offspring's phenotype, measured in childhood, uses the child genetic value
in place of the adult one:
    G_child = r*G_adult + sqrt(1 - r^2)*G_U,   r = 1 - u,
    Yc = Y1 - G_adult + G_child = r*(AO1 + AL1) + sqrt(1 - r^2)*G_U + F1 + MT1 + E1,
with G_U = genotype @ beta_U. beta_U is drawn once per population (from its own RNG stream). It is made orthogonal to
the observed (PGS) and the latent parts of the adult effects in the 2p(1-p)-weighted inner product, then scaled to
base-population variance delta^2 + a^2. In the base population, cor(G_adult, G_child) is then exactly r, and the
child's genetic variance equals the adult's.

Designs (--design):
    independent (default)  every u level has its own REPS populations: u level k (k = 0 for u = 0, ..., 5 for u = .5)
                           uses population reps k*REPS + 1 ... (k + 1)*REPS, and each population exports only that u.
    paired                 every population exports every u level (reps 1 ... REPS).

Outputs (under output/):
    data/Eq_{mech}_rep###_u###.tsv   7-column trio data, Yo1 = the child phenotype at that u
    data/Eq_{mech}_pooled_u###.tsv   the populations of that u pooled
    sim/moments.csv                  offspring moments per mechanism x rep x u (model units; every u is computed)
    sim/diagnostics.csv              checks of beta_U and G_U per mechanism x rep
    sim/rep_info.csv                 trios, seed, run time, and whether the u = 0 file equals Validation's part

Run from anywhere:  python 01_simulate.py [--design independent] [--reps 30] [--workers 3] [--mech primary,genetic,social]
"""
import argparse
import contextlib
import filecmp
import importlib.util
import io
import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
_spec = importlib.util.spec_from_file_location(
    "cascade_validation_sim", os.path.join(REPO, "OpenMxScripts_Cascade", "Validation", "01_simulate.py"))
V = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V)          # Validation's simulation: CascadeValidationSim, make_cv, build_sim, seeds, ...

U_VALUES = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]     # unique share of each generation's genetic variance; = conditions.R
REPS = 30                                      # populations per u level (independent) or in total (paired)
OUT_DATA = os.path.join(HERE, "output", "data")
OUT_SIM = os.path.join(HERE, "output", "sim")
PARTS = os.path.join(OUT_SIM, "parts")
OBS_COLS = ["Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1"]
fcol, cov = V.fcol, V.cov


def utag(u):
    return f"u{round(100 * u):03d}"


def ycol(u):
    return f"Yc_{utag(u)}"


def unique_effects(cv, prm, seed):
    """Child-specific effects beta_U: orthogonal (2p(1-p)-weighted) to the observed and the latent adult effects,
    with base-population variance delta^2 + a^2."""
    rng = np.random.default_rng([seed, 1])     # own stream: the population and the trio sample stay Validation's
    p = cv["maf"].to_numpy()
    het = 2 * p * (1 - p)
    alpha, mask = cv["alpha1"].to_numpy(), cv["mask_obs1"].to_numpy()
    b = rng.standard_normal(len(cv))
    for x in (alpha * mask, alpha * (1 - mask)):   # disjoint supports, so one pass leaves b orthogonal to both
        b = b - np.sum(het * b * x) / np.sum(het * x * x) * x
    return b * np.sqrt((prm["delta"] ** 2 + prm["a"] ** 2) / np.sum(het * b * b))


class GxAgeSim(V.CascadeValidationSim):
    """Validation's simulation; the final generation's offspring also get the child phenotype at every u."""
    beta_u = None

    def _post_offspring_hook(self, offspring_data, gen_abs_idx):
        offspring_data = super()._post_offspring_hook(offspring_data, gen_abs_idx)
        if gen_abs_idx == self.last_gen:
            off = self.final[1]                    # the copy Validation's hook keeps for the trio export
            xo = np.asarray(offspring_data["XO"], dtype=float)
            g_adult = fcol(off, "AO1") + fcol(off, "AL1")
            assert np.allclose(xo @ self.cv_info["alpha1"].to_numpy(), g_adult)
            g_u = xo @ self.beta_u
            y1 = fcol(off, "Y1")
            off["GU1"] = g_u
            for u in U_VALUES:
                r = 1.0 - u
                off[ycol(u)] = y1 - (1 - r) * g_adult + np.sqrt(1 - r ** 2) * g_u
            assert np.array_equal(fcol(off, ycol(0.0)), y1)
        return offspring_data


def parents_of(par, off):
    ids = par.assign(ID=par["ID"].astype(np.int64)).set_index("ID")
    return (fcol(ids.loc[off["Father.ID"].astype(np.int64).to_numpy()], "Y1"),
            fcol(ids.loc[off["Mother.ID"].astype(np.int64).to_numpy()], "Y1"))


def offspring_moments(par, off, prm):
    """VY_o, thetaT, thetaNT (per haplotype) and cov(Yo, Yp) of the child phenotype at every u, in model units and
    defined as in Validation's generation_stats (every offspring with its own parents)."""
    d = prm["delta"]
    yf, ym = parents_of(par, off)
    t_hap = [fcol(off, "TPO1") / d, fcol(off, "TMO1") / d]
    nt_hap = [fcol(off, "NTPO1") / d, fcol(off, "NTMO1") / d]
    rows = []
    for u in U_VALUES:
        yo = fcol(off, ycol(u))
        m = {"VY_o": np.var(yo, ddof=1), "thetaT": np.mean([cov(yo, h) for h in t_hap]),
             "thetaNT": np.mean([cov(yo, h) for h in nt_hap]), "Yo_Yp": np.mean([cov(yo, yf), cov(yo, ym)])}
        rows += [{"u": u, "quantity": k, "value": float(v)} for k, v in m.items()]
    return rows


def diagnostics(cv, beta_u, par, off, prm):
    """beta_U in the base population (orthogonality, variance) and G_U in the final generation's offspring."""
    d = prm["delta"]
    p = cv["maf"].to_numpy()
    het = 2 * p * (1 - p)
    alpha, mask = cv["alpha1"].to_numpy(), cv["mask_obs1"].to_numpy()
    g_adult, g_u = fcol(off, "AO1") + fcol(off, "AL1"), fcol(off, "GU1")
    yf, ym = parents_of(par, off)
    t_hap = [fcol(off, "TPO1") / d, fcol(off, "TMO1") / d]
    nt_hap = [fcol(off, "NTPO1") / d, fcol(off, "NTMO1") / d]
    out = {"VG0": prm["delta"] ** 2 + prm["a"] ** 2,
           "base_dot_obs": float(np.sum(het * beta_u * alpha * mask)),
           "base_dot_lat": float(np.sum(het * beta_u * alpha * (1 - mask))),
           "base_VU": float(np.sum(het * beta_u ** 2)),
           "VU": float(np.var(g_u, ddof=1)), "VG_adult": float(np.var(g_adult, ddof=1)),
           "cor_GU_Gadult": float(np.corrcoef(g_u, g_adult)[0, 1]),
           "cov_GU_T": float(np.mean([cov(g_u, h) for h in t_hap])),
           "cov_GU_NT": float(np.mean([cov(g_u, h) for h in nt_hap])),
           "cov_GU_Yp": float(np.mean([cov(g_u, yf), cov(g_u, ym)])),
           "cov_GU_F": cov(g_u, fcol(off, "F1"))}
    for u in U_VALUES:
        r = 1.0 - u
        out[f"rg_realized_{utag(u)}"] = float(np.corrcoef(g_adult, r * g_adult + np.sqrt(1 - r ** 2) * g_u)[0, 1])
    return out


def export_gxage(par, off, prm, rng):
    """Validation's export_trios (one random offspring per family, centred on this replicate's means) with Yo1 the
    child phenotype at each u. The family sample is the same rng draw, so u = 0 is Validation's SameTrait_Observed."""
    d = prm["delta"]
    ids = par.assign(ID=par["ID"].astype(np.int64)).set_index("ID")
    fam = off.sample(frac=1.0, random_state=int(rng.integers(2 ** 31))).drop_duplicates(["Father.ID", "Mother.ID"])
    fa = ids.loc[fam["Father.ID"].astype(np.int64).to_numpy()]
    mo = ids.loc[fam["Mother.ID"].astype(np.int64).to_numpy()]
    hap_mean = np.mean([fcol(par, "TPO1").mean(), fcol(par, "TMO1").mean()]) / d
    y1_par = fcol(par, "Y1").mean()
    base = pd.DataFrame({
        "Yp1": fcol(fa, "Y1") - y1_par, "Ym1": fcol(mo, "Y1") - y1_par,
        "Tp1": fcol(fam, "TPO1") / d - hap_mean, "NTp1": fcol(fam, "NTPO1") / d - hap_mean,
        "Tm1": fcol(fam, "TMO1") / d - hap_mean, "NTm1": fcol(fam, "NTMO1") / d - hap_mean})
    return {u: base.assign(Yo1=fcol(fam, ycol(u)) - fcol(off, ycol(u)).mean())[OBS_COLS] for u in U_VALUES}


def run_rep(args):
    mech, rep, prm, us = args                  # us: the u levels this population exports
    tag = f"Eq_{mech}_rep{rep:03d}"
    if os.path.exists(os.path.join(PARTS, f"{tag}_info.csv")):
        return f"{tag} (exists)"
    seed = V.SEED_BASE + 1000 * V.MECHS.index(mech) + rep          # Validation's run_eq seed
    t0 = time.time()
    rng = np.random.default_rng(seed)
    cv = V.make_cv(prm, rng)
    sim = V.build_sim(prm, cv, V.EQ_GENS, [prm["am"]] * V.EQ_GENS, seed, V.mating_paths(prm))
    sim.__class__ = GxAgeSim        # build_sim constructs Validation's class; the subclass only overrides the hook
    sim.beta_u = unique_effects(cv, prm, seed)
    sim.setup_recording(prm, last_gen=V.EQ_GENS - 1)
    with contextlib.redirect_stdout(io.StringIO()):
        sim.run_simulation()
    par, off = sim.final

    trios = export_gxage(par, off, prm, rng)
    for u in us:
        trios[u].to_csv(os.path.join(OUT_DATA, f"{tag}_{utag(u)}.tsv"), sep="\t", index=False)
    ref = os.path.join(V.PARTS, f"{tag}_SameTrait_Observed.tsv")
    same = (filecmp.cmp(os.path.join(OUT_DATA, f"{tag}_{utag(0.0)}.tsv"), ref, shallow=False)
            if 0.0 in us and os.path.exists(ref) else None)

    pd.DataFrame(offspring_moments(par, off, prm)).assign(mechanism=mech, rep=rep).to_csv(
        os.path.join(PARTS, f"{tag}_moments.csv"), index=False)
    pd.DataFrame([{"mechanism": mech, "rep": rep, **diagnostics(cv, sim.beta_u, par, off, prm)}]).to_csv(
        os.path.join(PARTS, f"{tag}_diag.csv"), index=False)
    n = len(trios[0.0])
    pd.DataFrame([{"mechanism": mech, "rep": rep, "seed": seed, "n_trios": n, "u0_equals_validation": same,
                   "time_sec": round(time.time() - t0, 1)}]).to_csv(os.path.join(PARTS, f"{tag}_info.csv"), index=False)
    return f"{tag}: n={n} u0 == Validation part: {same} ({time.time() - t0:.0f}s)"


def merge():
    parts = sorted(os.listdir(PARTS))
    for suffix, name in (("_moments.csv", "moments.csv"), ("_diag.csv", "diagnostics.csv"), ("_info.csv", "rep_info.csv")):
        pd.concat([pd.read_csv(os.path.join(PARTS, f)) for f in parts if f.endswith(suffix)], ignore_index=True
                  ).to_csv(os.path.join(OUT_SIM, name), index=False)
    for mech in V.MECHS:
        for u in U_VALUES:
            files = sorted(f for f in os.listdir(OUT_DATA) if f.startswith(f"Eq_{mech}_rep") and f.endswith(f"_{utag(u)}.tsv"))
            if files:
                pd.concat([pd.read_csv(os.path.join(OUT_DATA, f), sep="\t") for f in files], ignore_index=True
                          ).to_csv(os.path.join(OUT_DATA, f"Eq_{mech}_pooled_{utag(u)}.tsv"), sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", choices=["independent", "paired"], default="independent")
    ap.add_argument("--reps", type=int, default=REPS)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--mech", default=",".join(V.MECHS), help="comma-separated subset of primary,genetic,social")
    args = ap.parse_args()
    mechs = args.mech.split(",")
    assert set(mechs) <= set(V.MECHS), mechs
    for dname in (OUT_DATA, OUT_SIM, PARTS):
        os.makedirs(dname, exist_ok=True)
    gp = pd.read_csv(os.path.join(V.HERE, "output", "truth", "generating_params.csv"))
    prm = {r.mechanism: r._asdict() for r in gp.itertuples(index=False) if r.regime == "Eq"}
    if args.design == "paired":
        tasks = [(m, rep, prm[m], U_VALUES) for rep in range(1, args.reps + 1) for m in mechs]
    else:
        tasks = [(m, k * args.reps + i, prm[m], [u]) for k, u in enumerate(U_VALUES)
                 for i in range(1, args.reps + 1) for m in mechs]
    print(f"POP_SIZE={V.POP_SIZE} EQ_GENS={V.EQ_GENS} design={args.design} reps={args.reps} mechanisms={mechs} "
          f"u={U_VALUES} populations={len(tasks)} workers={args.workers}", flush=True)
    t0 = time.time()
    with Pool(args.workers) as pool:
        for msg in pool.imap_unordered(run_rep, tasks):
            print(f"[{time.time() - t0:6.0f}s] {msg}", flush=True)
    merge()
    print(f"merged -> {OUT_DATA}, {OUT_SIM} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
