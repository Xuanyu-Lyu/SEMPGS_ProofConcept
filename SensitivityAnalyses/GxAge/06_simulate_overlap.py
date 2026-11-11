#!/usr/bin/env python3
"""Step 6: the G x Age design with overlapping SNP sets and a child trait expressed in every generation.

Alternative to 01_simulate.py, which adds the child phenotype in the final generation only, as child effects at the
same SNPs (beta_child = r*beta_adult + sqrt(1 - r^2)*beta_U). Here the child expression is GeneEvolve's trait 2, which
the engine computes in every generation, and the genetic correlation comes from shared SNPs (--snps A S C):
    adult-only   A SNPs with an adult effect only (child effect 0)
    shared       S SNPs, a random S of the A + S adult SNPs, with the same effect at both ages
    child-only   C further SNPs with a child effect only
The adult trait is built as in Validation's make_cv on its A + S SNPs (MAF ~ U(.05, .5), effects N(0, 1)): the
Validation's share of them, N_OBS / N_CV = 40%, form the PGS and are rescaled to delta^2, the rest to a^2. The
child-only effects are scaled so that the child's base-population genetic variance is exactly the adult's,
VG0 = delta^2 + a^2. So r_g = var(shared part) / VG0, about S / (A + S), and each age has the unique share
u = 1 - r_g. Everything else is the Validation simulation: spouses assort on the adult mating phenotype (trait 1); the
child trait gets F2 = f*(Y1_father + Y1_mother) and its own residual with variance VE.
The trios are Validation's DiffTrait_Observed export: the parents' adult phenotypes (Y1) and the offspring's child
phenotype (Y2).

Runs:
    08_results_snpsets.ipynb   all three mechanisms, 50 populations each, output_snpsets/; 1,000 SNPs at every r_g level,
                               each level with its own seeds:
                                 r_g = .33  --snps 400 200 400  --seed-offset 210000   (u = .67)
                                 r_g = .67  --snps 250 500 250  --seed-offset 200000   (u = .33; the defaults)
                                 r_g = .87  --snps 115 770 115  --seed-offset 220000   (u = .13)
                                 r_g = 1    --snps 0 1000 0     --seed-offset 230000   (no G x Age)
    07_overlap_check.ipynb     --snps 500 500 500 (u = .5) --mech primary --reps 30 --seed-offset 100000
                               --out output_overlap

Outputs (under --out):
    data/Eq_{mech}_rep###_u###.tsv   7-column trio data, u = 1 - S / (A + S) in percent; data/Eq_{mech}_pooled_u###.tsv
                                     the populations pooled
    sim/trajectory.csv               generation_stats per generation (child moments: VY_off_o, thetaT_o, thetaNT_o, Yo_Yp_o)
    sim/diagnostics.csv              SNP counts, base-population r_g, final-generation cor(G_adult, G_child)
    sim/rep_info.csv                 trios, seed, run time

Run from anywhere:  python 06_simulate_overlap.py [--snps 250 500 250] [--mech primary genetic social] [--reps 50]
                                                  [--seed-offset 200000] [--out output_snpsets] [--workers 3]
"""
import argparse
import contextlib
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
_spec.loader.exec_module(V)          # Validation's simulation: CascadeValidationSim, build_sim, export_trios, ...

fcol = V.fcol


def utag(u):
    return f"u{round(100 * u):03d}"


def unique_share(snps):
    """u = 1 - S / (A + S): the expected share of each age's genetic variance that the other age does not share."""
    a, s, _ = snps
    return 1 - s / (a + s)


def paths(out):
    root = os.path.join(HERE, out)
    return {"data": os.path.join(root, "data"), "sim": os.path.join(root, "sim"),
            "parts": os.path.join(root, "sim", "parts")}


def make_cv_overlap(prm, snps, rng):
    """SNP table with adult (alpha1) and child (alpha2) effects on overlapping SNP sets; returns it and a summary."""
    n_adult_only, n_shared, n_c = snps
    m = n_adult_only + n_shared                                   # adult SNPs
    n_obs = round(m * V.N_OBS / V.N_CV)                           # PGS SNPs, the Validation's share
    maf = rng.uniform(V.MAF_MIN, V.MAF_MAX, m)
    raw = rng.standard_normal(m)
    mask = np.zeros(m, dtype=int)
    mask[rng.choice(m, n_obs, replace=False)] = 1
    het = 2 * maf * (1 - maf)
    obs, lat = raw * mask, raw * (1 - mask)                      # the adult trait, as in Validation's make_cv
    obs = obs * prm["delta"] / np.sqrt(np.sum(het * obs ** 2))
    lat = lat * prm["a"] / np.sqrt(np.sum(het * lat ** 2))
    adult = obs + lat
    shared = np.zeros(m, dtype=int)
    shared[rng.choice(m, n_shared, replace=False)] = 1
    child = adult * shared
    vg0 = prm["delta"] ** 2 + prm["a"] ** 2
    vg_s = float(np.sum(het * child ** 2))
    maf_c = rng.uniform(V.MAF_MIN, V.MAF_MAX, n_c)
    child_only = rng.standard_normal(n_c)
    if n_c > 0:
        child_only *= np.sqrt((vg0 - vg_s) / np.sum(2 * maf_c * (1 - maf_c) * child_only ** 2))
    zeros = np.zeros(n_c)
    cv = pd.DataFrame({"maf": np.concatenate([maf, maf_c]),
                       "alpha1": np.concatenate([adult, zeros]),
                       "alpha2": np.concatenate([child, child_only]),
                       "mask_obs1": np.concatenate([mask, zeros]).astype(int),
                       "mask_obs2": np.concatenate([mask * shared, zeros]).astype(int)})
    het_all = 2 * cv.maf.to_numpy() * (1 - cv.maf.to_numpy())
    arch = {"n_adult": m, "n_pgs": n_obs, "n_shared": n_shared, "n_adult_only": n_adult_only, "n_child_only": n_c,
            "VG0": vg0, "VG_adult_base": float(np.sum(het_all * cv.alpha1.to_numpy() ** 2)),
            "VG_child_base": float(np.sum(het_all * cv.alpha2.to_numpy() ** 2)), "rg_base": vg_s / vg0}
    return cv, arch


def run_rep(args):
    mech, rep, snps, prm, seed_offset, out = args
    p = paths(out)
    u = unique_share(snps)
    tag = f"Eq_{mech}_rep{rep:03d}"
    if os.path.exists(os.path.join(p["parts"], f"{tag}_{utag(u)}_info.csv")):
        return f"{tag} {utag(u)} (exists)"
    seed = V.SEED_BASE + 1000 * V.MECHS.index(mech) + rep + seed_offset
    prm = dict(prm, VE_o=prm["VE"], delta_o=prm["delta"], a_o=prm["a"])   # trait 2 = the child expression
    t0 = time.time()
    rng = np.random.default_rng(seed)
    cv, arch = make_cv_overlap(prm, snps, rng)
    sim = V.build_sim(prm, cv, V.EQ_GENS, [prm["am"]] * V.EQ_GENS, seed, V.mating_paths(prm))
    sim.setup_recording(prm, last_gen=V.EQ_GENS - 1)
    with contextlib.redirect_stdout(io.StringIO()):
        sim.run_simulation()
    par, off = sim.final
    trios = V.export_trios(par, off, prm, rng)["DiffTrait_Observed"]
    trios.to_csv(os.path.join(p["data"], f"{tag}_{utag(u)}.tsv"), sep="\t", index=False)
    g_adult = fcol(off, "AO1") + fcol(off, "AL1")
    g_child = fcol(off, "AO2") + fcol(off, "AL2")
    diag = {"mechanism": mech, "rep": rep, "u": u, **arch,
            "rg_final_generation": float(np.corrcoef(g_adult, g_child)[0, 1]),
            "VG_child_final": float(np.var(g_child, ddof=1)), "VG_adult_final": float(np.var(g_adult, ddof=1))}
    pd.DataFrame(sim.traj).assign(mechanism=mech, rep=rep, u=u).to_csv(
        os.path.join(p["parts"], f"{tag}_{utag(u)}_trajectory.csv"), index=False)
    pd.DataFrame([diag]).to_csv(os.path.join(p["parts"], f"{tag}_{utag(u)}_diag.csv"), index=False)
    pd.DataFrame([{"mechanism": mech, "rep": rep, "u": u, "seed": seed, "n_trios": len(trios),
                   "time_sec": round(time.time() - t0, 1)}]).to_csv(
        os.path.join(p["parts"], f"{tag}_{utag(u)}_info.csv"), index=False)
    return f"{tag} {utag(u)}: n={len(trios)} r_g(base)={arch['rg_base']:.3f} ({time.time() - t0:.0f}s)"


def merge(out):
    p = paths(out)
    parts = sorted(os.listdir(p["parts"]))
    for suffix, name in (("_trajectory.csv", "trajectory.csv"), ("_diag.csv", "diagnostics.csv"), ("_info.csv", "rep_info.csv")):
        pd.concat([pd.read_csv(os.path.join(p["parts"], f)) for f in parts if f.endswith(suffix)], ignore_index=True
                  ).to_csv(os.path.join(p["sim"], name), index=False)
    files = [f for f in os.listdir(p["data"]) if "_rep" in f]
    for key in sorted({(f.split("_rep")[0], f.rsplit("_", 1)[1]) for f in files}):
        sel = sorted(f for f in files if f.startswith(key[0] + "_rep") and f.endswith("_" + key[1]))
        pd.concat([pd.read_csv(os.path.join(p["data"], f), sep="\t") for f in sel], ignore_index=True
                  ).to_csv(os.path.join(p["data"], f"{key[0]}_pooled_{key[1]}"), sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snps", type=int, nargs=3, default=[250, 500, 250], metavar=("ADULT_ONLY", "SHARED", "CHILD_ONLY"))
    ap.add_argument("--mech", nargs="+", default=V.MECHS, choices=V.MECHS)
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--seed-offset", type=int, default=200000)   # new populations: seeds differ from the other runs'
    ap.add_argument("--out", default="output_snpsets")
    ap.add_argument("--workers", type=int, default=3)
    args = ap.parse_args()
    for dname in paths(args.out).values():
        os.makedirs(dname, exist_ok=True)
    gp = pd.read_csv(os.path.join(V.HERE, "output", "truth", "generating_params.csv"))
    prms = {r.mechanism: r._asdict() for r in gp.itertuples(index=False) if r.regime == "Eq"}
    tasks = [(mech, rep, tuple(args.snps), prms[mech], args.seed_offset, args.out)
             for mech in args.mech for rep in range(1, args.reps + 1)]
    print(f"overlapping-SNP design: SNPs adult-only/shared/child-only = {'/'.join(map(str, args.snps))} "
          f"(u = {unique_share(args.snps):.3f}) mechanisms={','.join(args.mech)} populations={len(tasks)} "
          f"POP_SIZE={V.POP_SIZE} EQ_GENS={V.EQ_GENS} workers={args.workers} -> {args.out}", flush=True)
    t0 = time.time()
    with Pool(args.workers) as pool:
        for msg in pool.imap_unordered(run_rep, tasks):
            print(f"[{time.time() - t0:6.0f}s] {msg}", flush=True)
    merge(args.out)
    print(f"merged -> {args.out} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
