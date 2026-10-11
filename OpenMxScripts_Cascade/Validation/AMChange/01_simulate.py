#!/usr/bin/env python3
"""Step 1: simulate the AMChange conditions (a change in assortment strength) with GeneEvolve.

For every mating mechanism and replicate, one population mates on the latent mating phenotype gamma~ with spousal
correlation am_pp (= .4) for HIST_GENS generations, which brings it to the am_pp equilibrium: these are the
grandparents (_pp). The population is then branched, and in each branch the parents (_p) mate ONCE with their own
spousal correlation am_p (weaker .2, stopped 0, stronger .6; output/truth/generating_params.csv from 00_truth.R).
Each branch exports
    output/data/{branch}_{mechanism}_{dataset}[_binary].tsv   pooled trio data for the OpenMx scripts
    output/sim/moments_final.csv                              moments of the parents' mating, per rep
    output/sim/trajectory.csv                                 the same moments, generation by generation
    output/sim/rep_info.csv                                   trios, seed, run time per rep
The simulation, the moments and the trio export are those of ../01_simulate.py (loaded from there). Model units:
haplotype scores are divided by delta (observed) or a (latent), so k = j = .5 in the founders.

Run from anywhere:  python 01_simulate.py [--reps 12] [--workers 3]
"""
import argparse
import contextlib
import copy
import importlib.util
import io
import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HIST_GENS = 30            # generations at am_pp before the parents' mating (the math converges in 22-33)
REPS = 12
SEED_BASE = 20261007
OUT_DATA = os.path.join(HERE, "output", "data")
OUT_SIM = os.path.join(HERE, "output", "sim")
PARTS = os.path.join(OUT_SIM, "parts")
MECHS = ["primary", "genetic", "social"]
DATASETS = ["SameTrait_Observed", "SameTrait_Latent", "DiffTrait_Observed", "DiffTrait_Latent"]


def load_sim_module():
    spec = importlib.util.spec_from_file_location("simulate", os.path.join(os.path.dirname(HERE), "01_simulate.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.PARTS = PARTS          # save_part() writes here
    return mod


def init_worker():
    global SIM
    SIM = load_sim_module()


def run_rep(args):
    mech, rep, prms = args                       # prms: branch -> generating values of this mechanism
    branches = list(prms)
    tags = [f"{b}_{mech}_rep{rep:03d}" for b in branches]
    if all(os.path.exists(os.path.join(PARTS, f"{t}_info.csv")) for t in tags):
        return f"{mech} rep{rep:03d} (exists)"
    seed = SEED_BASE + 1000 * MECHS.index(mech) + rep
    t0 = time.time()
    rng = np.random.default_rng(seed)
    p0 = prms[branches[0]]                       # the history does not depend on the branch
    for b in branches:
        assert all(abs(prms[b][k] - p0[k]) < 1e-12 for k in ("delta", "a", "VE", "f", "am_pp", "AM_G", "AM_E"))
    am_list = [p0["am_pp"]] * HIST_GENS + [p0["am_p"]]
    sim = SIM.build_sim(p0, SIM.make_cv(p0, rng), HIST_GENS + 1, am_list, seed, SIM.mating_paths(p0))
    sim.setup_recording(p0, last_gen=HIST_GENS)
    with contextlib.redirect_stdout(io.StringIO()):
        for g in range(HIST_GENS):               # the grandparental regime: am_pp every generation
            sim._run_generation(gen_abs_idx=g, am_idx=g, pop_size_target=SIM.POP_SIZE, save=False)
    out = []
    for i, b in enumerate(branches):             # branch: the parents' one mating event at am_p
        br = copy.deepcopy(sim)
        br.prm = prms[b]
        br.am_list[HIST_GENS] = prms[b]["am_p"]
        br.cove_mat = np.diag([prms[b]["VE"], prms[b]["VE_o"]])   # VE_o (DiffTrait offspring trait) is per branch
        np.random.seed(seed * 10 + i)
        with contextlib.redirect_stdout(io.StringIO()):
            br._run_generation(gen_abs_idx=HIST_GENS, am_idx=HIST_GENS, pop_size_target=SIM.POP_SIZE, save=False)
        out.append(SIM.save_part(b, mech, rep, seed, br.traj, br.final, prms[b], rng, t0))
    return ", ".join(out)


def merge(branches):
    parts = sorted(os.listdir(PARTS))
    for b in branches:
        for m in MECHS:
            for name in DATASETS:
                for suffix in ("", "_binary"):
                    files = [f for f in parts if f.startswith(f"{b}_{m}_rep") and f.endswith(f"_{name}{suffix}.tsv")]
                    if files:
                        pd.concat([pd.read_csv(os.path.join(PARTS, f), sep="\t") for f in files], ignore_index=True
                                  ).to_csv(os.path.join(OUT_DATA, f"{b}_{m}_{name}{suffix}.tsv"), sep="\t", index=False)
    traj = pd.concat([pd.read_csv(os.path.join(PARTS, f)) for f in parts if f.endswith("_trajectory.csv")],
                     ignore_index=True).rename(columns={"regime": "branch"})
    traj.to_csv(os.path.join(OUT_SIM, "trajectory.csv"), index=False)
    last = traj.loc[traj.groupby(["branch", "mechanism", "rep"])["gen"].idxmax()]
    last.melt(id_vars=["branch", "mechanism", "rep", "gen"], var_name="quantity", value_name="value").dropna(
        subset=["value"]).to_csv(os.path.join(OUT_SIM, "moments_final.csv"), index=False)
    pd.concat([pd.read_csv(os.path.join(PARTS, f)) for f in parts if f.endswith("_info.csv")], ignore_index=True
              ).rename(columns={"regime": "branch"}).to_csv(os.path.join(OUT_SIM, "rep_info.csv"), index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=REPS)
    ap.add_argument("--workers", type=int, default=3)
    args = ap.parse_args()
    for dname in (OUT_DATA, OUT_SIM, PARTS):
        os.makedirs(dname, exist_ok=True)
    gp = pd.read_csv(os.path.join(HERE, "output", "truth", "generating_params.csv"))
    branches = list(dict.fromkeys(gp["branch"]))
    prm = {(r.mechanism, r.branch): r._asdict() for r in gp.itertuples(index=False)}
    tasks = [(m, rep, {b: prm[(m, b)] for b in branches}) for rep in range(1, args.reps + 1) for m in MECHS]
    sim_mod = load_sim_module()
    print(f"POP_SIZE={sim_mod.POP_SIZE} REPS={args.reps} HIST_GENS={HIST_GENS} N_CV={sim_mod.N_CV} "
          f"branches={branches} workers={args.workers}", flush=True)
    t0 = time.time()
    with Pool(args.workers, initializer=init_worker) as pool:
        for msg in pool.imap_unordered(run_rep, tasks):
            print(f"[{time.time() - t0:6.0f}s] {msg}", flush=True)
    merge(branches)
    print(f"merged -> {OUT_DATA}, {OUT_SIM} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
