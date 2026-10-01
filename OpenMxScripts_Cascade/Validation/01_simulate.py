#!/usr/bin/env python3
"""Step 1: simulate the SEM-PGS Cascade conditions with GeneEvolve (latent mating phenotype, `mating_paths`).

For every regime x mating mechanism (generating values from output/truth/generating_params.csv, written by
00_truth.R) this runs REPS independent populations and exports
    output/data/{regime}_{mechanism}_{dataset}[_binary].tsv   pooled trio data for the OpenMx scripts
    output/sim/moments_final.csv                              population moments of the final generation, per rep
    output/sim/trajectory.csv                                 the same moments, generation by generation
    output/sim/rep_info.csv                                   trios, seed, run time per rep

Trait 1 is the parent / same trait: spouses assort on its latent mating phenotype
    gamma~ = AM_G*(AO1 + AL1) + AM_E*(F1 + E1)
and it is vertically transmitted (F = f*(Y1_father + Y1_mother)). Trait 2 is the DiffTrait offspring trait:
it loads delta_o / a_o on the same CVs (same observed/latent split) and receives the same F.
    Eq     assortment on gamma~ every generation for EQ_GENS generations (the iterative math converges in ~12)
    DisEq  random mating with VT for DISEQ_BURN generations, then ONE generation of assortment on gamma~; the
           random-mating history is shared, and the population is branched into the three mechanisms.
Model units: haplotype scores are divided by delta (observed) or a (latent), so k = j = .5 in the founders.

Run from anywhere:  python 01_simulate.py [--reps 12] [--workers 6] [--only Eq|DisEq]
(The validation used 12 Eq and 36 DisEq populations: python 01_simulate.py --reps 36 --only DisEq.)
"""
import argparse
import contextlib
import copy
import io
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "GeneEvolve-Python"))
from core_simulation import AssortativeMatingSimulation  # noqa: E402

POP_SIZE = 10000          # individuals per generation (~0.43 * POP_SIZE independent trios per population)
REPS = 12                 # independent populations per condition; trios are pooled across them for fitting
EQ_GENS = 30              # generations of AM + VT for the equilibrium regime
DISEQ_BURN = 20           # generations of random mating + VT before the one AM generation (DisEq)
N_CV = 1000               # causal variants; the observed PGS uses N_OBS of them
N_OBS = 400               # = N_CV * delta^2 / (delta^2 + a^2) for the base values (.2 / .5)
MAF_MIN, MAF_MAX = 0.05, 0.50
SEED_BASE = 20260928
PREV = {"parent": 0.25, "offspring_diff": 0.40}     # binary prevalences; must match conditions.R

OUT_DATA = os.path.join(HERE, "output", "data")
OUT_SIM = os.path.join(HERE, "output", "sim")
PARTS = os.path.join(OUT_SIM, "parts")
MECHS = ["primary", "genetic", "social"]


def fcol(df, name):
    return df[name].to_numpy(dtype=float)


def cov(x, y):
    return float(np.cov(x, y)[0, 1])


def norm_ppf(p):
    from statistics import NormalDist
    return NormalDist().inv_cdf(p)


THRESH = {"parent": norm_ppf(1 - PREV["parent"]), "offspring_diff": norm_ppf(1 - PREV["offspring_diff"])}


class CascadeValidationSim(AssortativeMatingSimulation):
    """Records population moments at every mating and keeps the last parents/offspring for trio export."""

    def setup_recording(self, prm, last_gen):
        self.prm = prm
        self.last_gen = last_gen
        self.traj = []
        self.final = None

    def _post_offspring_hook(self, offspring_data, gen_abs_idx):
        par, off = self.phen_df, offspring_data["PHEN"]
        self.traj.append({"gen": gen_abs_idx, **generation_stats(par, off, self.prm, self._mating_weights[0])})
        if gen_abs_idx == self.last_gen:
            self.final = (par.copy(), off.copy())
        return offspring_data


def generation_stats(par, off, prm, weights):
    """Moments of one mating event in model units: parents (within person), spouses, offspring."""
    d, a = prm["delta"], prm["a"]
    wG = weights.get("AO", 0.0)
    wE = weights.get("F", 0.0)
    s = {}
    Y, F, E = fcol(par, "Y1"), fcol(par, "F1"), fcol(par, "E1")
    A = fcol(par, "AO1") + fcol(par, "AL1")
    g = wG * A + wE * (F + E)
    s.update(VY=np.var(Y, ddof=1), VF=np.var(F, ddof=1), Vgamma=np.var(g, ddof=1),
             zeta=cov(g, F), tau=cov(Y, g), VA=np.var(A, ddof=1))
    have_haps = not par["TPO1"].isna().any()
    if have_haps:
        ho = [fcol(par, "TPO1") / d, fcol(par, "TMO1") / d]
        hl = [fcol(par, "TPL1") / a, fcol(par, "TML1") / a]
        s.update(Omega=np.mean([cov(Y, h) for h in ho]), Gamma=np.mean([cov(Y, h) for h in hl]),
                 Omega_td=np.mean([cov(g, h) for h in ho]), Gamma_td=np.mean([cov(g, h) for h in hl]),
                 gc=cov(*ho), hc=cov(*hl), ic=np.mean([cov(ho[0], hl[1]), cov(ho[1], hl[0])]),
                 hapvar_obs=np.mean([np.var(h, ddof=1) for h in ho]),
                 hapvar_lat=np.mean([np.var(h, ddof=1) for h in hl]),
                 w=2 * np.mean([cov(F, h) for h in ho]), v=2 * np.mean([cov(F, h) for h in hl]))

    # spouses: the (father, mother) pairs that produced offspring
    ids = par.assign(ID=par["ID"].astype(np.int64)).set_index("ID")
    pairs = off[["Father.ID", "Mother.ID"]].astype(np.int64).drop_duplicates()
    fa, mo = ids.loc[pairs["Father.ID"].to_numpy()], ids.loc[pairs["Mother.ID"].to_numpy()]
    gf = wG * (fcol(fa, "AO1") + fcol(fa, "AL1")) + wE * (fcol(fa, "F1") + fcol(fa, "E1"))
    gm = wG * (fcol(mo, "AO1") + fcol(mo, "AL1")) + wE * (fcol(mo, "F1") + fcol(mo, "E1"))
    Yf, Ym, Ff, Fm = fcol(fa, "Y1"), fcol(mo, "Y1"), fcol(fa, "F1"), fcol(mo, "F1")
    s.update(n_pairs=len(pairs), r_gamma=float(np.corrcoef(gf, gm)[0, 1]),
             mu=cov(gf, gm) / s["Vgamma"] ** 2, Yp_Ym=cov(Yf, Ym),
             Yp_Fm=np.mean([cov(Yf, Fm), cov(Ym, Ff)]), Fp_Fm=cov(Ff, Fm))
    if have_haps:
        hof = [fcol(fa, "TPO1") / d, fcol(fa, "TMO1") / d]
        hom = [fcol(mo, "TPO1") / d, fcol(mo, "TMO1") / d]
        hlf = [fcol(fa, "TPL1") / a, fcol(fa, "TML1") / a]
        hlm = [fcol(mo, "TPL1") / a, fcol(mo, "TML1") / a]
        s.update(gt=np.mean([cov(x, y) for x in hof for y in hom]),
                 ht=np.mean([cov(x, y) for x in hlf for y in hlm]),
                 it=np.mean([cov(x, y) for x in hof for y in hlm] + [cov(x, y) for x in hlf for y in hom]),
                 Yp_PGSm=np.mean([cov(Yf, h) for h in hom] + [cov(Ym, h) for h in hof]),
                 Yp_LGSm=np.mean([cov(Yf, h) for h in hlm] + [cov(Ym, h) for h in hlf]))

    # offspring (every offspring with its own parents)
    offp = ids.loc[off["Father.ID"].astype(np.int64).to_numpy()]
    offm = ids.loc[off["Mother.ID"].astype(np.int64).to_numpy()]
    Yfo, Ymo = fcol(offp, "Y1"), fcol(offm, "Y1")
    To = [fcol(off, "TPO1") / d, fcol(off, "TMO1") / d]
    NTo = [fcol(off, "NTPO1") / d, fcol(off, "NTMO1") / d]
    Fo = fcol(off, "F1")
    for tag, col in (("", "Y1"), ("_o", "Y2")):          # "" = same trait, "_o" = DiffTrait offspring trait
        Yo = fcol(off, col)
        s[f"VY_off{tag}"] = np.var(Yo, ddof=1)
        s[f"Yo_Yp{tag}"] = np.mean([cov(Yo, Yfo), cov(Yo, Ymo)])
        s[f"thetaT{tag}"] = np.mean([cov(Yo, h) for h in To])
        s[f"thetaNT{tag}"] = np.mean([cov(Yo, h) for h in NTo])
    s.update(Y_Fo=np.mean([cov(Yfo, Fo), cov(Ymo, Fo)]),
             VF_off=np.var(Fo, ddof=1), w_off=2 * np.mean([cov(Fo, h) for h in To]),
             Omega_off=np.mean([cov(fcol(off, "Y1"), h) for h in To]),
             gc_off=cov(*To), hapvar_off=np.mean([np.var(h, ddof=1) for h in To]))
    return {k: float(v) for k, v in s.items()}


def make_cv(prm, rng):
    """CV table: the observed part explains delta^2 and the latent part a^2 of generation-0 variance
    (trait 1); trait 2 has the same CVs and split, rescaled to delta_o / a_o."""
    maf = rng.uniform(MAF_MIN, MAF_MAX, N_CV)
    raw = rng.standard_normal(N_CV)
    mask = np.zeros(N_CV, dtype=int)
    mask[rng.choice(N_CV, N_OBS, replace=False)] = 1
    het = 2 * maf * (1 - maf)
    obs, lat = raw * mask, raw * (1 - mask)
    obs = obs * prm["delta"] / np.sqrt(np.sum(het * obs ** 2))
    lat = lat * prm["a"] / np.sqrt(np.sum(het * lat ** 2))
    return pd.DataFrame({"maf": maf, "alpha1": obs + lat,
                         "alpha2": obs * prm["delta_o"] / prm["delta"] + lat * prm["a_o"] / prm["a"],
                         "mask_obs1": mask, "mask_obs2": mask})


def mating_paths(prm):
    return {"AO": prm["AM_G"], "AL": prm["AM_G"], "F": prm["AM_E"], "E": prm["AM_E"]}


def build_sim(prm, cv, n_gens, am_list, seed, paths):
    kw = dict(cove_mat=np.diag([prm["VE"], prm["VE_o"]]),
              f_mat=np.array([[prm["f"], 0.0], [prm["f"], 0.0]]),   # F = f*(Y1_father + Y1_mother) for both traits
              s_mat=np.zeros((2, 2)), a_mat=np.diag([prm["a"], prm["a_o"]]),
              d_mat=np.diag([prm["delta"], prm["delta_o"]]), covy_mat=np.eye(2), k2_matrix=np.zeros((2, 2)))
    with contextlib.redirect_stdout(io.StringIO()):
        return CascadeValidationSim(cv_info=cv, num_generations=n_gens, pop_size=POP_SIZE, mate_on_trait=1,
                                    mating_paths=paths, am_list=am_list, seed=seed,
                                    save_each_gen=False, save_covs=False, **kw)


def export_trios(par, off, prm, rng):
    """One random offspring per family, centred on the population means of this replicate; plus the binary
    versions, dichotomized at the true thresholds of the VY = 1 liability scale."""
    d = prm["delta"]
    ids = par.assign(ID=par["ID"].astype(np.int64)).set_index("ID")
    fam = off.sample(frac=1.0, random_state=int(rng.integers(2 ** 31))).drop_duplicates(["Father.ID", "Mother.ID"])
    fa = ids.loc[fam["Father.ID"].astype(np.int64).to_numpy()]
    mo = ids.loc[fam["Mother.ID"].astype(np.int64).to_numpy()]
    hap_mean = np.mean([fcol(par, "TPO1").mean(), fcol(par, "TMO1").mean()]) / d
    y1_par, y1_off, y2_off = fcol(par, "Y1").mean(), fcol(off, "Y1").mean(), fcol(off, "Y2").mean()
    base = pd.DataFrame({
        "Yp1": fcol(fa, "Y1") - y1_par, "Ym1": fcol(mo, "Y1") - y1_par,
        "Tp1": fcol(fam, "TPO1") / d - hap_mean, "NTp1": fcol(fam, "NTPO1") / d - hap_mean,
        "Tm1": fcol(fam, "TMO1") / d - hap_mean, "NTm1": fcol(fam, "NTMO1") / d - hap_mean})
    obs_cols = ["Yp1", "Ym1", "Yo1", "Tp1", "NTp1", "Tm1", "NTm1"]
    lat_cols = ["Yo1", "Tp1", "NTp1", "Tm1", "NTm1"]
    same = base.assign(Yo1=fcol(fam, "Y1") - y1_off)[obs_cols]
    diff = base.assign(Yo1=fcol(fam, "Y2") - y2_off)[obs_cols]
    out = {"SameTrait_Observed": same, "SameTrait_Latent": same[lat_cols],
           "DiffTrait_Observed": diff, "DiffTrait_Latent": diff[lat_cols]}
    for name in list(out):
        b = out[name].copy()
        for c in ("Yp1", "Ym1"):
            if c in b:
                b[c] = (b[c] > THRESH["parent"]).astype(int)
        b["Yo1"] = (b["Yo1"] > THRESH["offspring_diff" if name.startswith("Diff") else "parent"]).astype(int)
        out[name + "_binary"] = b
    return out


def save_part(regime, mech, rep, seed, sim_traj, final, prm, rng, t0):
    tag = f"{regime}_{mech}_rep{rep:03d}"
    par, off = final
    trios = export_trios(par, off, prm, rng)
    for name, df in trios.items():
        df.to_csv(os.path.join(PARTS, f"{tag}_{name}.tsv"), sep="\t", index=False)
    traj = pd.DataFrame(sim_traj).assign(regime=regime, mechanism=mech, rep=rep)
    traj.to_csv(os.path.join(PARTS, f"{tag}_trajectory.csv"), index=False)
    pd.DataFrame([{"regime": regime, "mechanism": mech, "rep": rep, "seed": seed,
                   "n_trios": len(trios["SameTrait_Observed"]), "time_sec": round(time.time() - t0, 1)}]
                 ).to_csv(os.path.join(PARTS, f"{tag}_info.csv"), index=False)
    return tag


def run_eq(args):
    mech, rep, prm = args
    seed = SEED_BASE + 1000 * MECHS.index(mech) + rep
    if os.path.exists(os.path.join(PARTS, f"Eq_{mech}_rep{rep:03d}_info.csv")):
        return f"Eq_{mech}_rep{rep:03d} (exists)"
    t0 = time.time()
    rng = np.random.default_rng(seed)
    sim = build_sim(prm, make_cv(prm, rng), EQ_GENS, [prm["am"]] * EQ_GENS, seed, mating_paths(prm))
    sim.setup_recording(prm, last_gen=EQ_GENS - 1)
    with contextlib.redirect_stdout(io.StringIO()):
        sim.run_simulation()
    return save_part("Eq", mech, rep, seed, sim.traj, sim.final, prm, rng, t0)


def run_diseq(args):
    rep, prms = args
    seed = SEED_BASE + 5000 + rep
    if all(os.path.exists(os.path.join(PARTS, f"DisEq_{m}_rep{rep:03d}_info.csv")) for m in MECHS):
        return f"DisEq rep{rep:03d} (exists)"
    t0 = time.time()
    rng = np.random.default_rng(seed)
    p0 = prms["primary"]    # the random-mating history does not depend on the mechanism
    for m in MECHS:
        assert all(abs(prms[m][k] - p0[k]) < 1e-12 for k in ("delta", "a", "VE", "f"))
    # (VE_o, the DiffTrait offspring trait's residual, is mechanism-specific; it only matters in the offspring
    # generation, so each branch sets it just before its AM generation)
    sim = build_sim(p0, make_cv(p0, rng), DISEQ_BURN + 1, [0.0] * DISEQ_BURN + [p0["am"]], seed, mating_paths(p0))
    sim.setup_recording(p0, last_gen=DISEQ_BURN)
    with contextlib.redirect_stdout(io.StringIO()):
        for g in range(DISEQ_BURN):                     # random mating (am = 0) with vertical transmission
            sim._run_generation(gen_abs_idx=g, am_idx=g, pop_size_target=POP_SIZE, save=False)
    tags = []
    for i, m in enumerate(MECHS):                       # branch: one AM generation on gamma~ per mechanism
        br = copy.deepcopy(sim)
        br.prm = prms[m]
        br.mating_paths = mating_paths(prms[m])
        br._mating_weights = br._resolve_mating_paths(br.mating_paths)
        br.cove_mat = np.diag([prms[m]["VE"], prms[m]["VE_o"]])
        np.random.seed(seed * 10 + i)
        with contextlib.redirect_stdout(io.StringIO()):
            br._run_generation(gen_abs_idx=DISEQ_BURN, am_idx=DISEQ_BURN, pop_size_target=POP_SIZE, save=False)
        tags.append(save_part("DisEq", m, rep, seed, br.traj, br.final, prms[m], rng, t0))
    return ", ".join(tags)


def merge():
    parts = sorted(os.listdir(PARTS))
    for regime in ("Eq", "DisEq"):
        for m in MECHS:
            for name in ("SameTrait_Observed", "SameTrait_Latent", "DiffTrait_Observed", "DiffTrait_Latent"):
                for suffix in ("", "_binary"):
                    files = [f for f in parts if f.startswith(f"{regime}_{m}_rep") and f.endswith(f"_{name}{suffix}.tsv")]
                    if files:
                        pd.concat([pd.read_csv(os.path.join(PARTS, f), sep="\t") for f in files], ignore_index=True
                                  ).to_csv(os.path.join(OUT_DATA, f"{regime}_{m}_{name}{suffix}.tsv"), sep="\t", index=False)
    traj = pd.concat([pd.read_csv(os.path.join(PARTS, f)) for f in parts if f.endswith("_trajectory.csv")], ignore_index=True)
    traj.to_csv(os.path.join(OUT_SIM, "trajectory.csv"), index=False)
    last = traj.loc[traj.groupby(["regime", "mechanism", "rep"])["gen"].idxmax()]
    last.melt(id_vars=["regime", "mechanism", "rep", "gen"], var_name="quantity", value_name="value").dropna(
        subset=["value"]).to_csv(os.path.join(OUT_SIM, "moments_final.csv"), index=False)
    pd.concat([pd.read_csv(os.path.join(PARTS, f)) for f in parts if f.endswith("_info.csv")], ignore_index=True
              ).to_csv(os.path.join(OUT_SIM, "rep_info.csv"), index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=REPS)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--only", choices=["Eq", "DisEq"], default=None, help="simulate one regime only (existing parts are kept)")
    args = ap.parse_args()
    for dname in (OUT_DATA, OUT_SIM, PARTS):
        os.makedirs(dname, exist_ok=True)
    gp = pd.read_csv(os.path.join(HERE, "output", "truth", "generating_params.csv"))
    prm = {(r.regime, r.mechanism): r._asdict() for r in gp.itertuples(index=False)}
    tasks = [(m, rep, prm[("Eq", m)]) for rep in range(1, args.reps + 1) for m in MECHS]
    dtasks = [(rep, {m: prm[("DisEq", m)] for m in MECHS}) for rep in range(1, args.reps + 1)]
    print(f"POP_SIZE={POP_SIZE} REPS={args.reps} EQ_GENS={EQ_GENS} DISEQ_BURN={DISEQ_BURN} N_CV={N_CV} "
          f"workers={args.workers}", flush=True)
    t0 = time.time()
    with Pool(args.workers) as pool:
        if args.only != "Eq":
            for msg in pool.imap_unordered(run_diseq, dtasks):
                print(f"[{time.time() - t0:6.0f}s] {msg}", flush=True)
        if args.only != "DisEq":
            for msg in pool.imap_unordered(run_eq, tasks):
                print(f"[{time.time() - t0:6.0f}s] {msg}", flush=True)
    merge()
    print(f"merged -> {OUT_DATA}, {OUT_SIM} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
