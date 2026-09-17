#!/usr/bin/env python3
"""Calibrate each trait's assortment strength so that both traits reach the target phenotypic
spousal correlation R_SPOUSE_Y in the parental generation of an N_GENERATIONS run (not just at
generation 0, where the closed-form calibration ignores AM-induced variance inflation).

Bisection on each trait's mating-variable correlation; one simulation per iteration evaluates
both traits at once. Same founders (seed, CV effects) every iteration, so only AM changes.

Run from anywhere:  python 00_calibrate_am.py [--iters 8] [--pop-size 8000] [--seed 1]
Then paste the printed AM values into sim_conditions.py.
"""
import argparse
import contextlib
import importlib
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
S = importlib.import_module("01_simulate")
C = S.C
from core_simulation import AssortativeMatingSimulation  # noqa: E402


def parental_stats(res, par, i):
    off = res["PHEN"].drop_duplicates(["Father.ID", "Mother.ID"])
    par = par.copy()
    par["ID"] = par["ID"].astype(np.int64)
    par = par.set_index("ID")
    fa = par.loc[off["Father.ID"].astype(np.int64).values]
    mo = par.loc[off["Mother.ID"].astype(np.int64).values]
    Y = S.fcol(par, f"Y{i}")
    A = S.fcol(par, f"AO{i}") + S.fcol(par, f"AL{i}")
    return dict(r_Y=float(np.corrcoef(S.fcol(fa, f"Y{i}"), S.fcol(mo, f"Y{i}"))[0, 1]),
                VY=float(np.var(Y, ddof=1)), h2=float(np.var(A, ddof=1) / np.var(Y, ddof=1)))


def run(am, pop_size, seed):
    np.random.seed(seed)
    cv = AssortativeMatingSimulation.prepare_CV_random_selection(
        C.N_CV, 0.0, C.MAF_MIN, C.MAF_MAX, C.PROP_OBS, C.PROP_OBS)
    cv[["alpha1", "alpha2"]] *= np.sqrt(C.VG)
    with contextlib.redirect_stdout(io.StringIO()):
        sim = S.TrioExportSimulation(
            cv_info=cv, num_generations=C.N_GENERATIONS, pop_size=pop_size,
            mating_type=(C.TRAITS[1], C.TRAITS[2]), mate_on_trait=None,
            am_list=[np.diag(am)] * C.N_GENERATIONS,
            seed=seed, save_each_gen=False, save_covs=False, **S.model_matrices())
        res = sim.run_simulation()
    return [parental_stats(res, sim.parent_phen_df, i) for i in (1, 2)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iters", type=int, default=8)
    ap.add_argument("--pop-size", type=int, default=8000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    target = C.R_SPOUSE_Y
    names = [C.TRAITS[1], C.TRAITS[2]]
    lo, hi = np.array([0.05, 0.05]), np.array([0.95, 0.95])
    print(f"Target r(Y_p, Y_m) = {target:.3f} at generation {C.N_GENERATIONS}; "
          f"pop {args.pop_size}, seed {args.seed}, {args.iters} bisection steps\n")
    print(f"{'iter':>4}  " + "   ".join(f"{n:>10}: AM    r_Y    VY    h2 " for n in names))
    mid = None
    for it in range(1, args.iters + 1):
        mid = (lo + hi) / 2
        out = run(mid, args.pop_size, args.seed)
        print(f"{it:>4}  " + "   ".join(
            f"{names[j]:>10}: {mid[j]:.3f} {o['r_Y']:.3f} {o['VY']:.3f} {o['h2']:.3f}" for j, o in enumerate(out)))
        for j, o in enumerate(out):
            if o["r_Y"] < target:
                lo[j] = mid[j]
            else:
                hi[j] = mid[j]
    final = (lo + hi) / 2
    print("\nPaste into sim_conditions.py:")
    print("AM = {" + ", ".join(f'"{n}": {v:.3f}' for n, v in zip(names, final)) + "}")


if __name__ == "__main__":
    main()
