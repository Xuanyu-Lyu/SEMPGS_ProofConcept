#!/usr/bin/env python3
"""Step 6: a replication set for the one shift that survived Step 4 (Eq, social homogamy, SameTrait_Observed).

Simulates populations 13..36 of the Eq social-homogamy condition exactly as 01_simulate.py does (same code, same
seed scheme, so rep r here is the population 01_simulate.py --reps 36 would have made), but writes them to
output/replication/parts, so the pooled data and fits of the main validation stay as they are.
06_replicate_eq_social.R then fits them.

Run with: python 06_replicate_eq_social.py [--workers 8]
"""
import argparse
import importlib.util
import os
import time
from multiprocessing import Pool

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REP_PARTS = os.path.join(HERE, "output", "replication", "parts")
FIRST, LAST = 13, 36


def load_sim_module():
    spec = importlib.util.spec_from_file_location("simulate", os.path.join(HERE, "01_simulate.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.PARTS = REP_PARTS          # save_part() and run_eq() write here instead of output/sim/parts
    return mod


def init_worker():
    global SIM
    SIM = load_sim_module()


def run(args):
    return SIM.run_eq(args)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(REP_PARTS, exist_ok=True)
    gp = pd.read_csv(os.path.join(HERE, "output", "truth", "generating_params.csv"))
    prm = [r._asdict() for r in gp.itertuples(index=False) if r.regime == "Eq" and r.mechanism == "social"][0]
    tasks = [("social", rep, prm) for rep in range(FIRST, LAST + 1)]
    t0 = time.time()
    with Pool(args.workers, initializer=init_worker) as pool:
        for msg in pool.imap_unordered(run, tasks):
            print(f"[{time.time() - t0:6.0f}s] {msg}", flush=True)
    print(f"done ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
