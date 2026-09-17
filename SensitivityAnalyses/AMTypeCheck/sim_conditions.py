"""Simulation conditions for the AM-type sensitivity analysis. Edit values here only.

Two genetically and environmentally independent traits are simulated in one population:
trait 1 mates on the total genetic value (genetic AM), trait 2 on the familial + unique
environment (social AM). The traits share every generating value except the strength of
assortment on the mating variable, which is calibrated per trait so that both traits have
(approximately) the same phenotypic spousal correlation R_SPOUSE_Y.
"""

N_REPS = 100
SEED_BASE = 20260917          # replicate r uses seed SEED_BASE + r

POP_SIZE = 4000               # individuals per generation; ~0.43 * POP_SIZE independent trios are exported
N_GENERATIONS = 15            # generations of AM + VT before the trio sample is taken (no burn-in)
N_CV = 1000                   # causal variants; int(N_CV * PROP_OBS) of them form the observed PGS
MAF_MIN, MAF_MAX = 0.01, 0.50

VG = 0.5                      # generation-0 additive genetic variance per trait (CV effects are rescaled to this)
PROP_OBS = 0.4                # observed-PGS share of VG  ->  delta = sqrt(VG*PROP_OBS), a = sqrt(VG*(1-PROP_OBS))
VE = 0.5                      # residual variance -> generation-0 VY = VG + VE = 1, h2 = .5
F = 0.15                      # vertical-transmission path from each parental phenotype

R_SPOUSE_Y = 0.40             # target phenotypic spousal correlation, same for both traits
TRAITS = {1: "genotypic", 2: "social"}   # trait index -> what spouses assort on

# Correlation of the mating variable per AM type, calibrated with 00_calibrate_am.py so that the
# phenotypic spousal correlation in the parental generation of an N_GENERATIONS run is ~R_SPOUSE_Y
# for both traits (a generation-0 formula fails because genetic AM inflates VA over generations).
# Re-run the calibration if VG, PROP_OBS, VE, F, N_GENERATIONS or R_SPOUSE_Y change.
AM = {"genotypic": 0.547, "social": 0.649, "phenotypic": R_SPOUSE_Y}
