"""Step 08: exploratory comparisons by original randomized assignment.

Net-change differences resample participant pairs independently within each arm.
Same-assay persistence differences use Newcombe confidence intervals.
"""
import numpy as np
import pandas as pd
from statsmodels.stats.proportion import confint_proportions_2indep
from dream_common import NOVEL, fraction, require, run_cli

ARMS = ("ACTIVE", "PLACEBO")
BOOTSTRAP_SEED = 20261110
BOOTSTRAP_DRAWS = 100_000

def estimate_arm(first, second, rng):
    """Resample whole participant pairs via their three possible change scores."""
    baseline, followup = first.any(axis=1), second.any(axis=1)
    gain = int((~baseline & followup).sum())
    loss = int((baseline & ~followup).sum())
    retained = int((first & second).any(axis=1).sum())
    n = len(baseline)
    base_n = int(baseline.sum())
    boot = rng.multinomial(n, np.array([loss, gain, n - gain - loss]) / n, size=BOOTSTRAP_DRAWS)
    draws = 100 * (boot[:, 1] - boot[:, 0]) / n
    low, high = np.quantile(draws, [0.025, 0.975])
    return {"paired_n": n, "baseline_meets": base_n, "month12_meets": int(followup.sum()),
            "gain": gain, "loss": loss, "net_change_pp": 100 * (gain - loss) / n,
            "net_ci_low_pp": low, "net_ci_high_pp": high}, fraction(retained, base_n), draws

def compute(data):
    baseline, followup = data.pairs()
    gi = data.read("gi.xlsx", ["group", "treatment"])
    require(gi.treatment.eq(gi.group.map({0: "PLACEBO", 1: "ACTIVE"})).all(),
            "Randomized arm labels and codes do not agree")
    require(gi.treatment.value_counts().to_dict() == {"ACTIVE": 349, "PLACEBO": 186},
            "Randomization counts differ from the original DREAM trial report")
    complete_at_baseline = data.baseline().index
    require(len(complete_at_baseline) == len(data.baseline())
            and baseline.index.isin(complete_at_baseline).all(), "Baseline assay coverage differs from Step 02")
    coverage = []
    paired_arm = gi.loc[baseline.index, "treatment"]
    for arm in ARMS:
        randomized_n = int(gi.treatment.eq(arm).sum())
        baseline_n = int(gi.loc[complete_at_baseline, "treatment"].eq(arm).sum())
        paired_n = int(paired_arm.eq(arm).sum())
        coverage.append({"arm": arm, "randomized_n": randomized_n, "baseline_complete_n": baseline_n,
                         "paired_n": paired_n, "paired_pct_of_randomized": 100 * paired_n / randomized_n,
                         "paired_pct_of_baseline_complete": 100 * paired_n / baseline_n})

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    arm_rows, contrasts = [], []
    for rule, first, second in (("GT25", baseline[NOVEL].gt(25), followup[NOVEL].gt(25)),
                                ("GE20", baseline[NOVEL].ge(20), followup[NOVEL].ge(20))):
        results = [estimate_arm(first.loc[paired_arm.eq(arm)], second.loc[paired_arm.eq(arm)], rng)
                   for arm in ARMS]
        for arm, (net, same, _) in zip(ARMS, results):
            arm_rows.append({"rule": rule, "arm": arm, **net,
                             **{f"same_assay_{key}": value for key, value in same.items()}})
        (active, active_same, active_draws), (placebo, placebo_same, placebo_draws) = results
        low, high = np.quantile(active_draws - placebo_draws, [0.025, 0.975])
        contrasts.append({"rule": rule, "estimate": "panel_net_change", "ci_method": "arm_stratified_pair_bootstrap",
                          "active_minus_placebo_pp": active["net_change_pp"] - placebo["net_change_pp"],
                          "ci_low_pp": low, "ci_high_pp": high})
        low, high = confint_proportions_2indep(active_same["numerator"], active_same["denominator"],
                                              placebo_same["numerator"], placebo_same["denominator"],
                                              compare="diff", method="newcomb")
        contrasts.append({"rule": rule, "estimate": "same_assay_persistence", "ci_method": "Newcombe",
                          "active_minus_placebo_pp": active_same["percent"] - placebo_same["percent"],
                          "ci_low_pp": 100 * low, "ci_high_pp": 100 * high})
    return {"arm_coverage": pd.DataFrame(coverage), "arm_estimates": pd.DataFrame(arm_rows),
            "arm_contrasts": pd.DataFrame(contrasts)}

if __name__ == "__main__":
    run_cli(compute, "STEP08_randomized_arms.xlsx", __doc__)
