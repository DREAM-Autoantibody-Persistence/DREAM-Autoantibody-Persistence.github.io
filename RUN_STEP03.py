"""Step 03: transitions, persistence and paired prevalence changes.

Net-change intervals use 100,000 participant-pair bootstrap draws. A multinomial
draw of {-1,0,+1} scores is equivalent to resampling the observed whole pairs.
"""
import numpy as np
import pandas as pd
from dream_common import ASSAYS, fraction, transition, run_cli

def compute(data):
    m0, m12 = data.pairs()
    novel = ASSAYS[2:]
    ssa0, ssa12 = m0.ss_a_ro.gt(25), m12.ss_a_ro.gt(25)
    ssa_ssb_nonpositive = ~(ssa0 | m0.ss_b_la.gt(25))
    transitions, rates, patterns, differences = [], [], [], []
    rng = np.random.default_rng(20261108)
    for rule, first, second in (("strict", m0[novel].gt(25), m12[novel].gt(25)),
                                ("incl_borderline", m0[novel].ge(20), m12[novel].ge(20))):
        p0, p12 = first.any(axis=1), second.any(axis=1)
        same = (first & second).any(axis=1)
        ssa_ref = transition(rule, "anti_SSA", ssa0, ssa12)
        panel = transition(rule, "novel_any", p0, p12)
        subgroup = transition(rule, "novel_any_in_baseline_SSA_SSB_nonpositive",
                              p0[ssa_ssb_nonpositive], p12[ssa_ssb_nonpositive])
        same_ref = transition(rule, "novel_same_isotype_in_baseline_panel_positive", p0[p0], same[p0])
        transitions.extend([ssa_ref, panel, subgroup, same_ref])
        base_n = panel["repeat_pos"] + panel["pos_to_nonpos"]
        pattern = {
            "rule": rule, "baseline_panel_positive": base_n,
            "retain_at_least_one_same_assay": int(same.sum()),
            "panel_positive_with_only_different_assays": int((p0 & p12 & ~same).sum()),
            "no_longer_panel_positive": panel["pos_to_nonpos"],
            "new_assay_among_baseline_panel_positive": int((p0 & (~first & second).any(axis=1)).sum()),
            "all_nine_classifications_unchanged_among_baseline_panel_positive":
                int((p0 & first.eq(second).all(axis=1)).sum()),
            "same_assay_in_baseline_ssa_ssb_nonpositive": int((same & ssa_ssb_nonpositive).sum()),
        }
        patterns.append(pattern)
        for label, num, denom in (
            ("panel_persistence", panel["repeat_pos"], base_n),
            ("same_assay_persistence", pattern["retain_at_least_one_same_assay"], base_n),
            ("panel_new_positivity", panel["nonpos_to_pos"], len(m0) - base_n),
            ("anti_ssa_persistence", ssa_ref["repeat_pos"], ssa_ref["repeat_pos"] + ssa_ref["pos_to_nonpos"]),
            ("same_assay_persistence_in_baseline_ssa_ssb_nonpositive",
             pattern["same_assay_in_baseline_ssa_ssb_nonpositive"], int((p0 & ssa_ssb_nonpositive).sum())),
        ):
            rates.append({"rule": rule, "estimate": label, **fraction(num, denom)})
        n = panel["n"]
        gain, loss = panel["nonpos_to_pos"], panel["pos_to_nonpos"]
        paired = rng.multinomial(n, np.array([loss, gain, n - loss - gain]) / n, size=100_000)
        lo, hi = np.quantile((paired[:, 1] - paired[:, 0]) / n, [0.025, 0.975])
        differences.append({"rule": rule, "n": n, "baseline_positive": int(p0.sum()),
                            "followup_positive": int(p12.sum()), "gain": gain, "loss": loss,
                            "difference_pp": 100 * (gain - loss) / n, "ci_low_pp": float(100 * lo),
                            "ci_high_pp": float(100 * hi)})

    transition_table = pd.DataFrame(transitions)
    difference_table = pd.DataFrame(differences)
    return {"main_transitions": transition_table, "main_proportions_wilson": pd.DataFrame(rates),
            "panel_patterns": pd.DataFrame(patterns), "paired_panel_change": difference_table}

if __name__ == "__main__":
    run_cli(compute, "STEP03_main_results.xlsx", __doc__)
