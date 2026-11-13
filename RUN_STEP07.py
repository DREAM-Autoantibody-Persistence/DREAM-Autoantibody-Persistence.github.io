"""Step 07: exact two-sided McNemar tests of panel-state transitions."""
import pandas as pd
from scipy.stats import binomtest
from dream_common import NOVEL, transition, run_cli

def compute(data):
    base, follow = data.pairs()
    rows = []
    for rule, a, b in (("GT25", base[NOVEL].gt(25), follow[NOVEL].gt(25)),
                       ("GE20", base[NOVEL].ge(20), follow[NOVEL].ge(20))):
        row = transition(rule, "candidate_panel", a.any(axis=1), b.any(axis=1))
        gain, loss = row["nonpos_to_pos"], row["pos_to_nonpos"]
        row["discordant_pairs"] = gain+loss
        row["net_change_pp"] = 100*(gain-loss)/len(base)
        row["p_exact_two_sided"] = binomtest(gain, gain+loss, p=.5).pvalue if gain+loss else 1.
        rows.append(row)
    return {"exact_mcnemar": pd.DataFrame(rows)}

if __name__ == "__main__":
    run_cli(compute, "STEP07_mcnemar.xlsx", __doc__)
