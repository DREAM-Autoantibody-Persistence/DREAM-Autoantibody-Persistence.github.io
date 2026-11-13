"""Step 05: assemble the main manuscript tables as readable Excel worksheets."""
import pandas as pd
from dream_common import run_cli
from RUN_STEP02 import compute as cohort_tables
from RUN_STEP03 import compute as main_tables

def compute(data):
    main = main_tables(data)
    rates = main["main_proportions_wilson"]
    changes = main["paired_panel_change"]
    display, persistence = [], []
    for rule, label in (("strict", ">25 EU/mL"), ("incl_borderline", ">=20 EU/mL")):
        row = rates.loc[rates.rule.eq(rule) & rates.estimate.eq("panel_new_positivity")].iloc[0]
        change = changes.loc[changes.rule.eq(rule)].iloc[0]
        display.extend([
            dict(threshold=label, measure="Newly met panel threshold at month 12",
                 result=f"{int(row.numerator)}/{int(row.denominator)} ({row.percent:.1f}%)",
                 ci95=f"{row.ci_low_percent:.1f}% to {row.ci_high_percent:.1f}%"),
            dict(threshold=label, measure="Paired net change, percentage points",
                 result=f"{change.difference_pp:.1f}",
                 ci95=f"{change.ci_low_pp:.1f} to {change.ci_high_pp:.1f}"),
        ])
        for _, item in rates.loc[rates.rule.eq(rule)].iterrows():
            persistence.append(dict(threshold=label, measure=item.estimate,
                n=int(item.numerator), denominator=int(item.denominator), percent=item.percent,
                ci_low_percent=item.ci_low_percent, ci_high_percent=item.ci_high_percent))
    return {"Table1_panel_change": pd.DataFrame(display),
            "main_persistence": pd.DataFrame(persistence),
            "cohort_flow": cohort_tables(data)["cohort_flow"]}

if __name__ == "__main__":
    run_cli(compute, "STEP05_manuscript_tables.xlsx", __doc__)
