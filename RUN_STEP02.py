"""Step 02: reconstruct the complete paired antibody cohort from raw data."""
import pandas as pd
from dream_common import ASSAYS, KEY, require, run_cli

def compute(data):
    gi = data.read_sheet("gi.xlsx", [KEY])
    raw = data.numeric_lab()
    complete = raw.dropna(subset=ASSAYS)
    records = {v: set(raw.loc[raw.month.eq(v), KEY]) for v in ("00", "12", "24")}
    usable = {v: set(complete.loc[complete.month.eq(v), KEY]) for v in ("00", "12", "24")}
    paired = usable["00"] & usable["12"]
    base, follow = data.pairs()
    require(set(base.index) == paired, "Paired membership differs between reconstructions")
    n0, n12, both = len(usable["00"]), len(usable["12"]), len(paired)
    flow = {
        "trial_participants_in_GI": len(gi), "SJO_records_all_visits": len(raw),
        "people_with_any_SJO_record": raw[KEY].nunique(),
        "baseline_records": len(records["00"]), "baseline_records_with_all_11_assays": n0,
        "baseline_records_without_usable_results": len(records["00"])-n0,
        "month12_records": len(records["12"]), "month12_records_with_all_11_assays": n12,
        "month12_records_without_usable_results": len(records["12"])-n12,
        "people_with_records_at_both_visits": len(records["00"] & records["12"]),
        "complete_baseline_and_month12_pairs": both,
        "baseline_complete_without_complete_month12": n0-both,
        "month12_complete_without_complete_baseline": n12-both,
        "neither_visit_complete": len(gi)-len(usable["00"] | usable["12"]),
        "month24_records_outside_primary_pair": len(records["24"]),
    }
    days = data.collection_intervals()
    intervals = dict(n=len(days), median_days=days.median(), q1_days=days.quantile(.25),
                     q3_days=days.quantile(.75), min_days=days.min(), max_days=days.max())
    return {"cohort_flow": pd.DataFrame(flow.items(), columns=["item", "n"]),
            "collection_interval": pd.DataFrame([intervals])}

if __name__ == "__main__":
    run_cli(compute, "STEP02_cohort.xlsx", __doc__)
