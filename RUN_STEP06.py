"""Step 06: nested persistence, baseline levels and ocular sensitivity models.

Logistic GLMs use uncorrected sandwich covariance (HC0). Ocular linear models
use leverage-adjusted HC3 covariance and t inference. Holm correction is applied
to the 12 ocular contrasts within each threshold and adjustment specification.
These models describe exploratory associations, not validated prediction.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from dream_common import ASSAYS, NOVEL, ANTIGENS, OUTCOMES, fraction, require, run_cli

def compute(data):
    m0, m12 = data.pairs()
    baseline = data.baseline()
    require(len(baseline) == 492 and baseline.loc[m0.index, ASSAYS].eq(m0).all().all(),
            "Baseline source or pairing does not agree with the frozen cohort")
    clinical = data.read("gi.xlsx", ["age", "sex", "group", "treatment"]).join(
        data.read("mh.xlsx", ["mhsjog", "mhrheu", "mhthyr", "mhdiab"])).join(data.eyes("00")).join(data.eyes("12"))
    require(set(clinical.sex) == {"F", "M"} and clinical.treatment.eq(
        clinical.group.map({0: "PLACEBO", 1: "ACTIVE"})).all(), "Unexpected clinical coding")
    clinical["female"] = clinical.sex.eq("F").astype(int)
    clinical["age10"] = clinical.age / 10
    for field in ("mhsjog", "mhrheu", "mhthyr", "mhdiab"):
        require(clinical[field].isin([1, 2, 3]).all(), f"Unexpected disease history in {field}")
        clinical[field] = clinical[field].eq(3).astype(int)
    clinical["log2_max"] = np.log2(baseline[NOVEL].max(axis=1))
    clinical["ssa"] = baseline.ss_a_ro.gt(25).where(baseline.ss_a_ro.notna()).astype(float)
    clinical["panel_strict"] = baseline[NOVEL].gt(25).any(axis=1).astype(int)
    retained = clinical.loc[m0.index]
    missing = clinical.loc[baseline.index.difference(m0.index)]
    balance = []
    for field in ("age", "female", "group", "mhsjog", "mhrheu", "mhthyr", "mhdiab", "ssa", "panel_strict",
                  "osdi_00", "schirmer_00", "conjunctival_staining_00", "corneal_staining_00", "log2_max"):
        a, b = retained[field].dropna(), missing[field].dropna()
        sd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        balance.append({"feature": field, "paired_available": len(a), "no_complete_followup_available": len(b),
                        "paired_mean": a.mean(), "no_complete_followup_mean": b.mean(),
                            "paired_sd": a.std(ddof=1), "no_complete_followup_sd": b.std(ddof=1),
                        "smd_paired_minus_missing": (a.mean() - b.mean()) / sd if sd else np.nan})
    tiers, ssb, arm_rows, titer_rows, ocular, changes = [], [], [], [], [], []
    for rule, a, b in (("GT25", m0[NOVEL].gt(25), m12[NOVEL].gt(25)),
                        ("GE20", m0[NOVEL].ge(20), m12[NOVEL].ge(20))):
        a, b = a.astype(bool), b.astype(bool)
        pos0, pos12 = a.any(axis=1), b.any(axis=1)
        same_assay = (a & b).any(axis=1)
        same_antigen = pd.DataFrame({antigen: a[fields].any(axis=1) & b[fields].any(axis=1)
                                     for antigen, fields in ANTIGENS.items()}).any(axis=1)
        require((~same_assay | same_antigen).all() and (~same_antigen | (pos0 & pos12)).all(),
                "Same-assay, antigen and panel categories are not nested")
        base_n = int(pos0.sum())
        for label, event in (("same_antigen_and_isotype", same_assay), ("same_antigen_any_isotype", same_antigen),
                             ("any_candidate_panel", pos0 & pos12)):
            tiers.append({"rule": rule, "definition": label, **fraction(int(event.sum()), base_n)})
        for ssb_rule, ssb_positive in (("DREAM_GT25", m0.ss_b_la.gt(25)), ("IMMCO_046_GE50", m0.ss_b_la.ge(50))):
            selected = ~(m0.ss_a_ro.gt(25) | ssb_positive)
            ssb.append({"rule": rule, "ssb_rule": ssb_rule, "ssb_positive_all_pairs": int(ssb_positive.sum()),
                        "subgroup_n": int(selected.sum()), "baseline_panel_positive": int((selected & pos0).sum()),
                        "same_assay_retained": int((selected & same_assay).sum()),
                        "same_antigen_retained": int((selected & same_antigen).sum()),
                        "panel_retained": int((selected & pos0 & pos12).sum())})
        for arm, arm_data in retained.groupby("treatment"):
            select = pos0.index.isin(arm_data.index)
            arm_rows.append({"rule": rule, "arm": arm, "pairs": len(arm_data), "base_positive": int((pos0 & select).sum()),
                             "same_assay_retained": int((same_assay & select).sum()),
                             "same_antigen_retained": int((same_antigen & select).sum()),
                             "panel_retained": int((pos0 & pos12 & select).sum()),
                             "new_positive": int((~pos0 & pos12 & select).sum())})

        strata = pd.Series("both_nonpositive", index=m0.index)
        strata.loc[~pos0 & pos12] = "new_positive"
        strata.loc[pos0 & ~pos12] = "lost"
        strata.loc[pos0 & pos12 & ~same_assay] = "changed_assay"
        strata.loc[same_assay] = "retained_same"
        fit_data = retained.copy()
        fit_data["multiple"] = a.sum(axis=1).ge(2).astype(int)
        fit_data["trajectory"] = strata
        fit_data["same_assay"] = same_assay.astype(int)
        for label, selected in (("all_baseline_positive", pos0),
                                ("exclude_right_censored_max", pos0 & m0[NOVEL].max(axis=1).lt(160))):
            x = fit_data.loc[selected]
            terms = ["log2_max", "multiple", "group"]
            fit = sm.GLM(x.same_assay, sm.add_constant(x[terms]), family=sm.families.Binomial()).fit(cov_type="HC0")
            require(fit.converged, "Titer association model did not converge")
            for term in ("log2_max", "multiple"):
                lo, hi = fit.conf_int().loc[term]
                titer_rows.append({"rule": rule, "analysis_set": label, "term": term, "n": len(x),
                                   "events": int(x.same_assay.sum()), "odds_ratio": np.exp(fit.params[term]),
                                   "ci_low": np.exp(lo), "ci_high": np.exp(hi)})
        baseline_max = m0[NOVEL].max(axis=1)
        for band, selected in (("one_positive_item", pos0 & ~fit_data.multiple.astype(bool)),
                               ("two_or_more_positive_items", pos0 & fit_data.multiple.astype(bool)),
                               ("max_25_to_35", pos0 & baseline_max.gt(25) & baseline_max.le(35)),
                               ("max_35_to_50", pos0 & baseline_max.gt(35) & baseline_max.le(50)),
                               ("max_above_50", pos0 & baseline_max.gt(50))):
            changes.append({"rule": rule, "band": band,
                            **fraction(int((selected & same_assay).sum()), int(selected.sum()))})
        for outcome in OUTCOMES:
            base, final = f"{outcome}_00", f"{outcome}_12"
            for origin, control, target in (("baseline_positive", "lost", ["retained_same", "changed_assay"]),
                                            ("baseline_nonpositive", "both_nonpositive", ["new_positive"])):
                x = fit_data.loc[fit_data.trajectory.isin([control, *target])]
                for model in ("baseline_and_arm", "plus_demographics_and_antibody"):
                    terms = [base, "group"]
                    if model == "plus_demographics_and_antibody":
                        terms.extend(["age10", "female", "log2_max", "ssa"])
                        if origin == "baseline_positive":
                            terms.append("multiple")
                    observed = x[[*terms, final, "trajectory"]].dropna()
                    design = observed[terms].astype(float).copy()
                    for label in target:
                        design[label] = observed.trajectory.eq(label).astype(int)
                    fit = sm.OLS(observed[final].astype(float), sm.add_constant(design)).fit(cov_type="HC3", use_t=True)
                    for label in target:
                        lo, hi = fit.conf_int().loc[label]
                        ocular.append({"rule": rule, "outcome": outcome, "model": model, "origin": origin,
                                       "contrast": label, "reference": control,
                                       "n_total": len(observed), "n_contrast": int(observed.trajectory.eq(label).sum()),
                                       "n_reference": int(observed.trajectory.eq(control).sum()),
                                       "adjusted_difference": fit.params[label], "ci_low": lo, "ci_high": hi,
                                       "p_nominal": fit.pvalues[label]})
    ocular = pd.DataFrame(ocular)
    for _, group in ocular.groupby(["rule", "model"]):
        ocular.loc[group.index, "p_holm_across_12"] = multipletests(group.p_nominal, method="holm")[1]
    return {"persistence_tiers": pd.DataFrame(tiers), "ssb_cutoff_sensitivity": pd.DataFrame(ssb),
            "by_randomized_arm": pd.DataFrame(arm_rows), "titer_and_count_models": pd.DataFrame(titer_rows),
            "titer_and_count_bands": pd.DataFrame(changes), "eye_outcome_contrasts": ocular,
            "baseline_missingness_profile": pd.DataFrame(balance)}

if __name__ == "__main__":
    run_cli(compute, "STEP06_models_and_sensitivity.xlsx", __doc__)
