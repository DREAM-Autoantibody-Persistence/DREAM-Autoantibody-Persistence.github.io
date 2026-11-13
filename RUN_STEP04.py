"""Step 04: assay-specific, missing-data and clinical-context summaries.

Historical counts are a published aggregate comparison, not data used to infer
the longitudinal results. Concentration summaries preserve censoring intervals.
Context and directional clinical checks are exploratory supporting outputs.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.proportion import proportion_confint, confint_proportions_2indep
from dream_common import ASSAYS, NOVEL, ANTIGENS, run_cli

def flatten_results(mapping):
    rows = []
    def visit(value, path):
        if isinstance(value, dict):
            for key, item in value.items(): visit(item, path+[str(key)])
        elif isinstance(value, list):
            for key, item in enumerate(value): visit(item, path+[str(key)])
        else: rows.append({"measure": ".".join(path), "value": value})
    visit(mapping, [])
    return pd.DataFrame(rows)

def compute(data):
    m0, m12 = data.pairs()
    novel, markers, groups = NOVEL, ASSAYS, ANTIGENS
    raw = data.numeric_lab()
    read_sheet = data.read_sheet
    tables = {}
    base = raw.query("month == '00'").set_index("alpha_code")
    baseline = base.dropna(subset=markers)
    history = read_sheet("mh.xlsx", ["alpha_code", "mhsjog"]).set_index("alpha_code")
    assert history.index.is_unique
    assert base.index.isin(history.index).all()
    history_flag = history.mhsjog.reindex(base.index).eq(3)

    # Bunya et al., 14th International Symposium on Sjogren's, P-5, Table I.
    # The poster separates 375 without disease history, 66 with other autoimmune
    # disease and 53 with Sjogren history. Here the first two strata are pooled.
    # Values are: SJO records, complete sera, any candidate, SP1, CA6, PSP,
    # SSA >25, SSB >25, and the distribution for 0/1/2/3 positive antigens.
    reference = {
        "all": (494, 492, 167, 99, 86, 56, 57, 18, 325, 102, 56, 9),
        "history": (53, 52, 22, 11, 9, 5, 31, 13, 30, 19, 3, 0),
        "no_history": (441, 440, 145, 88, 77, 51, 26, 5, 295, 83, 53, 9),
    }
    baseline_rows = []
    for label, mask in (
        ("all", pd.Series(True, index=base.index)),
        ("history", history_flag),
        ("no_history", ~history_flag),
    ):
        sub = baseline.loc[mask, markers]
        for rule, compare in (("ge20", sub.ge(20)), ("gt20", sub.gt(20)), ("gt25", sub.gt(25))):
            antigen = pd.DataFrame({name: compare[fields].any(axis=1) for name, fields in groups.items()})
            number = antigen.sum(axis=1)
            row = {
                "stratum": label, "candidate_rule": rule, "lab_records": int(mask.sum()), "complete_n": len(sub),
                "any_candidate": int(number.gt(0).sum()),
                **{k: int(v) for k, v in antigen.sum().items()},
                "ssa_gt25": int(sub.ss_a_ro.gt(25).sum()), "ssb_gt25": int(sub.ss_b_la.gt(25).sum()),
                **{f"antigens_{i}": int(number.eq(i).sum()) for i in range(4)},
            }
            row["matches_2018_poster"] = tuple(list(row.values())[2:]) == reference[label]
            assert row["matches_2018_poster"] == (rule == "ge20" or (label == "history" and rule == "gt20"))
            baseline_rows.append(row)
    tables["baseline_2018_poster"] = pd.DataFrame(baseline_rows)

    def transition(a, b):
        pp, pn = int((a & b).sum()), int((a & ~b).sum())
        np_, nn = int((~a & b).sum()), int((~a & ~b).sum())
        lo, hi = proportion_confint(pp, pp + pn, method="wilson")
        return {"positive_at_baseline": pp + pn, "persistent_positive": pp,
                "positive_to_nonpositive": pn, "nonpositive_to_positive": np_,
                "persistent_nonpositive": nn, "persistence": pp / (pp + pn),
                "persistence_ci_low": lo, "persistence_ci_high": hi}

    rows = []
    for rule, baseline_flag, followup_flag in (("gt25", m0.gt(25), m12.gt(25)), ("ge20", m0.ge(20), m12.ge(20))):
        candidate_fields = [c for fields in groups.values() for c in fields]
        for assay in candidate_fields:
            rows.append({"rule": rule, "assay": assay, **transition(baseline_flag[assay], followup_flag[assay])})
        for isotype in ("igg", "iga", "igm"):
            fields = [c for c in candidate_fields if not c.endswith("_" + isotype)]
            rows.append({"rule": rule, "assay": f"panel_excluding_{isotype}", **transition(
                baseline_flag[fields].any(axis=1), followup_flag[fields].any(axis=1))})
    tables["paired_antibody_isotype_audit"] = pd.DataFrame(rows)

    candidate_fields = [c for fields in groups.values() for c in fields]
    max0, max12 = m0[candidate_fields].max(axis=1), m12[candidate_fields].max(axis=1)
    strength = []
    for band, inside in (("gt25_to35", max0.gt(25) & max0.le(35)),
                         ("gt35_to50", max0.gt(35) & max0.le(50)),
                         ("gt50", max0.gt(50))):
        strength.append({"baseline_max_band": band, "baseline_n": int(inside.sum()),
                         "followup_gt25": int((inside & max12.gt(25)).sum()),
                         "followup_20_to25": int((inside & max12.between(20, 25)).sum()),
                         "followup_lt20": int((inside & max12.lt(20)).sum()),
                         "any_same_candidate_isotype_gt25": int((inside & (
                             m0[candidate_fields].gt(25) & m12[candidate_fields].gt(25)).any(axis=1)).sum())})
    tables["paired_antibody_baseline_strength"] = pd.DataFrame(strength)

    results = {}
    positive0 = m0[novel].gt(25)
    positive12 = m12[novel].gt(25)
    ssa0, ssa12 = m0.ss_a_ro.gt(25), m12.ss_a_ro.gt(25)
    panel0, panel12 = positive0.any(axis=1), positive12.any(axis=1)
    rand = np.random.default_rng(20261108)
    idx = rand.integers(0, len(m0), size=(5000, len(m0)))
    def boot_fraction(base, follow):
        return ((base & follow).to_numpy()[idx].sum(axis=1) /
                base.to_numpy()[idx].sum(axis=1))
    delta = boot_fraction(panel0, panel12) - boot_fraction(ssa0, ssa12)
    results["panel_minus_ssa_persistence"] = {
        "estimate": float(panel12[panel0].mean() - ssa12[ssa0].mean()),
        "participant_bootstrap_95_ci": np.quantile(delta, [0.025, 0.975]).tolist(),
    }
    for name, values0, values12 in (
        ("novel_panel", m0[novel].max(axis=1), m12[novel].max(axis=1)),
        ("anti_ssa", m0.ss_a_ro, m12.ss_a_ro),
    ):
        base = values0.gt(25)
        follow = values12.gt(25)
        off = int((base & ~follow).sum())
        on = int((~base & follow).sum())
        ranges = {}
        for label, select in (
            ("25_to_35", values0.gt(25) & values0.le(35)),
            ("35_to_50", values0.gt(35) & values0.le(50)),
            ("above_50", values0.gt(50)),
        ):
            n = int(select.sum())
            yes = int((select & follow).sum())
            ranges[label] = {"baseline_positive": n, "still_positive": yes}
        results[name] = {
            "baseline_positive": int(base.sum()), "month12_positive": int(follow.sum()),
            "pos_to_nonpos": off, "nonpos_to_pos": on,
            "positive_to_borderline": int((base & values12.between(20, 25)).sum()),
            "positive_to_negative": int((base & values12.lt(20)).sum()),
            "negative_to_positive": int((values0.lt(20) & follow).sum()),
            "borderline_to_positive": int((values0.between(20, 25) & follow).sum()),
            "initial_level_persistence": ranges,
        }

    baseline = raw.query("month == '00'").dropna(subset=markers).set_index("alpha_code")
    retained = baseline.index.isin(m0.index)
    results["baseline_followup_availability"] = {}
    for label, has_marker, persistent in (
        ("novel_panel_positive", baseline[novel].gt(25).any(axis=1), panel0 & panel12),
        ("anti_ssa_positive", baseline.ss_a_ro.gt(25), ssa0 & ssa12),
    ):
        n_positive = int(has_marker.sum())
        n_missing = int((has_marker & ~retained).sum())
        n_persistent = int(persistent.sum())
        results["baseline_followup_availability"][label] = {
            "baseline_positive": n_positive,
            "missing_paired_test_among_baseline_positive": n_missing,
            "extreme_missing_status_bounds_for_persistence": [
                n_persistent / n_positive, (n_persistent + n_missing) / n_positive
            ],
            "paired_positive_at_baseline": int((has_marker & retained).sum()),
            "baseline_nonpositive": int((~has_marker).sum()),
            "paired_nonpositive_at_baseline": int((~has_marker & retained).sum()),
        }
    clinical = read_sheet("gi.xlsx", ["alpha_code","group"]).set_index("alpha_code")
    arm = clinical.group.reindex(m0.index)
    assert arm.notna().all()
    results["persistence_by_randomized_arm_descriptive"] = {
        str(k): {"n": int(arm.eq(k).sum()),
                 "novel_baseline_pos": int((arm.eq(k) & panel0).sum()),
                 "novel_persistent_pos": int((arm.eq(k) & panel0 & panel12).sum()),
                 "novel_incident_pos": int((arm.eq(k) & ~panel0 & panel12).sum())}
        for k in sorted(arm.unique())
    }
    results["per_marker_direction_of_raw_numeric_values"] = {
        col: {"increased": int(m12[col].gt(m0[col]).sum()),
              "decreased": int(m12[col].lt(m0[col]).sum()),
              "equal_or_same_censor_bound": int(m12[col].eq(m0[col]).sum())}
        for col in novel
    }

    date_fields = ["collection_date", "received_date", "srbldtc", "srblshdt"]
    dates = read_sheet("sjo-results.xlsx", ["alpha_code", "month", *date_fields])
    for col in date_fields:
        dates[col] = pd.to_numeric(dates[col], errors="raise")
    assert dates.collection_date.eq(dates.srbldtc).all()
    pair_dates = dates.pivot(index="alpha_code", columns="month", values=date_fields).reindex(m0.index)
    days = pair_dates[("collection_date", "12")] - pair_dates[("collection_date", "00")]
    assert days.notna().all() and days.between(240, 510).all()
    within_window = days.between(335, 395)
    results["collection_date_pair_check"] = {
        "basis": "difference of relative collection day codes, corroborated by srbldtc",
        "n_usable": int(days.notna().sum()),
        "median_days": float(days.median()),
        "iqr_days": [float(days.quantile(0.25)), float(days.quantile(0.75))],
        "min_days": int(days.min()),
        "max_days": int(days.max()),
        "n_335_to_395_days": int(within_window.sum()),
        "novel_baseline_pos_in_window": int((within_window & panel0).sum()),
        "novel_persistent_pos_in_window": int((within_window & panel0 & panel12).sum()),
        "ssa_baseline_pos_in_window": int((within_window & ssa0).sum()),
        "ssa_persistent_pos_in_window": int((within_window & ssa0 & ssa12).sum()),
    }
    results["specimen_collection_to_receipt_days"] = {}
    for visit in ("00", "12"):
        elapsed = pair_dates[("received_date", visit)] - pair_dates[("collection_date", visit)]
        results["specimen_collection_to_receipt_days"][visit] = {
            "median": float(elapsed.median()),
            "iqr": [float(elapsed.quantile(0.25)), float(elapsed.quantile(0.75))],
            "n_negative_or_missing": int((elapsed.lt(0) | elapsed.isna()).sum()),
        }

    # Independent directional check of one published prognostic association.
    # CA6 IgG >=20 and >25 are not used to select outcome or covariates here.
    # The endpoint is the 12-month two-eye mean of the anesthetized Schirmer result.
    # Baseline Schirmer, age, sex, randomized arm and traditional SSA are covariates.
    demo = read_sheet("gi.xlsx", ["alpha_code","group","age","sex"]).set_index("alpha_code")
    e0 = read_sheet("eo.xlsx",["alpha_code","month","mmwetrt","mmwetle"]).query("month == '00'").set_index("alpha_code")
    e12 = read_sheet("oa.xlsx",["alpha_code","month","mmwetrt","mmwetle"]).query("month == '12'").set_index("alpha_code")
    for tbl in (e0,e12):
        tbl["schirmer"] = tbl[["mmwetrt","mmwetle"]].apply(pd.to_numeric, errors="coerce").mean(axis=1,skipna=False)
    m = baseline[["ca_vi_igg","ss_a_ro"]].join(demo).join(
        e0.schirmer.rename("schirmer_0")).join(e12.schirmer.rename("schirmer_12")).dropna()
    x = pd.DataFrame({"schirmer_0":m.schirmer_0,"age":m.age,"female":m.sex.eq("F").astype(float),
                      "group":m.group,"ssa":m.ss_a_ro.gt(25).astype(float),
                      "ca6_igg_log2":np.log2(m.ca_vi_igg)})
    x = x.astype(float)
    assert np.isfinite(x.to_numpy()).all()
    fit = sm.OLS(m.schirmer_12.astype(float), sm.add_constant(x)).fit(cov_type="HC3", use_t=True)
    results["ca6_directional_prognostic_check"] = {
        "n":len(m),"coefficient_mm_per_ca6_igg_doubling":float(fit.params.iloc[-1]),
        "ci95":fit.conf_int().iloc[-1].tolist(),
        "both_eyes_required":True,
        "below_limit_numeric_proxy":1.0,"above_limit_numeric_proxy":160.0,
    }

    tables["supporting_checks"] = flatten_results(results)

    key = "alpha_code"
    def aligned(name, fields, visit=None):
        frame = data.read(name, fields, visit)
        return frame.reindex(m0.index)
    conditions = {
        "reported_sjogren": "mhsjog", "reported_rheumatoid_arthritis": "mhrheu",
        "reported_thyroid_dysfunction": "mhthyr", "reported_diabetes": "mhdiab",
        "reported_hypertension": "mhhypert", "reported_depression": "mhdepc",
    }
    history = aligned("mh.xlsx", [*conditions.values(), "cmcort", "cmcortnow"])
    assert history[list(conditions.values())].isin([1, 2, 3]).all().all()
    features = pd.DataFrame({name: history[field].eq(3) for name, field in conditions.items()})
    assert not (history.cmcort.eq(0) & history.cmcortnow.eq(1)).any()
    steroid = history.cmcortnow.where(history.cmcortnow.notna(), history.cmcort.where(history.cmcort.eq(0)))
    features["current_corticosteroids_or_immunosuppressants"] = steroid.eq(1).where(steroid.notna())
    work = aligned("wp.xlsx", ["opwork", "opwprod", "opacth"], "00")
    assert work.opwork.isin([0, 1]).all()
    features["employed_for_pay"] = work.opwork.eq(1)
    smoke = aligned("br.xlsx", ["sucigs", "susmcur"], "00")
    assert smoke.sucigs.isin([0, 1]).all()
    features["current_cigarette_smoking"] = smoke.susmcur.eq(1).where(smoke.susmcur.notna() | smoke.sucigs.eq(0))
    lab0 = aligned("sjo-results.xlsx", ["ana_titer"], "00")
    ana = lab0.ana_titer.str.extract(r":(\d+)", expand=False).astype(float)
    assert (ana.notna() | lab0.ana_titer.eq("-")).all()
    features["ana_titer_ge_1_320"] = ana.ge(320)
    features["anti_ssa_gt25"] = m0.ss_a_ro.gt(25)
    features["anti_ssb_gt25"] = m0.ss_b_la.gt(25)
    features["both_ssa_ssb_nonpositive"] = ~(features["anti_ssa_gt25"] | features["anti_ssb_gt25"])
    work_scores = {
        "work_impairment_0_to_10_if_employed": work.opwprod.where(work.opwork.eq(1)),
        "regular_activity_impairment_0_to_10": work.opacth,
    }
    feature_rows, difference_rows, work_rows, target_rows = [], [], [], []
    for cut, first, second in (("gt25", m0[novel].gt(25), m12[novel].gt(25)),
                               ("ge20", m0[novel].ge(20), m12[novel].ge(20))):
        base, follow = first.any(axis=1), second.any(axis=1)
        retain = (first & second).any(axis=1)
        groups = {
            "same_positive_item": retain,
            "different_item_only": base & follow & ~retain,
            "no_longer_panel_positive": base & ~follow,
            "new_panel_positive": ~base & follow,
            "both_nonpositive": ~base & ~follow,
        }
        for group, select in groups.items():
            for name in features:
                v = features.loc[select, name].dropna().astype(bool)
                feature_rows.append({
                    "rule": cut, "group": group, "feature": name, "group_n": int(select.sum()),
                    "available": len(v), "positive": int(v.sum()), "proportion": v.mean(),
                })
            for name, values in work_scores.items():
                v = values.loc[select].dropna()
                work_rows.append({"rule": cut, "group": group, "feature": name, "available": len(v),
                                  "median": v.median(), "q1": v.quantile(.25), "q3": v.quantile(.75)})
        for name in features:
            v = features[name].dropna().astype(bool)
            for endpoint, eligible, event in (("same_item_retention", base, retain),
                                               ("new_panel_positive", ~base, follow)):
                mask = eligible.reindex(v.index)
                yes, no = v & mask, ~v & mask
                denom_yes, denom_no = int(yes.sum()), int(no.sum())
                endpoint_status = event.reindex(v.index)
                n_yes = int((yes & endpoint_status).sum())
                n_no = int((no & endpoint_status).sum())
                if denom_yes and denom_no:
                    lo, hi = confint_proportions_2indep(n_yes, denom_yes, n_no, denom_no, method="newcomb")
                    difference_rows.append({"rule": cut, "feature": name, "endpoint": endpoint,
                                            "events_with_feature": n_yes, "n_with_feature": denom_yes,
                                            "events_without_feature": n_no, "n_without_feature": denom_no,
                                            "unadjusted_risk_difference": n_yes / denom_yes - n_no / denom_no,
                                            "ci_low": lo, "ci_high": hi})
        sub = features.both_ssa_ssb_nonpositive
        target_rows.append({"rule": cut, "n_ssa_ssb_nonpositive": int(sub.sum()),
                            "candidate_positive_at_baseline": int((sub & base).sum()),
                            "panel_persistent": int((sub & base & follow).sum()),
                            "same_item_persistent": int((sub & retain).sum()),
                            "new_panel_positive": int((sub & ~base & follow).sum())})

    # Exact numeric values can be compared directly. For left/right censoring, use
    # interval bounds (0,1) and (160,infinity), never treating bounds as measured values.
    raw = read_sheet("sjo-results.xlsx", [key, "month", "ss_a_ro", "ss_b_la", *novel])
    raw = raw.loc[raw[key].isin(m0.index) & raw.month.isin(["00", "12"])]
    assert not raw.duplicated([key, "month"]).any()
    pair = raw.pivot(index=key, columns="month").reindex(m0.index)

    def bounds(series):
        text = series.astype(str).str.strip()
        lo_flag, hi_flag = text.eq("<1.0"), text.eq(">160")
        v = pd.to_numeric(text.str.removeprefix("<").str.removeprefix(">"), errors="raise")
        return v.mask(lo_flag, 0.0), v.mask(hi_flag, np.inf), ~(lo_flag | hi_flag), v

    assay_rows = []
    for assay in [*novel, "ss_a_ro", "ss_b_la"]:
        lo0, hi0, exact0, proxy0 = bounds(pair[(assay, "00")])
        lo1, hi1, exact1, proxy1 = bounds(pair[(assay, "12")])
        delta_lo, delta_hi = lo1 - hi0, hi1 - lo0
        # Endpoints on the open censoring limits also define direction, e.g. <1 to 1.
        direction_up = delta_lo.gt(0) | (delta_lo.eq(0) & ~(exact0 & exact1))
        direction_down = delta_hi.lt(0) | (delta_hi.eq(0) & ~(exact0 & exact1))
        exact_tie = exact0 & exact1 & delta_lo.eq(0)
        known = direction_up | direction_down | exact_tie
        assert not (direction_up & direction_down).any()
        exact = exact0 & exact1
        inc = proxy0.le(25) & proxy1.gt(25)
        dec = proxy0.gt(25) & proxy1.le(25)
        novel_no_baseline = ~m0[novel].gt(25).any(axis=1)
        assay_rows.append({"assay": assay, "n": len(exact), "exact_pairs": int(exact.sum()),
                           "baseline_median_lower": lo0.median(), "baseline_median_upper": hi0.median(),
                           "month12_median_lower": lo1.median(), "month12_median_upper": hi1.median(),
                           "paired_median_change_lower": delta_lo.median(),
                           "paired_median_change_upper": delta_hi.median(),
                           "median_change_exact_pairs": (proxy1 - proxy0)[exact].median(),
                           "increased": int(direction_up.sum()), "decreased": int(direction_down.sum()),
                           "ties_exact": int(exact_tie.sum()), "unknown_direction": int((~known).sum()),
                           "new_strict_positive": int(inc.sum()), "no_longer_strict_positive": int(dec.sum()),
                           "new_strict_positive_in_initial_candidate_panel_nonpositive": int((inc & novel_no_baseline).sum())})
    # Keep the same 307 initially nonpositive participants when excluding assays.
    base_strict = m0[novel].gt(25).any(axis=1)
    month12_strict = m12[novel].gt(25)
    max0 = m0[novel].astype(float).max(axis=1)
    igg_iga = [assay for assay in novel if not assay.endswith("igm")]
    two_igm = ["sp_1_igm", "ca_vi_igm"]
    gain_context_rows = []
    for band, inside in (
        ("all_initially_nonpositive", ~base_strict),
        ("max_below_15", ~base_strict & max0.lt(15)),
        ("max_15_to_below_20", ~base_strict & max0.ge(15) & max0.lt(20)),
        ("max_20_to_25", ~base_strict & max0.between(20, 25)),
    ):
        any_gain = month12_strict.any(axis=1)
        has_two_igm = month12_strict[two_igm].any(axis=1)
        has_other = month12_strict.drop(columns=two_igm).any(axis=1)
        gain_context_rows.append({
            "band": band, "n": int(inside.sum()), "any_new_panel_positive": int((inside & any_gain).sum()),
            "positive_on_at_least_one_of_two_igm": int((inside & has_two_igm).sum()),
            "positive_only_on_one_or_both_of_two_igm": int((inside & has_two_igm & ~has_other).sum()),
            "still_positive_excluding_two_igm": int((inside & has_other).sum()),
            "still_positive_excluding_all_igm": int((inside & month12_strict[igg_iga].any(axis=1)).sum()),
        })
    tables["clinical_context_new_positivity_composition"] = pd.DataFrame(gain_context_rows)

    for name, rows in (("clinical_context_by_trajectory", feature_rows), ("clinical_context_unadjusted", difference_rows),
                       ("employment_activity_by_trajectory", work_rows), ("clinical_context_ssa_ssb_nonpositive", target_rows),
                       ("paired_raw_assay_changes", assay_rows)):
        tables[name] = pd.DataFrame(rows)

    groups = ANTIGENS
    def status(df, novel_cut):
        return pd.concat([df[["ss_a_ro", "ss_b_la"]].gt(25),
                          df[novel].gt(25) if novel_cut == "strict" else df[novel].ge(20)], axis=1)

    def summary(rule, marker, first, second):
        repeat = int((first & second).sum())
        loss = int((first & ~second).sum())
        gain = int((~first & second).sum())
        neg = int((~first & ~second).sum())
        lo, hi = proportion_confint(repeat, repeat + loss, method="wilson")
        return dict(rule=rule, marker=marker, n=len(first), repeat_pos=repeat, pos_to_nonpos=loss,
                    nonpos_to_pos=gain, repeat_nonpos=neg, persistence=repeat/(repeat+loss),
                    persistence_ci_low=lo, persistence_ci_high=hi)

    rows = []
    for rule in ("strict", "incl_borderline"):
        a, b = status(m0, rule), status(m12, rule)
        for marker, fields in {"anti_SSA": ["ss_a_ro"], "anti_SSB": ["ss_b_la"],
                               **groups, "novel_any": novel}.items():
            rows.append(summary(rule, marker, a[fields].any(axis=1), b[fields].any(axis=1)))
        new_a, new_b = a[novel], b[novel]
        eligible = ~a[["ss_a_ro", "ss_b_la"]].any(axis=1)
        rows.append(summary(rule, "novel_any_in_baseline_SSA_SSB_nonpositive",
                            new_a.loc[eligible].any(axis=1), new_b.loc[eligible].any(axis=1)))
        eligible = a.ss_a_ro & new_a.any(axis=1)
        rows.append(summary(rule, "novel_any_among_copositive_with_SSA",
                            new_a.loc[eligible].any(axis=1), new_b.loc[eligible].any(axis=1)))
        rows.append(summary(rule, "anti_SSA_among_copositive_with_novel",
                            a.loc[eligible, "ss_a_ro"], b.loc[eligible, "ss_a_ro"]))
        base_pos = new_a.any(axis=1)
        rows.append(summary(rule, "novel_same_isotype_in_baseline_panel_positive",
                            base_pos.loc[base_pos], (new_a & new_b).any(axis=1).loc[base_pos]))
    out = pd.DataFrame(rows)

    tables["paired_antibody_pilot"] = out

    return tables

if __name__ == "__main__":
    run_cli(compute, "STEP04_supporting_results.xlsx", __doc__)
