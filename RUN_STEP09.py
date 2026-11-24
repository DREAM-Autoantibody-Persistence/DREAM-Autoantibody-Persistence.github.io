"""Step 09: quantitative longitudinal tracking and month-24 feasibility.

Uses the frozen official workbooks. Exact-value paired analyses exclude either
visit's <1.0 or >160 results; bound-proxy rank correlations are sensitivity
analyses only. These are longitudinal associations, not assay repeatability.
No participant-level results are exported. Existing Steps 01-08 are untouched.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest, rankdata
from statsmodels.stats.multitest import multipletests
from dream_common import ASSAYS, NOVEL, KEY, fraction, run_cli, require

SEED = 20261123
B = 2000


def rank_boot(x, y, rng):
    n = len(x)
    if n < 10 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return np.nan, np.nan
    ix = rng.integers(0, n, size=(B, n))
    a = rankdata(x[ix], axis=1)
    b = rankdata(y[ix], axis=1)
    a -= a.mean(axis=1, keepdims=True)
    b -= b.mean(axis=1, keepdims=True)
    den = np.sqrt((a*a).sum(axis=1)*(b*b).sum(axis=1))
    z = np.divide((a*b).sum(axis=1), den, out=np.full(B, np.nan), where=den > 0)
    return np.nanquantile(z, [.025, .975])


def compute(data):
    m0, m12 = data.pairs()
    raw = data.read_sheet('sjo-results.xlsx')
    r0 = raw.loc[raw.month.eq('00')].set_index(KEY).loc[m0.index]
    r12 = raw.loc[raw.month.eq('12')].set_index(KEY).loc[m0.index]
    rng = np.random.default_rng(SEED)
    quantitative, directions, cross, extflow, exttrans, trajectories, checks = [], [], [], [], [], [], []
    for a in ASSAYS:
        token0 = r0[a].astype(str).str.strip().str.upper()
        token12 = r12[a].astype(str).str.strip().str.upper()
        exact0 = pd.to_numeric(token0, errors='coerce')
        exact12 = pd.to_numeric(token12, errors='coerce')
        ok = exact0.notna() & exact12.notna() & exact0.gt(0) & exact12.gt(0)
        x, y = exact0[ok].to_numpy(), exact12[ok].to_numpy()
        rho, p = spearmanr(x, y)
        independent_rho=pd.Series(x).rank().corr(pd.Series(y).rank())
        require(np.isclose(rho,independent_rho),'Rank correlation reconstruction differs: '+a)
        checks.append(dict(check='Spearman_independent_rank_Pearson_'+a,passed=True))
        lo, hi = rank_boot(x, y, rng)
        delta = y-x
        logratio = np.log2(y/x)
        draws = rng.integers(0, len(x), size=(B, len(x)))
        med_ci = np.quantile(np.median(delta[draws], axis=1), [.025, .975])
        quantitative.append(dict(assay=a, exact_pair_n=int(ok.sum()),
            excluded_bound_pair_n=int((~ok).sum()),
            baseline_exact_median=float(np.median(x)), month12_exact_median=float(np.median(y)),
            paired_delta_median=float(np.median(delta)), delta_median_ci_low=med_ci[0], delta_median_ci_high=med_ci[1],
            paired_log2ratio_median=float(np.median(logratio)), spearman_rho=rho,
            rho_boot_ci_low=lo, rho_boot_ci_high=hi, p_spearman=p,
            bound_proxy_rho_sensitivity=spearmanr(m0[a],m12[a]).statistic,
            exact_at_least_twofold_up=int((logratio>=1).sum()),
            exact_at_least_twofold_down=int((logratio<=-1).sum())))
        up, down, tied = (delta>0).sum(), (delta<0).sum(), (delta==0).sum()
        directions.append(dict(assay=a, exact_pair_n=len(x), increased=int(up),decreased=int(down), equal=int(tied),
            p_sign=binomtest(int(up),int(up+down),.5).pvalue if up+down else 1.0))
        for label, b0, b12 in [('gt25',m0[a]>25,m12[a]>25),('ge20',m0[a]>=20,m12[a]>=20)]:
            # Independently classify raw strings, including open censor bounds.
            def rawflag(tokens,cut,strict):
                out=pd.to_numeric(tokens,errors='coerce')
                out=(out>cut) if strict else (out>=cut)
                return out.mask(tokens.eq('>160'),True).mask(tokens.eq('<1.0'),False)
            require(rawflag(token0,25 if label=='gt25' else 20,label=='gt25').equals(b0), 'Raw baseline classification differs')
            require(rawflag(token12,25 if label=='gt25' else 20,label=='gt25').equals(b12), 'Raw followup classification differs')
            for state, mask in [('positive_to_nonpositive',b0&~b12),('nonpositive_to_positive',~b0&b12)]:
                sel = mask & ok
                ld = np.log2(exact12[sel]/exact0[sel])
                clear = ((m0[a]>25)&(m12[a]<20)) if state.startswith('positive') else ((m0[a]<20)&(m12[a]>25))
                cross.append(dict(rule=label,assay=a,transition=state,total_n=int(mask.sum()),exact_n=int(sel.sum()),
                    clear_lt20_to_gt25_or_reverse=int((mask&clear).sum()),
                    exact_median_log2ratio=ld.median(),exact_at_least_twofold=int((ld.abs()>=1).sum())))
    qt, dt = pd.DataFrame(quantitative), pd.DataFrame(directions)
    for tab, pcol, new in [(qt,'p_spearman','p_spearman_BH9'),(dt,'p_sign','p_sign_BH9')]:
        tab[new]=np.nan
        novel=tab.assay.isin(NOVEL)
        tab.loc[novel,new]=multipletests(tab.loc[novel,pcol],method='fdr_bh')[1]

    # Reconcile three-time-point availability, without treating extension as a
    # representative cohort or as independent external validation.
    labs=data.numeric_lab()
    complete=labs.dropna(subset=ASSAYS)
    visits={v:complete.loc[complete.month.eq(v)].set_index(KEY)[ASSAYS] for v in ['00','12','24']}
    ids=sorted(set(visits['00'].index)&set(visits['12'].index)&set(visits['24'].index))
    gi=data.read_sheet('gi.xlsx').set_index(KEY)
    for v in ['00','12','24']:
        extflow.append(dict(item='visit_'+v+'_records',n=int(raw.month.eq(v).sum())))
        extflow.append(dict(item='visit_'+v+'_complete_11',n=len(visits[v])))
    extflow.append(dict(item='complete_three_visits',n=len(ids)))
    ids1224=sorted(set(visits['12'].index)&set(visits['24'].index))
    extflow.append(dict(item='complete_month12_month24_pairs',n=len(ids1224)))
    for value,n in gi.loc[ids,'group'].value_counts(dropna=False).items():
        extflow.append(dict(item='three_visit_original_group_'+str(value),n=int(n)))
    if ids:
        tt={v:visits[v].loc[ids] for v in visits}
        for rule,threshold in [('gt25',25),('ge20',20)]:
            flag={v:(tt[v][NOVEL]>threshold if rule=='gt25' else tt[v][NOVEL]>=threshold) for v in tt}
            panel={v:f.any(axis=1) for v,f in flag.items()}
            codes=panel['00'].astype(int).astype(str)+panel['12'].astype(int).astype(str)+panel['24'].astype(int).astype(str)
            for pattern,n in codes.value_counts().sort_index().items():
                trajectories.append(dict(rule=rule,pattern_0_12_24='M0-'+pattern[0]+'_M12-'+pattern[1]+'_M24-'+pattern[2],n=int(n)))
            for start,end in [('00','12'),('12','24'),('00','24')]:
                base=panel[start]; follow=panel[end]
                same=(flag[start]&flag[end]).any(axis=1)
                for label,ev in [('panel',base&follow),('same_assay',same)]:
                    exttrans.append(dict(scope='three_timepoint_complete',rule=rule,start=start,end=end,measure=label,
                        **fraction(int(ev.sum()),int(base.sum()))))
            b=panel['00']
            allsame=(flag['00']&flag['12']&flag['24']).any(axis=1)
            exttrans.append(dict(scope='three_timepoint_complete',rule=rule,start='00',end='24',measure='identical_assay_at_all_three',
                **fraction(int(allsame.sum()),int(b.sum()))))
    if ids1224:
        for rule in ['gt25','ge20']:
            a=visits['12'].loc[ids1224,NOVEL]
            b=visits['24'].loc[ids1224,NOVEL]
            a=(a>25) if rule=='gt25' else (a>=20)
            b=(b>25) if rule=='gt25' else (b>=20)
            base,follow=a.any(axis=1),b.any(axis=1)
            for label,event in [('panel',base&follow),('same_assay',(a&b).any(axis=1))]:
                exttrans.append(dict(scope='all_complete_12_24_pairs',rule=rule,start='12',end='24',measure=label,
                    **fraction(int(event.sum()),int(base.sum()))))
            exttrans.append(dict(scope='all_complete_12_24_pairs',rule=rule,start='12',end='24',measure='newly_threshold_positive',
                    **fraction(int((~base&follow).sum()),int((~base).sum()))))
        for value,n in gi.loc[ids1224,'group_extension'].value_counts(dropna=False).items():
            extflow.append(dict(item='12_24_pair_extension_group_'+str(value),n=int(n)))
    # Both ends of a clear crossing are outside the borderline band.
    max0,max12=m0[NOVEL].max(axis=1),m12[NOVEL].max(axis=1)
    panel_clear=pd.DataFrame([
        dict(transition='baseline_gt25_to_month12_all_lt20',**fraction(int(((max0>25)&(max12<20)).sum()),int((max0>25).sum()))),
        dict(transition='baseline_all_lt20_to_month12_gt25',**fraction(int(((max0<20)&(max12>25)).sum()),int((max0<20).sum())))])
    methods=pd.DataFrame([
        ('cohort','407 original complete baseline/month12 pairs; verify against Step02'),
        ('quantitative','Each assay separately; primary quantitative summaries use exact positive numeric values at both visits'),
        ('censoring','<1.0 and >160 are not exact values; excluded from primary quantitative analyses; bound-proxy ranks sensitivity only'),
        ('correlation','Spearman rank tracking across visits, not test-retest reliability or diagnostic accuracy'),
        ('uncertainty',f'Participant-pair percentile bootstrap; {B} draws; seed {SEED}; nominal 95% intervals'),
        ('multiplicity','BH across nine candidate assays separately for correlations and exact-pair sign tests; SSA/SSB descriptive comparators'),
        ('twofold','Exploratory magnitude descriptor, not a validated analytical or clinical change threshold'),
        ('month24','Selected extension participants; original active group eligible and extension re-randomized; not external validation'),
        ('intervals','Month24 selection and re-randomization prevent direct natural-history interpretation'),
        ('clear_crossing','Outside 20-25 borderline band at both visits; does not exclude assay or batch variability'),
    ],columns=['item','definition'])
    sources=pd.DataFrame([
        ('DREAM','Official frozen source workbooks; 407 pairs','Direct longitudinal analysis; checksums in source_checks','https://www.med.upenn.edu/cpob/dream-data-export'),
        ('DREAM extension','Month24 complete triple subset n='+str(len(ids)), 'Selected same-cohort extension; original active-group enrollment and re-randomization; not independent validation','https://pmc.ncbi.nlm.nih.gov/articles/PMC7004875/'),
        ('GSE272409','Public minor salivary gland scRNA-seq: Sjogren/sicca','Can study antigen transcript cell localization; no matched SP1/CA6/PSP antibody trajectories verified; expression matrices not downloaded or analyzed','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272409'),
        ('GSE157278','Public PBMC scRNA-seq: Sjogren/healthy','Blood-cell expression; no antibody-specific repeat outcomes verified; not a validation of antibody persistence','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157278'),
        ('CA6 pQTL','Carbonic anhydrase 6 protein pQTL; rs3765967; discovery n=996','Instrument concerns CA6 protein abundance, not anti-CA6 antibodies or one-year persistence; no MR performed','https://metabolomics.helmholtz-munich.de/pgwas/index.php?loc=400&task=locuscard'),
        ('OpenGWAS CA6','prot-c-3352_80_3 catalog entry','Catalog located; automated download requires JWT; no authenticated summary statistics downloaded','https://opengwas.io/datasets/prot-c-3352_80_3'),
        ('Raolji 2026','DED candidate autoantibody study; data on reasonable request','Not openly downloadable repeated-antibody data; author request excluded by current project constraints','https://doi.org/10.1136/bmjophth-2026-002776'),
        ('SP1 mapping','Human SP1 is Sp1 transcription factor','Do not use human SP1 gene expression or pQTL as a surrogate for salivary protein-1 autoantibody','https://www.ncbi.nlm.nih.gov/gene/6667'),
        ('PSP mapping','PSP synonym maps to BPIFA2','Protein/transcript abundance still differs from anti-PSP antibody concentration and persistence','https://www.ncbi.nlm.nih.gov/datasets/gene/140683/'),
    ],columns=['resource','verified_information','fitness_for_this_question','source_url'])
    sources['checked_date']='2026-11-23'
    checks.append(dict(check='raw_token_threshold_reconstruction_all_11_at_both_visits',passed=True))
    checks.append(dict(check='original_pair_n_407',passed=len(m0)==407))
    for rule in ['gt25','ge20']:
        require(sum(r['n'] for r in trajectories if r['rule']==rule)==len(ids),'Month24 trajectory counts do not reconcile')
    checks.append(dict(check='three_timepoint_trajectories_reconcile',passed=True))
    require(len(m0)==407,'Original cohort count differs')
    return dict(methods=methods,quantitative_tracking=qt,paired_direction=dt,
                assay_transition_magnitude=pd.DataFrame(cross),clear_panel_crossings=panel_clear,
                month24_feasibility=pd.DataFrame(extflow),month24_persistence=pd.DataFrame(exttrans),
                month24_trajectories=pd.DataFrame(trajectories),resource_assessment=sources,
                source_checks=data.check_sources(),analysis_checks=pd.DataFrame(checks))


if __name__=='__main__':
    run_cli(compute,'STEP09_quantitative_and_month24.xlsx',__doc__)
