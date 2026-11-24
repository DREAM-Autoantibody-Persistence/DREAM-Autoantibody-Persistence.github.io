"""Step 10: person-specific persistence benchmark and two targeted ocular checks.

Randomize whole month-12 assay vectors, not individual assay columns. The
randomization reference distribution is NOT a confidence interval for observed
persistence. This is an internal correspondence benchmark, not external or
technical validation. Prognostic models examine two literature-motivated
directions with fixed endpoints; they do not reproduce another study's exact
progression definition and do not establish clinical prediction performance.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from dream_common import ASSAYS, NOVEL, ANTIGENS, KEY, require, fraction, run_cli

B = 20000
SEED = 20261123


def bitcode(frame):
    return (frame.to_numpy(dtype=np.uint16) * (1 << np.arange(frame.shape[1],dtype=np.uint16))).sum(axis=1,dtype=np.uint16)


def compute(data):
    m0,m12=data.pairs()
    demo=data.read('gi.xlsx',['age','sex','group'])
    d=demo.loc[m0.index]
    require(d.notna().all().all(),'Missing permutation stratification variables')
    band=pd.cut(d.age,[-np.inf,50,65,np.inf],right=False,labels=['lt50','50to64','ge65'])
    strata=d.sex.astype(str)+'_'+d.group.astype(str)+'_'+band.astype(str)
    groups=[np.flatnonzero(strata.to_numpy()==v) for v in sorted(strata.unique())]
    blocks=pd.DataFrame([dict(stratum=str(strata.iloc[ix[0]]),n=len(ix)) for ix in groups])
    require(sum(blocks.n)==407,'Permutation blocks do not cover cohort')
    rows,methodrows,checks=[],[],[]
    for rule,cut in [('GT25',25),('GE20',20)]:
        a=(m0[NOVEL]>cut) if rule=='GT25' else (m0[NOVEL]>=cut)
        b=(m12[NOVEL]>cut) if rule=='GT25' else (m12[NOVEL]>=cut)
        ca,cb=bitcode(a),bitcode(b)
        ag_a=pd.DataFrame({k:a[v].any(axis=1) for k,v in ANTIGENS.items()})
        ag_b=pd.DataFrame({k:b[v].any(axis=1) for k,v in ANTIGENS.items()})
        aa,ab=bitcode(ag_a),bitcode(ag_b)
        pos0=ca>0; pos12=cb>0;den=int(pos0.sum())
        observed=np.array([np.count_nonzero(pos0 & pos12),np.count_nonzero((aa&ab)>0),np.count_nonzero((ca&cb)>0)])
        direct=np.array([int((a.any(axis=1)&b.any(axis=1)).sum()),int((ag_a&ag_b).any(axis=1).sum()),int((a&b).any(axis=1).sum())])
        require(np.array_equal(observed,direct),'Bit-coded counts differ from direct definitions')
        state=pos0.astype(int)*2+pos12.astype(int)
        transition_groups=[np.flatnonzero(state==v) for v in np.unique(state)]
        schemes={'unrestricted':[np.arange(len(m0))],
                 'sex_ageband_originalarm':groups,
                 'panel_transition_conditioned':transition_groups}
        for scheme_no,(scheme,gs) in enumerate(schemes.items()):
            rng=np.random.default_rng(SEED+(0 if rule=='GT25' else 100)+scheme_no)
            null=np.zeros((B,3),dtype=int)
            expected=np.zeros(3)
            for ix in gs:
                expected[0]+=np.sum(pos0[ix,None] & pos12[ix][None,:])/len(ix)
                expected[1]+=np.sum((aa[ix,None]&ab[ix][None,:])>0)/len(ix)
                expected[2]+=np.sum((ca[ix,None]&cb[ix][None,:])>0)/len(ix)
            for k in range(B):
                order=np.arange(len(m0))
                for ix in gs:
                    order[ix]=rng.permutation(ix)
                null[k]=[np.count_nonzero(pos0&pos12[order]),np.count_nonzero((aa&ab[order])>0),np.count_nonzero((ca&cb[order])>0)]
                if k==0:
                    require(np.array_equal(np.sort(cb[order]),np.sort(cb)),'Whole-panel vector multiset changed')
                    if scheme=='sex_ageband_originalarm':
                        require(np.array_equal(strata.to_numpy()[order],strata.to_numpy()),'Stratified permutation crossed blocks')
                    if scheme=='panel_transition_conditioned':
                        require(np.array_equal(state[order],state),'Transition-conditioned shuffle crossed status groups')
            for j,label in enumerate(['panel','same_antigen','same_assay']):
                mcse=float(null[:,j].std(ddof=1)/np.sqrt(B))
                require(abs(null[:,j].mean()-expected[j])<max(.08,6*mcse),'Monte Carlo mean differs from analytical expectation')
                low,high=np.quantile(null[:,j],[.025,.975])
                rows.append(dict(rule=rule,scheme=scheme,definition=label,
                    **fraction(int(observed[j]),den),null_expected_count=expected[j],
                    null_expected_percent=100*expected[j]/den,
                    observed_minus_null_pp=100*(observed[j]-expected[j])/den,
                    null_reference_low_count=low,null_reference_high_count=high,
                    null_reference_low_percent=100*low/den,null_reference_high_percent=100*high/den,
                    p_upper_monte_carlo=(1+int((null[:,j]>=observed[j]).sum()))/(B+1),draws=B,
                    null_mean_mcse_count=mcse))
            checks.append(dict(check=rule+'_'+scheme+'_whole_vectors_and_expected_counts',passed=True))
    bench=pd.DataFrame(rows)
    # Six comparisons per scheme (3 definitions x 2 thresholds); show all,
    # with same-assay GT25 as the central correspondence question.
    for scheme,sub in bench.groupby('scheme'):
        bench.loc[sub.index,'p_holm_six']=multipletests(sub.p_upper_monte_carlo,method='holm')[1]

    # Reassess two specified prognostic directions using all complete baseline
    # laboratory participants with outcome/covariate availability, not only
    # participants with repeat antibody testing.
    base=data.baseline()
    clinical=demo.join(data.eyes('00')).join(data.eyes('12'))
    clinical['female']=clinical.sex.eq('F').where(clinical.sex.notna()).astype(float)
    clinical['age10']=clinical.age/10
    raw0=data.read_sheet('sjo-results.xlsx',[KEY,'month','ca_vi_igg','psp_igm','ss_a_ro'])
    raw0=raw0.loc[raw0.month.eq('00')].set_index(KEY)
    mods=[]
    for assay,outcome,expected_direction in [('ca_vi_igg','schirmer','negative'),('psp_igm','osdi','positive')]:
        z=clinical.join(base[[assay,'ss_a_ro']],how='inner')
        z['log2_antibody']=np.log2(z[assay])
        z['ssa']=(z.ss_a_ro>25).astype(float)
        tokens=raw0[assay].reindex(z.index).astype(str).str.strip()
        exact=pd.to_numeric(tokens,errors='coerce').notna()
        for subset in ['all_with_reporting_bound_proxies','exact_exposure_values_only']:
            use=z if subset.startswith('all_') else z.loc[exact]
            for model in ['baseline_age_sex_arm','plus_ssa']:
                cov=[outcome+'_00','age10','female','group']
                if model=='plus_ssa': cov.append('ssa')
                cov.append('log2_antibody')
                obs=use[[outcome+'_12',*cov]].dropna()
                require(np.isfinite(obs.to_numpy(dtype=float)).all(),'Nonfinite prognostic model input')
                fit=sm.OLS(obs[outcome+'_12'].astype(float),sm.add_constant(obs[cov].astype(float))).fit(cov_type='HC3',use_t=True)
                lo,hi=fit.conf_int().loc['log2_antibody']
                mods.append(dict(assay=assay,outcome=outcome,analysis_set=subset,model=model,n=len(obs),
                    censored_exposure_n=int((~exact.reindex(obs.index)).sum()),
                    coefficient_per_doubling=fit.params['log2_antibody'],ci_low=lo,ci_high=hi,
                    p_nominal=fit.pvalues['log2_antibody'],literature_direction=expected_direction,
                    endpoint='month12 continuous outcome adjusted for baseline; not original binary progression endpoint',
                    covariance='HC3; t-based inference'))
    mods=pd.DataFrame(mods)
    for _,sub in mods.groupby(['analysis_set','model']):
        require(len(sub)==2,'Two fixed endpoints expected per specification')
        mods.loc[sub.index,'p_holm_two']=multipletests(sub.p_nominal,method='holm')[1]
    # Match the previously implemented CA6 check without running unrelated models.
    ca=mods.loc[mods.assay.eq('ca_vi_igg')&mods.analysis_set.str.startswith('all_')&mods.model.eq('plus_ssa')].iloc[0]
    checks.append(dict(check='CA6_original_Step04_model_definition_reconstructed_n451',passed=int(ca.n)==451))
    require(int(ca.n)==451,'CA6 prognostic model sample differs from original audit')
    methods=pd.DataFrame([
        ('source','Eight frozen official workbooks; participant-level data not exported'),
        ('benchmark','Observed within-person persistence vs random reassignment of complete month12 nine-assay vectors; permutations may retain some original pairings'),
        ('central_comparison','GT25 same-assay correspondence; other nested definitions and GE20 sensitivity reported completely; no prospective registration implied'),
        ('null_reference','Central 95% randomization reference range, not confidence interval for observed persistence'),
        ('permutation','20000 draws per scheme/threshold; seed20261123 plus documented offsets'),
        ('stratification','Sex, original randomized arm, and age <50/50-64/>=65; conditional descriptive benchmark, not full confounding control'),
        ('panel_conditioning','Additional scheme shuffles within the observed four baseline/followup panel-state groups. Panel persistence fixed; assesses assay composition beyond panel-state correspondence'),
        ('exchangeability','Randomization benchmark assumes interchangeability within chosen blocks; persistent demographics or technical factors may contribute to excess'),
        ('interpretation','Excess correspondence does not establish assay repeatability, biology, diagnosis, or treatment utility'),
        ('ocular','Two fixed continuous baseline-adjusted endpoint directions; literature-guided exploration, not prospective preregistration or exact external replication'),
        ('ocular_population','Complete baseline eleven assays plus available baseline/month12 eye outcomes/covariates; does not require month12 antibody measurements'),
        ('censoring','Primary log models use bounds1/160 as in existing Step04; exact-exposure-only sensitivity reported separately'),
        ('multiplicity','Holm6 benchmark comparisons per scheme; Holm2 ocular endpoints separately per exposure-handling and adjustment specification'),
        ('CI','Nominal 95% Wilson observed persistence; nominal HC3/t model intervals; adjusted P values labeled separately'),
        ('literature','Raolji et al. 2026; doi10.1136/bmjophth-2026-002776; endpoints here differ from original progression definitions'),
    ],columns=['item','definition'])
    return dict(methods=methods,permutation_benchmark=bench,permutation_blocks=blocks,
                targeted_ocular_models=mods,source_checks=data.check_sources(),analysis_checks=pd.DataFrame(checks))


if __name__=='__main__':
    run_cli(compute,'STEP10_person_specificity_and_prognosis.xlsx',__doc__)
