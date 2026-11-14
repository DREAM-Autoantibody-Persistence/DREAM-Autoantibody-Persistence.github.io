# DREAM Autoantibody Persistence

**Reproducible analysis code for longitudinal changes in candidate Sjögren autoantibodies in the Dry Eye Assessment and Management (DREAM) trial.**

This repository contains Python scripts for a secondary analysis of publicly available DREAM trial data. It evaluates 12-month changes in nine candidate autoantibody assays: IgG, IgA, and IgM directed against salivary protein 1 (SP1), carbonic anhydrase VI (CA6), and parotid secretory protein (PSP). Anti-SSA/Ro and anti-SSB/La results are used for cohort eligibility and selected reference analyses.

**Scope:** This is an analysis of repeated laboratory results among DREAM participants with moderate-to-severe dry eye disease, **not** a validation study of Sjögren disease diagnosis or clinical monitoring. This repository does not contain individual-level DREAM data.

## Data source

The original DREAM trial (NCT02128763) enrolled 535 participants. Source data and the corresponding study documentation can be accessed from the University of Pennsylvania's official [DREAM Data and Documents page](https://www.med.upenn.edu/cpob/dream-data-export).

Download the following eight workbooks from that official page, retaining these filenames:

| Filename | DREAM dataset | Use in this analysis |
| --- | --- | --- |
| `sjo-results.xlsx` | SJO_RESULTS | Autoantibody concentrations and collection-day codes |
| `gi.xlsx` | GI | Participant characteristics and randomized assignment |
| `mh.xlsx` | MH | Self-reported medical histories |
| `eo.xlsx` | EO | Baseline ocular examinations |
| `oa.xlsx` | OA | Follow-up ocular examinations |
| `osi.xlsx` | OSI | Ocular Surface Disease Index |
| `br.xlsx` | BR | Baseline health review / supplementary checks |
| `wp.xlsx` | WP | Work and activity outcomes / supplementary checks |

The official website also provides case report forms and a laboratory results data dictionary. Please consult those documents for the original variable definitions and coding. Data are made available by the DREAM investigators; obtain them from the source and follow any applicable source-data reuse conditions rather than treating the files as part of this code repository.

Place all eight files in a folder called `data/` in the project directory. The script `dream_common.py` verifies each workbook against `source_manifest.csv` before analysis. The manifest records the SHA-256 checksums of the analysis inputs. On 14 November 2026, all eight files were downloaded again from the official University of Pennsylvania data site and their checksums matched the manifest. If a later official download fails verification, investigate the source-file version rather than bypassing the check.

## Repository structure

```text
DREAM-Autoantibody-Persistence/
├── README.md
├── requirements.txt
├── source_manifest.csv
├── .gitignore
├── dream_common.py
├── RUN_STEP01.py
├── RUN_STEP02.py
├── RUN_STEP03.py
├── RUN_STEP04.py
├── RUN_STEP05.py
├── RUN_STEP06.py
├── RUN_STEP07.py
├── RUN_STEP08.py
├── data/                # Create locally; not included in repository
│   ├── sjo-results.xlsx
│   ├── gi.xlsx
│   └── ... (six additional source workbooks)
└── results/             # Generated automatically; not included in repository
```

`dream_common.py` is required by all eight scripts; it standardizes source validation, laboratory measurements, participant pairing, ocular outcome construction, and Excel exports. The eight `RUN_STEP*.py` files can be executed independently against the same source data; they do **not** read intermediate analysis results.

## Software environment

Python **3.10** is recommended with the versions pinned in `requirements.txt`:

```text
numpy==2.2.6
pandas==2.3.3
openpyxl==3.1.5
scipy==1.15.3
statsmodels==0.15.0
```

On 14 November 2026, all eight public scripts from commit `cd739ea1d3ab58c179d6c602d81e793d3576a162` ran successfully with Python 3.10.0 and the five pinned dependency versions above, using freshly downloaded official workbooks. All 42 output worksheets matched the existing analysis outputs at an absolute numerical tolerance of 1e-10 after treating the provider's `www.med.upenn.edu` and `www.cceb.upenn.edu` source-host names as equivalent metadata. `RUN_STEP01.py` reports the actually installed versions; it checks package availability but does **not** reject all package-version mismatches.

To install dependencies in an isolated environment (commands shown for Windows PowerShell):

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

On macOS/Linux, use `python3.10 -m venv .venv` and `source .venv/bin/activate` before running the installation command.

## Reproduce the analyses

After downloading the eight files into `data/`, run from the repository root:

```powershell
python RUN_STEP01.py
python RUN_STEP02.py
python RUN_STEP03.py
python RUN_STEP04.py
python RUN_STEP05.py
python RUN_STEP06.py
python RUN_STEP07.py
python RUN_STEP08.py
```

The default output directory is `results/`. To use an alternative input or output directory, add `--data-dir PATH` and/or `--output-dir PATH` to an individual command. Each step produces a workbook containing a `Contents` worksheet and its analysis output worksheets.

| Script | Main purpose | Generated workbook |
| --- | --- | --- |
| `RUN_STEP01.py` | Validate data files, required fields, software versions | `STEP01_source_checks.xlsx` |
| `RUN_STEP02.py` | Reconstruct the baseline/month-12 paired cohort | `STEP02_cohort.xlsx` |
| `RUN_STEP03.py` | Panel transitions, persistence, paired net change, bootstrap intervals | `STEP03_main_results.xlsx` |
| `RUN_STEP04.py` | Assay-specific patterns and supporting sensitivity checks | `STEP04_supporting_results.xlsx` |
| `RUN_STEP05.py` | Assemble selected manuscript-facing summary tables | `STEP05_manuscript_tables.xlsx` |
| `RUN_STEP06.py` | Nested persistence, baseline-level models, ocular outcomes, missingness | `STEP06_models_and_sensitivity.xlsx` |
| `RUN_STEP07.py` | Exact two-sided McNemar tests | `STEP07_mcnemar.xlsx` |
| `RUN_STEP08.py` | Exploratory comparisons by original randomized assignment | `STEP08_randomized_arms.xlsx` |

## Analysis definitions

The primary strict cutoff is **>25 EU/mL**. A sensitivity analysis counts borderline results as meeting the cutoff (**≥20 EU/mL**). An individual is panel-positive if at least one of the nine candidate assays meets the relevant cutoff.

Three nested definitions are applied among participants whose candidate panel is positive at baseline:

1. **Panel persistence:** at least one candidate assay meets the cutoff at month 12.
2. **Same-antigen persistence:** at least one antigen that was positive at baseline remains positive at month 12, allowing a different immunoglobulin isotype.
3. **Same-assay persistence:** at least one *identical antigen–isotype combination* is positive at both visits.

Missing assay results (`TNP`, `QNS`) are treated as missing, and the recorded censoring tokens (`<1.0`, `>160`) are represented by bounds for binary cutoff classification. Paired comparisons use matched participant identifiers and baseline/month-12 visit codes. The script comments specify the bootstrap seeds, confidence-interval methods, model covariance estimators, and Holm adjustment families.

## Basic output checks

The intended paired analysis cohort contains **407 participants** with evaluable results for all eleven antibodies at both time points. The following counts provide a rapid verification of key outputs:

| Metric | >25 EU/mL | ≥20 EU/mL |
| --- | ---: | ---: |
| Baseline panel-positive participants | 100 | 147 |
| Panel persistence | 51/100 | 93/147 |
| Same-antigen persistence | 45/100 | 75/147 |
| Same-assay persistence | 38/100 | 69/147 |
| Newly meeting the panel cutoff | 92 | 110 |
| No longer meeting the panel cutoff | 49 | 54 |

These values are descriptive outputs for reproducibility checks. They do not establish diagnostic accuracy, biologic persistence of an antibody response between visits, or clinical utility.

## Reproducibility and limitations

- This is a secondary analysis of the original randomized trial; the present autoantibody analyses should not be interpreted as estimating a causal treatment effect.
- Antibody results were measured at two visits approximately one year apart; the data cannot distinguish analytic variation from biological fluctuation or prove continuous positivity between visits.
- Missing follow-up measurements and exploratory comparisons affect interpretation.
- Participant-level datasets, case report forms, and original data dictionaries belong to the DREAM source provider and are not redistributed here.
- External users should verify that the manifest hashes agree with the official files available to them and record the software versions printed by Step 01.

## Data availability and citation

**Original DREAM data:** [University of Pennsylvania DREAM Data and Documents](https://www.med.upenn.edu/cpob/dream-data-export).

**Original trial:** Dry Eye Assessment and Management Study Research Group. n-3 fatty acid supplementation for the treatment of dry eye disease. *New England Journal of Medicine*. 2018;378:1681–1690. DOI: [10.1056/NEJMoa1709691](https://doi.org/10.1056/NEJMoa1709691).

**Study title:** One-Year Changes in Candidate Sjögren-Related Serum Autoantibodies Among Patients With Moderate-to-Severe Dry Eye: A Longitudinal Secondary Analysis of the DREAM Trial
