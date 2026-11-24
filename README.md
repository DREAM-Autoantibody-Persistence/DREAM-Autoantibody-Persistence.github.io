# DREAM autoantibody persistence analysis

Reproducible Python analysis for **One-year changes in candidate Sjögren-related autoantibodies in the DREAM dry eye cohort**.

## Run

Tested with Python 3.10 and the versions in `requirements.txt`. Run the following commands from this folder in an activated Python environment:

```text
python -m pip install -r requirements.txt
python DOWNLOAD_DATA.py
python RUN_ALL.py
```

The downloader obtains eight workbooks from the [official DREAM data page](https://www.med.upenn.edu/cpob/dream-data-export) and checks their SHA256 values against `source_manifest.csv`. Existing matching files are reused. Source data are not included in this repository or the code package.

If the workbooks are already stored elsewhere:

```text
python RUN_ALL.py --data-dir /path/to/workbooks --output-dir /path/to/results
```

Each script can also be run separately, for example `python RUN_STEP09.py`. Outputs are ten Excel workbooks in `results/`; they contain aggregate results, not participant records. No previous analysis output is required as input.

| Step | Contents |
|---|---|
| 01 | Source files, software versions and data checks |
| 02 | Paired cohort and participant flow |
| 03 | Panel transitions and nested persistence estimates |
| 04 | Assay-level and other supporting analyses |
| 05 | Main and supplementary table summaries |
| 06 | Baseline concentration and ocular-outcome models |
| 07 | Paired McNemar tests |
| 08 | Exploratory randomised-arm contrasts |
| 09 | Quantitative rank correlations; additional month-24 feasibility summaries |
| 10 | Whole-panel reassignment benchmarks; additional targeted ocular checks |

## Additional analyses reported in the manuscript

Table S6 uses `quantitative_tracking` in Step 09. Exact-value correlations exclude results reported only as <1.0 or >160 EU/mL at either visit. Their intervals use 2,000 participant-pair bootstrap samples. Bound-proxy correlations are sensitivity estimates.

Table S7 uses `permutation_benchmark` in Step 10: 20,000 random reassignments of whole month-12 nine-assay vectors under each scheme and threshold. The observed-panel-transition scheme preserves all four panel-status groups. Its reference ranges describe reassignment distributions, not confidence intervals for observed persistence. These analyses assess longitudinal correspondence, not laboratory repeatability or external validation.

The supplementary month-24 feasibility summaries and targeted ocular checks remain available in the workbooks but are not reported as manuscript validation. Seeds and statistical definitions are recorded in each workbook's `methods` sheet and in the scripts.
