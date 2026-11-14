# DREAM autoantibody analysis

Python scripts for the baseline and month-12 candidate autoantibody analyses in the DREAM trial.

## Data

Download the following files from the [official DREAM data page](https://www.med.upenn.edu/cpob/dream-data-export) and place them in `data/` next to the scripts, keeping their original filenames:

```text
sjo-results.xlsx
gi.xlsx
mh.xlsx
eo.xlsx
oa.xlsx
osi.xlsx
br.xlsx
wp.xlsx
```

The scripts verify these files against `source_manifest.csv`. Participant-level data are obtained from the original provider and are not included in this repository.

## Run

Use Python 3.10 and the versions in `requirements.txt`. In Windows PowerShell, open the repository folder and run:

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe RUN_STEP01.py
.\.venv\Scripts\python.exe RUN_STEP02.py
.\.venv\Scripts\python.exe RUN_STEP03.py
.\.venv\Scripts\python.exe RUN_STEP04.py
.\.venv\Scripts\python.exe RUN_STEP05.py
.\.venv\Scripts\python.exe RUN_STEP06.py
.\.venv\Scripts\python.exe RUN_STEP07.py
.\.venv\Scripts\python.exe RUN_STEP08.py
```

On macOS/Linux, create the environment with `python3.10 -m venv .venv` and use `.venv/bin/python` in place of `.\.venv\Scripts\python.exe`.

Each step can also run independently. Keep `dream_common.py`, `requirements.txt` and `source_manifest.csv` in the same folder as the scripts. Paths default to that folder, regardless of the working directory.

## Outputs

Steps 01–08 write eight aggregate Excel workbooks to `results/`. Each workbook has a `Contents` sheet listing its tables. Close an output workbook before rerunning its step; the file will be overwritten.

Alternative folders can be supplied to any step, for example:

```powershell
.\.venv\Scripts\python.exe RUN_STEP01.py --data-dir "path/to/data" --output-dir "path/to/results"
```

## Analysis definitions

Candidate-antibody positivity is >25 EU/mL in the main analysis; the ≥20 EU/mL sensitivity threshold includes borderline results. Panel, same-antigen and same-assay persistence use the corresponding baseline panel-positive participants as the denominator. Same-assay persistence requires an identical antigen–isotype combination at both visits. Missing antibody results are not imputed.
