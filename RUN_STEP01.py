"""Step 01: validate official source files and report the software environment."""
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def check_dependencies():
    """Explain how to prepare this exact interpreter before scientific imports."""
    requirements = Path(__file__).resolve().parent / "requirements.txt"
    missing = []
    for line in requirements.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        package = line.split("==", 1)[0]
        try:
            version(package)
        except PackageNotFoundError:
            missing.append(package)
    if missing:
        command = f'& "{sys.executable}" -m pip install -r "{requirements}"'
        raise SystemExit(
            "Missing Python dependencies: " + ", ".join(missing)
            + "\nInterpreter: " + sys.executable
            + "\nRun this command in PowerShell, then rerun Step 01:\n" + command
        )


if __name__ == "__main__":
    check_dependencies()

import pandas as pd
from dream_common import ASSAYS, PACKAGES, run_cli

SCHEMAS = {
    "sjo-results.xlsx": ["alpha_code", "month", "collection_date", *ASSAYS, "ana_titer"],
    "gi.xlsx": ["alpha_code", "age", "sex", "group", "treatment"],
    "eo.xlsx": ["alpha_code", "month", "mmwetrt", "mmwetle"],
    "oa.xlsx": ["alpha_code", "month", "mmwetrt", "mmwetle"],
    "osi.xlsx": ["alpha_code", "month", "osdi_score"],
    "mh.xlsx": ["alpha_code", "mhsjog", "mhrheu", "mhthyr", "mhdiab", "cmcort", "cmcortnow"],
    "br.xlsx": ["alpha_code", "month", "sucigs", "susmcur"],
    "wp.xlsx": ["alpha_code", "month", "opwork", "opwprod", "opacth"],
}

def compute(data):
    sources = data.check_sources()
    schemas = []
    for name, fields in SCHEMAS.items():
        frame = data.read_sheet(name, fields)
        schemas.append(dict(file=name, rows=len(frame), required_fields_present=True))
    software = [dict(component="Python", version=platform.python_version())]
    software.extend(dict(component=n, version=version(n)) for n in PACKAGES)
    data.pairs()  # Validate coding, linkage and completeness before calculation.
    return {"source_checks": sources, "data_schemas": pd.DataFrame(schemas),
            "software_versions": pd.DataFrame(software)}

if __name__ == "__main__":
    run_cli(compute, "STEP01_source_checks.xlsx", __doc__)
