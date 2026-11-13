"""Shared data handling and Excel export for the DREAM antibody analyses.

All paths are relative to this code directory unless explicitly supplied.
Every analysis starts from the official workbooks; no previous results are inputs.
"""
import argparse
import hashlib
import math
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl.styles import Alignment, Font
from statsmodels.stats.proportion import proportion_confint

PROJECT_DIR = Path(__file__).resolve().parent
KEY = "alpha_code"
ASSAYS = ["ss_a_ro", "ss_b_la", "sp_1_igg", "sp_1_iga", "sp_1_igm",
          "ca_vi_igg", "ca_vi_iga", "ca_vi_igm", "psp_igg", "psp_iga", "psp_igm"]
NOVEL = ASSAYS[2:]
ANTIGENS = {"SP1": NOVEL[:3], "CA6": NOVEL[3:6], "PSP": NOVEL[6:]}
MISSING = {"TNP", "QNS"}
CENSORED = {"<1.0", ">160"}
OUTCOMES = ["conjunctival_staining", "corneal_staining", "osdi", "schirmer"]
PACKAGES = ["numpy", "pandas", "openpyxl", "scipy", "statsmodels"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fraction(n, denominator):
    """Wilson intervals; all percentage fields are on the 0--100 scale."""
    if denominator == 0:
        return dict(numerator=n, denominator=0, percent=np.nan,
                    ci_low_percent=np.nan, ci_high_percent=np.nan)
    low, high = proportion_confint(n, denominator, method="wilson")
    return dict(numerator=n, denominator=denominator, percent=100*n/denominator,
                ci_low_percent=float(100*low), ci_high_percent=float(100*high))


def transition(rule, marker, first, second):
    require(first.index.equals(second.index), "Unaligned paired observations")
    return dict(rule=rule, marker=marker, n=len(first),
                repeat_pos=int((first & second).sum()),
                pos_to_nonpos=int((first & ~second).sum()),
                nonpos_to_pos=int((~first & second).sum()),
                repeat_nonpos=int((~first & ~second).sum()))


class DreamData:
    """Read and validate the frozen public data without modifying the sources."""
    def __init__(self, data_dir=None):
        self.data_dir = Path(data_dir) if data_dir else PROJECT_DIR / "data"
        self.manifest = pd.read_csv(PROJECT_DIR / "source_manifest.csv", dtype=str)
        self._cache = {}
        self._checks = None
        self._lab = None

    def check_sources(self):
        if self._checks is not None:
            return self._checks.copy()
        rows = []
        for item in self.manifest.to_dict("records"):
            path = self.data_dir / item["File"]
            require(path.is_file(), f"Missing data/{item['File']}. See README.md or run DOWNLOAD_DATA.py.")
            actual = digest(path)
            require(actual == item["Hash"].lower(),
                    f"Checksum mismatch: {item['File']}. Use the documented frozen source version.")
            rows.append(dict(file=item["File"], sha256=actual, checksum_match=True,
                             source_url=item["SourceUrl"]))
        self._checks = pd.DataFrame(rows)
        return self._checks.copy()

    def read_sheet(self, name, columns=None):
        self.check_sources()
        if name not in self._cache:
            frame = pd.read_excel(self.data_dir / name, sheet_name=0,
                                  dtype={KEY: str, "month": str})
            frame = frame.dropna(subset=[KEY])
            self._cache[name] = frame
        frame = self._cache[name]
        if columns is not None:
            missing = set(columns) - set(frame.columns)
            require(not missing, f"Missing columns in {name}: {sorted(missing)}")
            frame = frame[columns]
        return frame.copy()

    def read(self, name, columns, visit=None):
        fields = [KEY, *columns]
        if visit is not None:
            fields.append("month")
        frame = self.read_sheet(name, fields)
        if visit is not None:
            frame = frame.loc[frame.month.eq(visit)].drop(columns="month")
        require(frame[KEY].notna().all() and frame[KEY].is_unique,
                f"Invalid participant/visit keys in {name}")
        return frame.set_index(KEY)

    def numeric_lab(self):
        if self._lab is None:
            frame = self.read_sheet("sjo-results.xlsx", [KEY, "month", *ASSAYS])
            require(not frame.duplicated([KEY, "month"]).any(), "Duplicate antibody participant/visit")
            require(set(frame.month) == {"00", "12", "24"}, "Unexpected antibody visit")
            for name in ASSAYS:
                token = frame[name].astype(str).str.strip().str.upper()
                number = pd.to_numeric(token, errors="coerce")
                require((number.between(0, 160) | token.isin(MISSING | CENSORED)).all(),
                        f"Unrecognized assay encoding: {name}")
                # These bound proxies preserve status at 20 and 25 EU/mL.
                # Continuous descriptive summaries use interval bounds instead.
                frame[name] = number.mask(token.eq("<1.0"), 1.0).mask(token.eq(">160"), 160.0)
            count = frame[ASSAYS].isna().sum(axis=1)
            require(((count == 0) | (count == len(ASSAYS))).all(),
                    "Partially missing antibody visit requires an explicit eligibility decision")
            ids = self.read_sheet("gi.xlsx", [KEY])[KEY]
            require(ids.is_unique and frame[KEY].isin(ids).all(), "Unlinked antibody participant")
            self._lab = frame
        return self._lab.copy()

    def baseline(self):
        frame = self.numeric_lab()
        return frame.loc[frame.month.eq("00")].dropna(subset=ASSAYS).set_index(KEY)[ASSAYS]

    def pairs(self):
        frame = self.numeric_lab().dropna(subset=ASSAYS)
        base = frame.loc[frame.month.eq("00")].set_index(KEY)[ASSAYS]
        follow = frame.loc[frame.month.eq("12")].set_index(KEY)[ASSAYS]
        ids = sorted(base.index.intersection(follow.index))
        require(len(ids) > 0, "No complete baseline/month-12 pairs")
        return base.loc[ids], follow.loc[ids]

    def collection_intervals(self):
        ids = self.pairs()[0].index
        raw = self.read_sheet("sjo-results.xlsx", [KEY, "month", "collection_date", "srbldtc"])
        require(raw.collection_date.eq(raw.srbldtc).all(), "Collection-day fields disagree")
        dates = raw.pivot(index=KEY, columns="month", values="collection_date").loc[ids]
        days = dates["12"] - dates["00"]
        require(days.notna().all() and days.gt(0).all(), "Invalid paired collection interval")
        return days

    def eyes(self, visit):
        left = ["fluotople", "fluocentle", "fluotemple", "fluobottle", "fluonasrtle"]
        right = ["fluotoprt", "fluocentrt", "fluotemprt", "fluobottrt", "fluonasrtrt"]
        cleft, cright = ["lisstemple", "lissnasle"], ["lisstemprt", "lissnasrt"]
        schirmer = ["mmwetle", "mmwetrt"]
        fields = [*left, *right, *cleft, *cright, *schirmer]
        frame = self.read("eo.xlsx" if visit == "00" else "oa.xlsx", fields, visit)
        require(frame[fields[:-2]].stack().isin([0, 1, 2, 3]).all(), "Invalid staining score")
        require(frame[schirmer].stack().between(0, 36).all(), "Invalid Schirmer value")
        return pd.DataFrame({
            "conjunctival_staining": (frame[cleft].sum(axis=1, min_count=2) +
                                      frame[cright].sum(axis=1, min_count=2))/2,
            "corneal_staining": (frame[left].sum(axis=1, min_count=5) +
                                  frame[right].sum(axis=1, min_count=5))/2,
            "schirmer": frame[schirmer].mean(axis=1, skipna=False),
            "osdi": self.read("osi.xlsx", ["osdi_score"], visit).osdi_score,
        }).add_suffix("_" + visit)


def export_excel(tables, path):
    """Write aggregate tables only. No participant identifiers are exported."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    contents, names = [], set()
    for key, frame in tables.items():
        require(isinstance(frame, pd.DataFrame), f"Not a result table: {key}")
        require(KEY not in frame.columns, "Participant identifiers cannot be exported")
        short = key[:31]
        require(short not in names, f"Excel worksheet name collision: {short}")
        names.add(short)
        contents.append(dict(table=key, worksheet=short, rows=len(frame)))
    temporary = path.with_name(path.stem + ".pending.xlsx")
    try:
        with pd.ExcelWriter(temporary, engine="openpyxl") as writer:
            pd.DataFrame(contents).to_excel(writer, sheet_name="Contents", index=False)
            for key, frame in tables.items():
                frame.to_excel(writer, sheet_name=key[:31], index=False, inf_rep="Infinity")
            for sheet in writer.book.worksheets:
                sheet.freeze_panes = "A2"
                sheet.auto_filter.ref = sheet.dimensions
                sheet.row_dimensions[1].height = 42
                for cell in sheet[1]:
                    cell.font = Font(name="Calibri", size=11, bold=True)
                    cell.alignment = Alignment(vertical="center", wrap_text=True)
                for col in sheet.columns:
                    values = list(col)
                    width = max(14, min(44, max(len(str(x.value or "")) for x in values[:50]) + 2))
                    sheet.column_dimensions[values[0].column_letter].width = width
                    for cell in values[1:]:
                        cell.font = Font(name="Calibri", size=11)
                        cell.alignment = Alignment(vertical="center", wrap_text=isinstance(cell.value, str))
                        if isinstance(cell.value, float):
                            cell.number_format = "0.000000E+00" if str(values[0].value).startswith("p_") else "0.000000"
                for row in sheet.iter_rows(min_row=2):
                    lines=max([1]+[math.ceil(len(cell.value)/max(10,sheet.column_dimensions[cell.column_letter].width-2))
                                   for cell in row if isinstance(cell.value,str)])
                    sheet.row_dimensions[row[0].row].height=max(18,lines*14+4)
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()
    print("Saved " + path.name)


def parser(description):
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--data-dir", type=Path, default=PROJECT_DIR / "data",
                   help="Directory containing the eight official source workbooks")
    p.add_argument("--output-dir", type=Path, default=PROJECT_DIR / "results",
                   help="Directory for aggregate Excel workbooks")
    return p


def run_cli(compute, filename, description):
    args = parser(description).parse_args()
    data = DreamData(args.data_dir)
    tables = compute(data)
    export_excel(tables, args.output_dir / filename)
