"""Reproduce steps 01--10 from the official source workbooks."""
from importlib import import_module
from dream_common import DreamData, export_excel, parser

OUTPUTS = ["STEP01_source_checks.xlsx", "STEP02_cohort.xlsx", "STEP03_main_results.xlsx",
           "STEP04_supporting_results.xlsx", "STEP05_manuscript_tables.xlsx",
           "STEP06_models_and_sensitivity.xlsx", "STEP07_mcnemar.xlsx",
           "STEP08_randomized_arms.xlsx",
           "STEP09_quantitative_and_month24.xlsx",
           "STEP10_person_specificity_and_prognosis.xlsx"]

def main():
    args = parser(__doc__).parse_args()
    data = DreamData(args.data_dir)
    for number, filename in enumerate(OUTPUTS, 1):
        module = import_module(f"RUN_STEP{number:02d}")
        export_excel(module.compute(data), args.output_dir / filename)

if __name__ == "__main__":
    main()
