"""Validate the frozen protocol configuration schema (no secrets, required keys).

Usage: python3 scripts/validate_config.py [path/to/config.yaml]
Exit code 0 = valid; 1 = missing/invalid keys.
"""
import sys
import yaml

REQUIRED_TOP = ["protocol", "study", "clinvar", "gnomad", "predictors",
                "thresholds", "evidence", "resources", "output_dirs"]

REQUIRED_STUDY = ["gene", "genome_build", "transcript"]
REQUIRED_CLINVAR = ["file", "include_significance", "include_consequence"]
REQUIRED_PREDICTOR_KEYS = {"name", "source", "range", "higher_is_damaging"}


def main() -> int:
    path = sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"
    with open(path) as fh:
        cfg = yaml.safe_load(fh)

    problems = []
    for k in REQUIRED_TOP:
        if k not in cfg:
            problems.append(f"missing top-level key: {k}")
    if "study" in cfg:
        for k in REQUIRED_STUDY:
            if k not in cfg["study"]:
                problems.append(f"missing study.{k}")
    if "clinvar" in cfg:
        for k in REQUIRED_CLINVAR:
            if k not in cfg["clinvar"]:
                problems.append(f"missing clinvar.{k}")
    if "predictors" in cfg:
        for p in cfg["predictors"]:
            missing = REQUIRED_PREDICTOR_KEYS - set(p.keys())
            if missing:
                problems.append(f"predictor '{p.get('name')}' missing {sorted(missing)}")
    if "thresholds" in cfg:
        for k in ("ba1_af", "brca1_bs1_faf", "revel_bp4", "revel_pp3"):
            if k not in cfg["thresholds"]:
                problems.append(f"missing thresholds.{k}")

    if problems:
        print("CONFIG INVALID:")
        for p in problems:
            print("  -", p)
        return 1
    print(f"CONFIG VALID: {path}")
    print("  gene =", cfg["study"]["gene"], "| transcript =", cfg["study"]["transcript"])
    print("  predictors =", [p["name"] for p in cfg["predictors"]])
    return 0


if __name__ == "__main__":
    sys.exit(main())
