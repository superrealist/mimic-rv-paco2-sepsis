"""Run the authorized-data primary extraction pipeline; does not fit models."""
import argparse
from . import extract_core, extract_raw_gas, culture_time, build_candidates
from . import extract_baseline, build_baseline, finalize_cohort


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["all", "core", "baseline", "course"], default="all")
    args = parser.parse_args()
    if args.stage in ("all", "core"):
        extract_core.main()
        extract_raw_gas.main()
        culture_time.extract()
        culture_time.summarize()
        build_candidates.main()
    if args.stage in ("all", "baseline"):
        extract_baseline.main()
        extract_baseline.audit_sofa_availability()
        build_baseline.main()
    if args.stage in ("all", "course"):
        finalize_cohort.prepare_input()
        finalize_cohort.extract_pressor_intervals()
        finalize_cohort.build_final()


if __name__ == "__main__":
    main()
