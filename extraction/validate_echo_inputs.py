"""Validate prepared echo-table CSVs without inventing study IDs or times."""
import argparse
from pathlib import Path
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("events", type=Path)
    p.add_argument("measurements", type=Path)
    args = p.parse_args()
    e = pd.read_csv(args.events)
    m = pd.read_csv(args.measurements)
    keys = ["measurement_id", "subject_id", "test_type", "measurement_datetime"]
    for d, required in ((e, keys), (m, keys + ["measurement", "result"])):
        missing = sorted(set(required) - set(d.columns))
        if missing:
            raise ValueError("Missing required columns: " + ", ".join(missing))
        if d[keys].isna().any().any():
            raise ValueError("Study linkage keys must be nonmissing")
        for col in ("measurement_id", "subject_id"):
            n = pd.to_numeric(d[col], errors="raise")
            if not n.mod(1).eq(0).all():
                raise ValueError(col + " must contain source-assigned integer identifiers")
        # Date-only strings are rejected, not silently converted to midnight.
        if not d.measurement_datetime.astype(str).str.match(r"^\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}").all():
            raise ValueError("Timestamp must retain source-recorded date and time")
        d["measurement_datetime"] = pd.to_datetime(d.measurement_datetime, errors="raise")
    if not e.measurement_id.is_unique:
        raise ValueError("measurement_events must have one row per measurement_id")
    check = m[keys].merge(e[keys], on=keys, how="left", indicator=True, validate="many_to_one")
    if not check._merge.eq("both").all():
        raise ValueError("Measurement linkage differs from the prepared event table")
    print(f"Validated {len(e)} prepared studies and {len(m)} measurement rows")


if __name__ == "__main__":
    main()
