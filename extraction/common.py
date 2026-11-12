"""Portable, read-only PostgreSQL export helpers.

Connection settings are supplied through libpq environment variables or a local
passfile. Nothing in this module reads the original workstation configuration.
"""
from pathlib import Path
import os
import subprocess

REPO = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("MIMIC_DATA_DIR", REPO / "data")).expanduser().resolve()
OUTPUT_DIR = Path(os.environ.get("MIMIC_OUTPUT_DIR", REPO / "outputs")).expanduser().resolve()
ROOT = DATA_DIR / "stages" / "primary"
for _name in ("data", "tables", "logs", "sql"):
    (ROOT / _name).mkdir(parents=True, exist_ok=True)


def copy(query, target):
    """Build a client-side CSV export inside the private primary stage folder."""
    path = (ROOT / target).resolve()
    if ROOT.resolve() not in path.parents:
        raise ValueError("Export target must be inside the primary stage folder")
    path.parent.mkdir(parents=True, exist_ok=True)
    quoted = path.as_posix().replace("'", "''")
    compact = " ".join(query.strip().rstrip(";").splitlines())
    return "\\copy (" + compact + ") TO '" + quoted + "' CSV HEADER\n"


def run_sql(name, body):
    """Run generated SELECT/COPY statements under a read-only connection."""
    path = ROOT / "sql" / (name + ".sql")
    path.write_text(
        "\\set ON_ERROR_STOP on\nSET default_transaction_read_only=on;\n"
        "SET statement_timeout='30min';\n" + body,
        encoding="utf-8",
    )
    command = [os.environ.get("MIMIC_PSQL", "psql"), "-X", "-w", "-v",
               "ON_ERROR_STOP=1", "--file", str(path)]
    result = subprocess.run(command, cwd=REPO, capture_output=True)
    # The local log may identify the user's own connection or rows. It remains
    # beneath MIMIC_DATA_DIR and must not be committed or shared.
    out = result.stdout.decode("utf-8", errors="replace")
    err = result.stderr.decode("utf-8", errors="replace")
    (ROOT / "logs" / (name + ".log")).write_text(out + err, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"Read-only export {name!r} failed; inspect the private stage log")
    print(name + ": " + "; ".join(x for x in out.splitlines() if x.startswith("COPY ")), flush=True)


def copy_query(query, target, name="export"):
    """Export a query to an explicit path below MIMIC_DATA_DIR (outcome API)."""
    path = Path(target).expanduser().resolve()
    if DATA_DIR not in path.parents:
        raise ValueError("Export target must be inside MIMIC_DATA_DIR")
    path.parent.mkdir(parents=True, exist_ok=True)
    compact = " ".join(query.strip().rstrip(";").splitlines())
    quoted = path.as_posix().replace("'", "''")
    run_sql(name, "\\copy (" + compact + ") TO '" + quoted + "' CSV HEADER\n")
