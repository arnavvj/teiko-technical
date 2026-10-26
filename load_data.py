"""Create cell_counts.db from cell-count.csv.

Schema (normalized; subject-level metadata is stored once per subject):

    projects     (project_id)
    subjects     (subject_id, project_id, condition, age, sex, treatment, response)
    samples      (sample_id, subject_id, sample_type, time_from_treatment_start)
    cell_counts  (sample_id, population, count)   -- long format, one row per population

Run with: python load_data.py
"""

import csv
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition  TEXT NOT NULL,
    age        INTEGER,
    sex        TEXT CHECK (sex IN ('M', 'F')),
    treatment  TEXT NOT NULL,
    response   TEXT CHECK (response IN ('yes', 'no'))  -- NULL for untreated subjects
);

CREATE TABLE samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type               TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL
);

CREATE TABLE cell_counts (
    sample_id  TEXT NOT NULL REFERENCES samples(sample_id),
    population TEXT NOT NULL,
    count      INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population)
);

CREATE INDEX idx_subjects_filter ON subjects(condition, treatment, response);
CREATE INDEX idx_samples_subject ON samples(subject_id);
CREATE INDEX idx_samples_filter  ON samples(sample_type, time_from_treatment_start);
"""


def load(conn, rows):
    projects, subjects, samples, counts = set(), {}, [], []
    for r in rows:
        projects.add(r["project"])
        subjects[r["subject"]] = (
            r["subject"], r["project"], r["condition"], int(r["age"]),
            r["sex"], r["treatment"], r["response"] or None,
        )
        samples.append((
            r["sample"], r["subject"], r["sample_type"],
            int(r["time_from_treatment_start"]),
        ))
        counts.extend((r["sample"], p, int(r[p])) for p in POPULATIONS)

    conn.executemany("INSERT INTO projects VALUES (?)", [(p,) for p in sorted(projects)])
    conn.executemany("INSERT INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)", subjects.values())
    conn.executemany("INSERT INTO samples VALUES (?, ?, ?, ?)", samples)
    conn.executemany("INSERT INTO cell_counts VALUES (?, ?, ?)", counts)


def main():
    # Rebuild from scratch so re-running the script is idempotent.
    DB_PATH.unlink(missing_ok=True)

    with open(CSV_PATH, newline="") as f:
        rows = list(csv.DictReader(f))

    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    with conn:
        conn.executescript(SCHEMA)
        load(conn, rows)

    for table in ("projects", "subjects", "samples", "cell_counts"):
        n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"{table:12} {n:>6} rows")
    conn.close()
    print(f"Wrote {DB_PATH.name}")


if __name__ == "__main__":
    main()
