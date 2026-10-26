"""Analyses over cell_counts.db (build it first with: python load_data.py).

Run with: python analysis.py
"""

import csv
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"
OUTPUT_DIR = ROOT / "outputs"


def connect():
    if not DB_PATH.exists():
        raise SystemExit(f"{DB_PATH.name} not found - run `python load_data.py` first.")
    return sqlite3.connect(DB_PATH)


# Part 2: relative frequency of each population in each sample.
FREQUENCY_SQL = """
SELECT
    sample_id                                   AS sample,
    SUM(count) OVER (PARTITION BY sample_id)    AS total_count,
    population,
    count,
    100.0 * count / SUM(count) OVER (PARTITION BY sample_id) AS percentage
FROM cell_counts
ORDER BY sample_id, population
"""
FREQUENCY_COLUMNS = ["sample", "total_count", "population", "count", "percentage"]


def cell_frequencies(conn):
    return conn.execute(FREQUENCY_SQL).fetchall()


def write_csv(path, columns, rows):
    path.parent.mkdir(exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(rows)


def print_table(columns, rows):
    fmt = "{:<12} {:>11}  {:<11} {:>7} {:>10}"
    print(fmt.format(*columns))
    for r in rows:
        print(fmt.format(r[0], r[1], r[2], r[3], f"{r[4]:.2f}"))


def main():
    with connect() as conn:
        freqs = cell_frequencies(conn)

    path = OUTPUT_DIR / "cell_frequencies.csv"
    write_csv(path, FREQUENCY_COLUMNS, freqs)
    print(f"Part 2: cell population frequencies ({len(freqs)} rows, first 10 shown)\n")
    print_table(FREQUENCY_COLUMNS, freqs[:10])
    print(f"\nFull table written to {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
