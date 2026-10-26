"""Analyses over cell_counts.db (build it first with: python load_data.py).

Run with: python analysis.py
"""

import sqlite3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell_counts.db"
OUTPUT_DIR = ROOT / "outputs"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
ALPHA = 0.05


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
"""


def cell_frequencies(conn):
    return pd.read_sql(FREQUENCY_SQL + " ORDER BY sample_id, population", conn)


# Part 3: frequencies with sample/subject metadata, for cohort comparisons.
COHORT_SQL = f"""
WITH freq AS ({FREQUENCY_SQL})
SELECT freq.*, sa.subject_id, sa.sample_type, sa.time_from_treatment_start,
       su.condition, su.treatment, su.response
FROM freq
JOIN samples  sa ON sa.sample_id  = freq.sample
JOIN subjects su ON su.subject_id = sa.subject_id
WHERE su.condition = ? AND su.treatment = ? AND sa.sample_type = ?
  AND su.response IS NOT NULL
"""


def cohort_frequencies(conn, condition="melanoma", treatment="miraclib", sample_type="PBMC"):
    return pd.read_sql(COHORT_SQL, conn, params=(condition, treatment, sample_type))


def _mwu(df):
    yes = df.loc[df.response == "yes", "percentage"]
    no = df.loc[df.response == "no", "percentage"]
    res = mannwhitneyu(yes, no, alternative="two-sided")
    return yes, no, res


def compare_responders(freqs):
    """Responders vs non-responders, per population.

    Primary test: two-sided Mann-Whitney U on all samples (frequencies are not
    assumed normal), with Benjamini-Hochberg correction across the 5 populations.
    Each subject contributes 3 samples (days 0/7/14), so as a robustness check
    the test is repeated on one value per subject (mean over timepoints), and on
    baseline (day 0) samples alone - the only timepoint usable to predict
    response before treatment.
    """
    by_subject = (freqs.groupby(["subject_id", "response", "population"], as_index=False)
                  .percentage.mean())
    baseline = freqs[freqs.time_from_treatment_start == 0]

    rows = []
    for pop in POPULATIONS:
        yes, no, res = _mwu(freqs[freqs.population == pop])
        rows.append({
            "population": pop,
            "n_responder_samples": len(yes),
            "n_nonresponder_samples": len(no),
            "median_responders": yes.median(),
            "median_nonresponders": no.median(),
            "median_difference": yes.median() - no.median(),
            # Rank-biserial correlation: -1..1, 0 = no separation between groups.
            "effect_size_rbc": 2 * res.statistic / (len(yes) * len(no)) - 1,
            "p_value": res.pvalue,
            "p_value_subject_mean": _mwu(by_subject[by_subject.population == pop])[2].pvalue,
            "p_value_baseline": _mwu(baseline[baseline.population == pop])[2].pvalue,
        })
    out = pd.DataFrame(rows)
    out["q_value_bh"] = multipletests(out.p_value, method="fdr_bh")[1]
    out["significant"] = out.p_value < ALPHA
    out["significant_after_fdr"] = out.q_value_bh < ALPHA
    return out


def summarize(stats):
    sig = stats[stats.significant]
    if sig.empty:
        return f"No population differs significantly between responders and non-responders (p < {ALPHA})."
    lines = []
    for r in sig.itertuples():
        direction = "higher" if r.median_difference > 0 else "lower"
        lines.append(
            f"{r.population}: {direction} in responders "
            f"(median {r.median_responders:.2f}% vs {r.median_nonresponders:.2f}%, "
            f"p = {r.p_value:.4f}, BH q = {r.q_value_bh:.4f}, "
            f"subject-level p = {r.p_value_subject_mean:.4f}, baseline p = {r.p_value_baseline:.4f})."
        )
    if not sig.significant_after_fdr.any():
        lines.append(f"Caveat: no population remains significant after FDR correction (q < {ALPHA}).")
    return "\n".join(lines)


def plot_responders(freqs, stats, path):
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(16, 4.5))
    p_values = stats.set_index("population").p_value
    for ax, pop in zip(axes, POPULATIONS):
        d = freqs[freqs.population == pop]
        ax.boxplot([d.loc[d.response == "yes", "percentage"], d.loc[d.response == "no", "percentage"]],
                   tick_labels=["Responders", "Non-responders"], widths=0.6)
        ax.set_title(f"{pop}\np = {p_values[pop]:.3g}", fontsize=10)
        ax.set_ylabel("Relative frequency (%)")
    fig.suptitle("Melanoma, miraclib, PBMC: responders vs non-responders")
    fig.tight_layout()
    path.parent.mkdir(exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    with connect() as conn:
        freqs = cell_frequencies(conn)
        cohort = cohort_frequencies(conn)

    path = OUTPUT_DIR / "cell_frequencies.csv"
    freqs.to_csv(path, index=False)
    print(f"Part 2: cell population frequencies ({len(freqs)} rows, first 10 shown)\n")
    print(freqs.head(10).to_string(index=False, float_format="%.2f"))
    print(f"\nFull table written to {path.relative_to(ROOT)}")

    stats = compare_responders(cohort)
    stats.to_csv(OUTPUT_DIR / "responder_stats.csv", index=False)
    plot_responders(cohort, stats, OUTPUT_DIR / "responder_boxplots.png")
    n_subjects = cohort.drop_duplicates("subject_id").response.value_counts()
    print(f"\nPart 3: melanoma + miraclib + PBMC, {cohort['sample'].nunique()} samples from "
          f"{n_subjects.sum()} subjects ({n_subjects['yes']} responders, {n_subjects['no']} non-responders)\n")
    cols = ["population", "median_responders", "median_nonresponders", "effect_size_rbc",
            "p_value", "q_value_bh", "p_value_subject_mean", "p_value_baseline"]
    print(stats[cols].to_string(index=False, float_format="%.4f"))
    print("\n" + summarize(stats))
    print("\nWritten outputs/responder_stats.csv and outputs/responder_boxplots.png")


if __name__ == "__main__":
    main()
