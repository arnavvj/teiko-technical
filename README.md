# Loblaw Bio: immune cell population analysis

Analysis of immune cell populations in a clinical trial of the drug candidate **miraclib**. The data is loaded from `cell-count.csv` (10,500 samples, 3,500 subjects, 5 cell populations) into SQLite, analysed in Python, and shown in a Streamlit dashboard.

**Dashboard:** <https://teiko-technical-arnav.streamlit.app/>. To run it locally, use `make dashboard` and open <http://localhost:8501>. In Codespaces, port 8501 is forwarded automatically.

## Data

Each row is one sample: metadata plus counts for `b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell` and `monocyte`. The code uses the file's column names:

| Brief | File column | Values |
|---|---|---|
| sample_id | `sample` | 10,500 unique |
| indication | `condition` | melanoma, carcinoma, healthy |
| gender | `sex` | M, F |
| response | `response` | yes, no; empty for the 1,422 untreated healthy samples |
| — | `project`, `subject`, `age`, `sample_type` | prj1–3; 3,500 subjects × 3 samples; integer; PBMC/WB |

`treatment` (miraclib, phauximab, none) and `time_from_treatment_start` (0, 7, 14 days) are named as in the brief.

## Running it

These commands work in Codespaces (`.devcontainer/` uses Python 3.12) or locally with Python 3.9+:

```bash
make setup      # create .venv, install requirements.txt
make pipeline   # build cell_counts.db, run Parts 2-4, write outputs/
make dashboard  # serve the dashboard on port 8501
```

`make clean` deletes generated files. The pipeline writes `cell_counts.db` and these outputs: `cell_frequencies.csv` (Part 2), `responder_stats.csv` and `responder_boxplots.png` (Part 3), and `baseline_subset.csv` (Part 4).

| Code | Role |
|---|---|
| `load_data.py` | Part 1: schema and CSV load (standard library only; rebuilds on each run) |
| `analysis.py` | Parts 2–4: SQL and statistics functions; writes `outputs/` |
| `dashboard.py` | Streamlit app calling the same `analysis.py` functions, so the two always agree |

## Part 1: Database schema

```
projects     (project_id PK)
subjects     (subject_id PK, project_id → projects, condition, age, sex, treatment, response)
samples      (sample_id PK, subject_id → subjects, sample_type, time_from_treatment_start)
cell_counts  (sample_id → samples, population, count)   PK (sample_id, population)
```

- **Each fact is stored once.** Condition, age, sex, treatment and response describe a subject, not a sample, so they live in `subjects` instead of being repeated on each of the subject's 3 samples.
- **Counts are in long format,** with one row per sample and population. Frequencies are then a single grouped sum, and a new population adds rows without a schema change.
- **Constraints keep the data valid:** foreign keys, `count >= 0`, `sex` in (M, F), and `response` in (yes, no, NULL), with NULL for untreated subjects.
- **Indexes** cover condition, treatment, response, sample type and time point.

**Scaling:** more projects and samples are just more rows. Natural next steps are lookup tables for populations and treatments, a per-sample QC/assay table, and PostgreSQL once several users write concurrently.

## Part 2: Cell population frequencies

`total_count` is the sum of a sample's 5 population counts, and `percentage = 100 × count / total_count`. The result has 52,500 rows (10,500 samples × 5 populations) and is computed with a SQL window function (`FREQUENCY_SQL`).

## Part 3: Responders vs non-responders

**Cohort:** melanoma patients on miraclib, PBMC samples only. That's 1,968 samples from 656 patients (331 responders, 325 non-responders).

**Method:**
- **Test:** two-sided Mann-Whitney U per population, which doesn't assume normality.
- **Multiple testing:** Benjamini–Hochberg correction across the 5 populations.
- **Effect size:** rank-biserial correlation.
- **Repeated samples:** each patient contributes 3 non-independent samples, so the test is repeated with one value per patient (mean of days 0, 7 and 14).
- **Baseline:** the test is repeated on day-0 samples only, the only timepoint that could predict response before treatment.

| population | median, responders | median, non-responders | p | BH q | p, per patient | p, day 0 |
|---|---|---|---|---|---|---|
| **cd4_t_cell** | **30.22%** | **29.66%** | **0.013** | 0.067 | 0.012 | 0.80 |
| b_cell | 9.43% | 9.79% | 0.056 | 0.139 | 0.35 | 0.55 |
| nk_cell | 14.51% | 14.80% | 0.12 | 0.20 | 0.13 | 0.89 |
| monocyte | 19.61% | 19.94% | 0.16 | 0.20 | 0.26 | 0.21 |
| cd8_t_cell | 24.73% | 24.60% | 0.64 | 0.64 | 0.62 | 0.51 |

**Conclusion:**
- **CD4 T cells are the only population that differs significantly.** They are slightly higher in responders (p = 0.013), and the result holds with one value per patient (p = 0.012).
- **The effect is small:** about 0.6 percentage points, with heavily overlapping distributions.
- **It doesn't survive correction for testing 5 populations** (q = 0.067).
- **There's no difference at baseline** (p = 0.80), so the shift appears during treatment and does not predict response beforehand. It's a lead to follow up, not a validated marker.

## Part 4: Baseline subset

**Subset:** melanoma PBMC samples at day 0 from patients on miraclib. That's 656 samples from 656 subjects.

| Samples per project | | Subjects by response | | Subjects by sex | |
|---|---|---|---|---|---|
| prj1 | 384 | responders | 331 | male | 344 |
| prj3 | 272 | non-responders | 325 | female | 312 |

prj2 doesn't appear because all of its samples are whole blood (WB).

## Dashboard

The dashboard has three tabs, all read live from `cell_counts.db`:
1. **Cell frequencies:** the Part 2 table, with filters and a CSV download.
2. **Responders vs non-responders:** the Part 3 boxplots, statistics and conclusion.
3. **Baseline subset:** the Part 4 counts, with dropdowns to explore other subsets.
