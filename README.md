# Loblaw Bio: immune cell population analysis

Analysis of immune cell populations in a clinical trial of the drug candidate **miraclib**, from `cell-count.csv` (10,500 samples, 3,500 subjects, 5 cell populations). The data is loaded into SQLite, analysed in Python, and shown in an interactive Streamlit dashboard.

**Dashboard:** TODO: add hosted link. Run it locally with `make dashboard`, then open <http://localhost:8501>. In GitHub Codespaces, port 8501 is forwarded and opens automatically.

## Data

`cell-count.csv` has one row per biological sample: its metadata plus the cell count for each of 5 immune populations (`b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, `monocyte`). The brief describes the metadata as `sample_id`, `indication`, `treatment`, `time_from_treatment_start`, `response` and `gender`. In the file some of these columns have different names, and there are four extra ones:

| Brief | Column in the file | Values |
|---|---|---|
| sample_id | `sample` | sample00000 … (10,500 unique) |
| indication | `condition` | melanoma, carcinoma, healthy |
| gender | `sex` | M, F |
| treatment | `treatment` | miraclib, phauximab, none |
| time_from_treatment_start | `time_from_treatment_start` | 0, 7, 14 (days) |
| response | `response` | yes, no; empty for the 1,422 untreated healthy samples |
| — | `project` | prj1, prj2, prj3 |
| — | `subject` | 3,500 subjects, 3 samples each |
| — | `age` | integer |
| — | `sample_type` | PBMC, WB |

The code uses the names in the file.

## Running it

The commands below work in GitHub Codespaces or locally with Python 3.9+:

```bash
make setup      # create .venv and install requirements.txt
make pipeline   # build cell_counts.db, run Parts 2-4, write outputs/
make dashboard  # serve the dashboard on port 8501
```

`make clean` deletes the database and outputs. The repo also includes a dev container (`.devcontainer/`, Python 3.12), so a Codespace opens with the right Python version and port 8501 forwarded.

`make pipeline` writes these files:

| File | Contents |
|---|---|
| `cell_counts.db` | SQLite database (Part 1) |
| `outputs/cell_frequencies.csv` | Relative frequency of each population in each sample (Part 2) |
| `outputs/responder_stats.csv` | Statistical tests for responders vs non-responders (Part 3) |
| `outputs/responder_boxplots.png` | Boxplots for responders vs non-responders (Part 3) |
| `outputs/baseline_subset.csv` | Melanoma, miraclib, PBMC samples at day 0 (Part 4) |

It also prints a summary of each part to the terminal.

## Code structure

| File | Role |
|---|---|
| `load_data.py` | Part 1: creates the schema and loads the CSV. Uses only the Python standard library. Each run rebuilds the database from scratch. |
| `analysis.py` | Parts 2–4: SQL queries and statistics as plain functions that take a database connection, plus a `main()` that writes `outputs/`. |
| `dashboard.py` | Streamlit app that calls the same `analysis.py` functions on `cell_counts.db`, so the dashboard and the pipeline cannot disagree. |
| `Makefile` | The `setup`, `pipeline` and `dashboard` targets. |

## Part 1: Database schema

```
projects     (project_id PK)
subjects     (subject_id PK, project_id → projects, condition, age, sex, treatment, response)
samples      (sample_id PK, subject_id → subjects, sample_type, time_from_treatment_start)
cell_counts  (sample_id → samples, population, count)   PK (sample_id, population)
```

**Why this design:**

- **Each fact is stored once.** Condition, age, sex, treatment and response describe a *subject*, not a sample. Every subject has 3 samples (days 0, 7 and 14), so these values live in `subjects` and are never repeated per sample. This also means a subject's samples can't disagree with each other.
- **Cell counts are stored in long format,** with one row per (sample, population) pair. Part 2 (each population's share of a sample's total) is then a single grouped sum in SQL. Adding a new population, for example from a larger flow panel, adds rows; no table needs to change.
- **Constraints keep the data valid.** Foreign keys link the tables, counts can't be negative, `sex` must be M or F, and `response` must be yes, no or NULL. Untreated (healthy) subjects have `response = NULL` rather than a made-up value.
- **There are indexes on the columns used for filtering:** condition, treatment, response, sample type and time point.

**How it scales.** With hundreds of projects and thousands of samples, the shape stays the same. New projects, subjects and samples are new rows, and the indexes keep cohort queries fast. Some natural next steps:

- a `populations` lookup table (marker definitions, lineage),
- a `treatments` table (dose, schedule) once treatments carry more information than a name,
- per-sample QC or assay metadata (panel, instrument, batch) in its own table keyed by `sample_id`,
- PostgreSQL instead of SQLite once several users write concurrently. The schema carries over unchanged.

## Part 2: Cell population frequencies

For each sample, `total_count` is the sum of the 5 population counts, and `percentage = 100 × count / total_count`. The result has 52,500 rows (10,500 samples × 5 populations) with the columns `sample, total_count, population, count, percentage`. It's computed in SQL with a window function (`FREQUENCY_SQL` in `analysis.py`).

## Part 3: Responders vs non-responders

**Cohort:** melanoma patients treated with miraclib, PBMC samples only. That's 1,968 samples from 656 patients: 331 responders and 325 non-responders.

**Method:**

- **Main test:** a two-sided Mann-Whitney U test per population. It doesn't assume the frequencies are normally distributed.
- **Multiple testing:** Benjamini–Hochberg correction across the 5 populations.
- **Effect size:** rank-biserial correlation, where 0 means the two groups don't separate at all.
- **Check for repeated samples:** each patient contributes 3 samples, which are not independent. So the test is repeated with one value per patient (their mean over days 0, 7 and 14).
- **Check at baseline:** the test is repeated on day-0 samples only. That's the only timepoint that could predict response before treatment starts.

| population | median, responders | median, non-responders | p | BH q | p, per patient | p, day 0 |
|---|---|---|---|---|---|---|
| **cd4_t_cell** | **30.22%** | **29.66%** | **0.013** | 0.067 | 0.012 | 0.80 |
| b_cell | 9.43% | 9.79% | 0.056 | 0.139 | 0.35 | 0.55 |
| nk_cell | 14.51% | 14.80% | 0.12 | 0.20 | 0.13 | 0.89 |
| monocyte | 19.61% | 19.94% | 0.16 | 0.20 | 0.26 | 0.21 |
| cd8_t_cell | 24.73% | 24.60% | 0.64 | 0.64 | 0.62 | 0.51 |

**Conclusion:**

- **CD4 T cells are the only population that differs significantly** (p = 0.013). They are slightly higher in responders.
- **The result holds with one value per patient** (p = 0.012), so it isn't caused by counting each patient three times.
- **The effect is small.** The medians differ by about 0.6 percentage points and the two distributions overlap heavily.
- **It does not survive correction for testing 5 populations** (q = 0.067).
- **There's no difference at baseline** (p = 0.80). The CD4 shift appears during treatment, so it is a lead worth following up rather than a marker that predicts response before treatment.

## Part 4: Baseline subset

**Subset:** melanoma PBMC samples at day 0 from patients treated with miraclib: 656 samples from 656 subjects.

| Samples per project | | Subjects by response | | Subjects by sex | |
|---|---|---|---|---|---|
| prj1 | 384 | responders | 331 | male | 344 |
| prj3 | 272 | non-responders | 325 | female | 312 |

prj2 doesn't appear because all of its samples are whole blood (WB), not PBMC. The filters are function arguments (`subset_samples()`), and the dashboard exposes them as dropdowns so other subsets can be explored too.

## Dashboard

The dashboard has three tabs, all read live from `cell_counts.db`:

1. **Cell frequencies:** the Part 2 table, filterable by population and sample, with a per-population summary and a CSV download.
2. **Responders vs non-responders:** interactive boxplots with p-values, the full statistics table, and the conclusion.
3. **Baseline subset:** selectors for condition, treatment, sample type and day. They default to the Part 4 subset and show sample and subject counts by project, response and sex.
