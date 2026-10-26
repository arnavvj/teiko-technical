"""Interactive dashboard for the Part 2-4 results, read from cell_counts.db.

Run with: make dashboard   (or: streamlit run dashboard.py)
"""

import sqlite3

import plotly.express as px
import streamlit as st

import analysis
import load_data

RESPONSE_COLORS = {"yes": "#2a78d6", "no": "#eb6834"}
RESPONSE_LABELS = {"yes": "Responders", "no": "Non-responders"}

st.set_page_config(page_title="Loblaw Bio: immune cell populations", layout="wide")


@st.cache_resource
def ensure_database():
    # The pipeline builds the database; build it here too so a fresh deploy still works.
    if not analysis.DB_PATH.exists():
        load_data.main()


def query(fn, *args):
    ensure_database()
    with sqlite3.connect(analysis.DB_PATH) as conn:
        return fn(conn, *args)


@st.cache_data
def frequencies():
    return query(analysis.cell_frequencies)


@st.cache_data
def responder_cohort():
    cohort = query(analysis.cohort_frequencies)
    return cohort, analysis.compare_responders(cohort)


@st.cache_data
def subset(condition, treatment, sample_type, time):
    return query(analysis.subset_samples, condition, treatment, sample_type, time)


@st.cache_data
def distinct(column, table):
    return query(lambda conn: [r[0] for r in conn.execute(
        f"SELECT DISTINCT {column} FROM {table} ORDER BY {column}")])


def default_index(options, value):
    return options.index(value) if value in options else 0


st.title("Immune cell populations in the miraclib trial")
tab_freq, tab_resp, tab_subset = st.tabs(
    ["Cell frequencies", "Responders vs non-responders", "Baseline subset"])

# --- Part 2 ---------------------------------------------------------------
with tab_freq:
    freqs = frequencies()
    st.subheader("Relative frequency of each population in each sample")
    st.caption(f"{freqs['sample'].nunique():,} samples x 5 populations = {len(freqs):,} rows. "
               "Percentage = population count / total count of the sample.")

    pops = st.multiselect("Populations", analysis.POPULATIONS, default=analysis.POPULATIONS)
    samples = st.multiselect("Samples (leave empty for all)", sorted(freqs["sample"].unique()))
    shown = freqs[freqs.population.isin(pops)]
    if samples:
        shown = shown[shown["sample"].isin(samples)]

    overview = (shown.groupby("population").percentage
                .agg(mean="mean", median="median", min="min", max="max")
                .reindex([p for p in analysis.POPULATIONS if p in pops]))
    st.dataframe(overview.style.format("{:.2f}%"), width="stretch")
    st.dataframe(shown, width="stretch", hide_index=True, height=400,
                 column_config={"percentage": st.column_config.NumberColumn(format="%.2f")})
    st.download_button("Download table (CSV)", shown.to_csv(index=False),
                       "cell_frequencies.csv", "text/csv")

# --- Part 3 ---------------------------------------------------------------
with tab_resp:
    cohort, stats = responder_cohort()
    n = cohort.drop_duplicates("subject_id").response.value_counts()
    st.subheader("Melanoma patients on miraclib, PBMC samples")
    st.caption(f"{cohort['sample'].nunique():,} samples from {n.sum()} patients "
               f"({n['yes']} responders, {n['no']} non-responders), all timepoints.")

    p_by_pop = stats.set_index("population").p_value
    plot = cohort.assign(
        group=cohort.response.map(RESPONSE_LABELS),
        label=cohort.population.map(lambda p: f"{p}<br>p = {p_by_pop[p]:.3g}"),
    )
    fig = px.box(
        plot, x="label", y="percentage", color="group",
        category_orders={"label": [f"{p}<br>p = {p_by_pop[p]:.3g}" for p in analysis.POPULATIONS],
                         "group": list(RESPONSE_LABELS.values())},
        color_discrete_map={RESPONSE_LABELS[k]: v for k, v in RESPONSE_COLORS.items()},
        labels={"label": "", "percentage": "Relative frequency (%)", "group": ""},
        hover_data={"sample": True, "time_from_treatment_start": True, "label": False},
    )
    fig.update_layout(legend=dict(orientation="h", y=1.08), margin=dict(t=40), height=480)
    st.plotly_chart(fig)

    st.markdown("**Mann-Whitney U test per population** (two-sided), "
                "Benjamini-Hochberg correction across the 5 populations.")
    st.dataframe(
        stats, width="stretch", hide_index=True,
        column_config={c: st.column_config.NumberColumn(format="%.4f") for c in stats.columns
                       if c.startswith(("median", "effect", "p_value", "q_value"))},
    )
    st.info(analysis.summarize(stats))
    st.markdown(
        "- **p_value_subject_mean** repeats the test with one value per patient (mean of days 0/7/14), "
        "since each patient contributes 3 samples.\n"
        "- **p_value_baseline** uses day-0 samples only: the timepoint that could predict response "
        "before treatment. No population differs at baseline, so the CD4 difference emerges during "
        "treatment rather than predicting it up front.")

# --- Part 4 ---------------------------------------------------------------
with tab_subset:
    st.subheader("Sample subset explorer")
    st.caption("Defaults: melanoma PBMC samples at baseline (day 0) from patients treated with miraclib.")
    conditions = distinct("condition", "subjects")
    treatments = distinct("treatment", "subjects")
    sample_types = distinct("sample_type", "samples")
    times = distinct("time_from_treatment_start", "samples")
    c1, c2, c3, c4 = st.columns(4)
    condition = c1.selectbox("Condition", conditions, default_index(conditions, "melanoma"))
    treatment = c2.selectbox("Treatment", treatments, default_index(treatments, "miraclib"))
    sample_type = c3.selectbox("Sample type", sample_types, default_index(sample_types, "PBMC"))
    time = c4.selectbox("Days from treatment start", times, default_index(times, 0))

    sub, counts = subset(condition, treatment, sample_type, time)
    if sub.empty:
        st.warning("No samples match this combination.")
    else:
        m1, m2 = st.columns(2)
        m1.metric("Samples", f"{len(sub):,}")
        m2.metric("Subjects", f"{sub.subject.nunique():,}")
        k1, k2, k3 = st.columns(3)
        k1.markdown("**Samples per project**")
        k1.dataframe(counts["samples_per_project"], hide_index=True, width="stretch")
        k2.markdown("**Subjects by response**")
        k2.dataframe(counts["subjects_per_response"].fillna("n/a (untreated)"),
                     hide_index=True, width="stretch")
        k3.markdown("**Subjects by sex**")
        k3.dataframe(counts["subjects_per_sex"], hide_index=True, width="stretch")
        st.markdown("**Matching samples**")
        st.dataframe(sub, hide_index=True, width="stretch", height=360)
        st.download_button("Download subset (CSV)", sub.to_csv(index=False),
                           "subset.csv", "text/csv")
