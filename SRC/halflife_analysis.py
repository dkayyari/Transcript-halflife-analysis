"""
Transcript Half-Life Analysis from mRNA Decay Timecourse Data
===============================================================

Estimates the half-life of each yeast (S. cerevisiae) mRNA transcript from a
60-minute decay time series, then ranks genes to identify the most and least
stable transcripts for downstream GO enrichment analysis.

Method
------
- Missing time points are filled in using linear interpolation in log-space,
  so all 9 time points contribute equally to the decay fit.
- A linear regression of ln(expression) vs. time, forced through the origin,
  is fit per replicate to estimate the decay rate (lambda).
  (Forcing through the origin is valid here because expression is normalized
  to 1 at t=0, and ln(1) = 0.)
- Half-life = ln(2) / lambda.
- The three replicates are averaged per gene to produce a final estimate.

Usage
-----
Run from anywhere inside the repo; paths are resolved relative to this file:

    python src/halflife_analysis.py

Outputs are written to ../results/ (all_halflives.csv, bottom10pct_genes.txt,
top10pct_genes.txt).
"""

import os
import warnings

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────
# PATHS  (relative to the repo root)
# ─────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)

INPUT_FILE = os.path.join(REPO_ROOT, "data", "DecayTimecourse.txt")
OUTPUT_DIR = os.path.join(REPO_ROOT, "results")

OUT_ALL = os.path.join(OUTPUT_DIR, "all_halflives.csv")
OUT_BOTTOM = os.path.join(OUTPUT_DIR, "bottom10pct_genes.txt")
OUT_TOP = os.path.join(OUTPUT_DIR, "top10pct_genes.txt")

TIME_POINTS = np.array([0, 5, 10, 15, 20, 30, 40, 50, 60], dtype=float)

# Columns: 0 = gene, 1-9 = rep1, 10-18 = rep2, 19-27 = rep3
REP_COL_RANGES = {
    "rep1": list(range(1, 10)),
    "rep2": list(range(10, 19)),
    "rep3": list(range(19, 28)),
}

TOP_BOTTOM_PCT = 0.10  # top/bottom 10%


# ─────────────────────────────────────────────
# INTERPOLATION
# ─────────────────────────────────────────────
def interpolate_timeseries(times, values):
    """Fill missing/non-positive time points via linear interpolation in log-space.

    Needs at least 2 valid (non-NaN, positive) observations; otherwise the
    series is returned unchanged and will be filtered out downstream.
    """
    values = values.copy().astype(float)
    valid_mask = ~np.isnan(values) & (values > 0)

    if valid_mask.sum() < 2:
        return values

    t_valid = times[valid_mask]
    ln_v_valid = np.log(values[valid_mask])

    interpolator = interp1d(
        t_valid, ln_v_valid,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
    )

    for i, t in enumerate(times):
        if np.isnan(values[i]) or values[i] <= 0:
            values[i] = np.exp(interpolator(t))

    return values


# ─────────────────────────────────────────────
# HALF-LIFE CALCULATION
# ─────────────────────────────────────────────
def calc_halflife(times, values):
    """Fit ln(expression) vs. time (through the origin) and return the half-life.

    Requires at least 3 valid points. A non-decaying series (slope >= 0)
    has an undefined half-life and returns NaN.
    """
    mask = ~np.isnan(values) & (values > 0)
    t, v = times[mask], values[mask]

    if len(t) < 3:
        return np.nan

    ln_v = np.log(v)
    slope = np.sum(t * ln_v) / np.sum(t ** 2)

    if slope >= 0:
        return np.nan

    return np.log(2) / -slope


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 60)
    print("TRANSCRIPT HALF-LIFE ANALYSIS  (with interpolation)")
    print("=" * 60)

    raw = pd.read_csv(INPUT_FILE, sep="\t", header=None)

    # Row 0: timecourse labels, Row 1: time points, Row 2+: gene data
    gene_data = raw.iloc[2:].copy()
    gene_data.columns = range(gene_data.shape[1])

    results = []
    for _, row in gene_data.iterrows():
        gene = row[0]
        if pd.isna(gene) or str(gene).strip() == "":
            continue

        half_lives = []
        for cols in REP_COL_RANGES.values():
            raw_vals = pd.to_numeric(row[cols].values, errors="coerce").astype(float)
            filled_vals = interpolate_timeseries(TIME_POINTS, raw_vals)
            half_lives.append(calc_halflife(TIME_POINTS, filled_vals))

        valid_hls = [h for h in half_lives if not np.isnan(h)]
        avg_hl = np.mean(valid_hls) if valid_hls else np.nan

        results.append({
            "gene": gene,
            "halflife_rep1": half_lives[0],
            "halflife_rep2": half_lives[1],
            "halflife_rep3": half_lives[2],
            "avg_halflife": avg_hl,
            "n_valid_reps": len(valid_hls),
        })

    df = pd.DataFrame(results)
    df_valid = df.dropna(subset=["avg_halflife"]).copy()
    df_valid = df_valid[df_valid["avg_halflife"] > 0].copy()
    df_valid = df_valid.sort_values("avg_halflife").reset_index(drop=True)

    print(f"\nTotal genes in file         : {len(df)}")
    print(f"Genes with valid half-lives : {len(df_valid)}")
    print("\nHalf-life summary (minutes):")
    print(df_valid["avg_halflife"].describe().round(2))

    # Top / bottom 10%
    n10 = int(np.ceil(len(df_valid) * TOP_BOTTOM_PCT))
    bottom10 = df_valid.head(n10).copy()
    top10 = df_valid.tail(n10).copy()

    print(f"\nBottom {TOP_BOTTOM_PCT:.0%} (shortest half-lives) — {len(bottom10)} genes")
    print(f"Top {TOP_BOTTOM_PCT:.0%} (longest half-lives) — {len(top10)} genes")

    # Save results
    df_valid.to_csv(OUT_ALL, index=False)
    bottom10["gene"].to_csv(OUT_BOTTOM, index=False, header=False)
    top10["gene"].to_csv(OUT_TOP, index=False, header=False)

    print(f"\nSaved: {OUT_ALL}")
    print(f"Saved: {OUT_BOTTOM}")
    print(f"Saved: {OUT_TOP}")


if __name__ == "__main__":
    main()
