#!/usr/bin/env python3
import argparse
import csv
import math
import os
import sys
from collections import defaultdict


def read_rows(path):
    rows = []
    with open(path, "r", newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def to_float(row, key, default=0.0):
    try:
        return float(row.get(key, default))
    except Exception:
        return default


def mean_ci95(values):
    if not values:
        return 0.0, 0.0
    n = len(values)
    mean = sum(values) / float(n)
    if n < 2:
        return mean, 0.0
    var = sum((x - mean) ** 2 for x in values) / float(n - 1)
    std = math.sqrt(max(var, 0.0))
    return mean, 1.96 * std / math.sqrt(float(n))


def aggregate_by_eps(rows):
    grouped = defaultdict(list)
    for r in rows:
        eps = to_float(r, "eps_force")
        grouped[eps].append(r)

    agg = []
    for eps, grows in grouped.items():
        qvals = [to_float(r, "query_time_s") for r in grows]
        rvals = [to_float(r, "rmse_phi") for r in grows]
        q_mean, q_ci = mean_ci95(qvals)
        r_mean, r_ci = mean_ci95(rvals)
        agg.append(
            {
                "eps_force": eps,
                "n": len(grows),
                "query_time_s_mean": q_mean,
                "query_time_s_ci95": q_ci,
                "rmse_phi_mean": r_mean,
                "rmse_phi_ci95": r_ci,
            }
        )
    agg.sort(key=lambda r: r["eps_force"], reverse=True)
    return agg


def extract_pareto(points):
    front = []
    for i, p in enumerate(points):
        dominated = False
        for j, q in enumerate(points):
            if i == j:
                continue
            better_or_equal_time = q["query_time_s_mean"] <= p["query_time_s_mean"]
            better_or_equal_err = q["rmse_phi_mean"] <= p["rmse_phi_mean"]
            strictly_better = (q["query_time_s_mean"] < p["query_time_s_mean"]) or (
                q["rmse_phi_mean"] < p["rmse_phi_mean"]
            )
            if better_or_equal_time and better_or_equal_err and strictly_better:
                dominated = True
                break
        if not dominated:
            front.append(p)
    front.sort(key=lambda r: r["query_time_s_mean"])
    return front


def main():
    parser = argparse.ArgumentParser(description="Plot M3 octree Pareto curves.")
    parser.add_argument("--csv", required=True, help="Path to sdf_octree_pareto.csv")
    parser.add_argument("--pareto-csv", default="", help="Path to sdf_octree_pareto_front.csv")
    parser.add_argument("--out", default="m3_octree_pareto.png", help="Output PNG")
    args = parser.parse_args()

    try:
        import matplotlib.pyplot as plt
    except Exception:
        print("matplotlib is required. Install via: pip install matplotlib", file=sys.stderr)
        return 2

    rows = read_rows(args.csv)
    if not rows:
        print("No rows in input CSV.", file=sys.stderr)
        return 3

    agg = aggregate_by_eps(rows)
    if not agg:
        print("No valid aggregated points in input CSV.", file=sys.stderr)
        return 3

    fig, ax = plt.subplots(1, 1, figsize=(6.8, 4.8))
    x = [r["query_time_s_mean"] for r in agg]
    y = [r["rmse_phi_mean"] for r in agg]
    xerr = [r["query_time_s_ci95"] for r in agg]
    yerr = [r["rmse_phi_ci95"] for r in agg]
    labels = [r["eps_force"] for r in agg]

    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt="o", color="steelblue", ecolor="lightsteelblue", capsize=3, label="mean +/- 95% CI")
    for xi, yi, lb, row in zip(x, y, labels, agg):
        ax.annotate(
            f"eps={lb:g} (n={row['n']})",
            (xi, yi),
            fontsize=8,
            xytext=(4, 4),
            textcoords="offset points",
        )

    pareto = extract_pareto(agg)
    if pareto:
        px = [r["query_time_s_mean"] for r in pareto]
        py = [r["rmse_phi_mean"] for r in pareto]
        ax.plot(px, py, "-o", color="crimson", label="pareto front (mean)", markersize=4)

    ax.set_title("M3 Octree Pareto (time vs RMSE_phi)")
    ax.set_xlabel("query_time_s")
    ax.set_ylabel("rmse_phi")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(args.out, dpi=160)
    plt.close(fig)

    print(f"Generated: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
