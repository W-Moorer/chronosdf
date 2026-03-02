#!/usr/bin/env python3
import argparse
import csv
import math
import os
import sys
from collections import defaultdict


def read_csv(path):
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
    ci95 = 1.96 * std / math.sqrt(float(n))
    return mean, ci95


def aggregate_rows(rows, group_keys, metric_keys):
    groups = defaultdict(list)
    for r in rows:
        key = tuple(r.get(k, "") for k in group_keys)
        groups[key].append(r)

    out = []
    for key, grows in groups.items():
        row = {k: v for k, v in zip(group_keys, key)}
        row["n"] = len(grows)
        for mk in metric_keys:
            vals = [to_float(r, mk, 0.0) for r in grows]
            mean, ci95 = mean_ci95(vals)
            row[f"{mk}_mean"] = mean
            row[f"{mk}_ci95"] = ci95
        out.append(row)
    return out


def plot_mesh_compare(rows, out_path):
    agg = aggregate_rows(
        rows,
        group_keys=["mode"],
        metric_keys=["wall_time_s", "avg_contacts", "max_penetration", "max_relative_energy_drift"],
    )
    agg.sort(key=lambda r: r["mode"])
    modes = [r["mode"] for r in agg]
    wall = [r["wall_time_s_mean"] for r in agg]
    wall_ci = [r["wall_time_s_ci95"] for r in agg]
    avg_contacts = [r["avg_contacts_mean"] for r in agg]
    avg_contacts_ci = [r["avg_contacts_ci95"] for r in agg]
    penetration = [r["max_penetration_mean"] for r in agg]
    penetration_ci = [r["max_penetration_ci95"] for r in agg]
    drift = [r["max_relative_energy_drift_mean"] for r in agg]
    drift_ci = [r["max_relative_energy_drift_ci95"] for r in agg]

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 7.2))
    axes = axes.flatten()

    axes[0].bar(modes, wall, yerr=wall_ci, capsize=3)
    axes[0].set_title("Wall Time (s)")
    axes[0].set_ylabel("seconds")
    axes[0].tick_params(axis="x", rotation=20)

    axes[1].bar(modes, avg_contacts, yerr=avg_contacts_ci, capsize=3)
    axes[1].set_title("Average Contacts")
    axes[1].tick_params(axis="x", rotation=20)

    axes[2].bar(modes, penetration, yerr=penetration_ci, capsize=3)
    axes[2].set_title("Max Penetration")
    axes[2].tick_params(axis="x", rotation=20)

    axes[3].bar(modes, drift, yerr=drift_ci, capsize=3)
    axes[3].set_title("Max Relative Energy Drift")
    axes[3].tick_params(axis="x", rotation=20)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_sdf_sdf_stability(rows, out_path):
    import matplotlib.pyplot as plt

    agg = aggregate_rows(
        rows,
        group_keys=["mode", "step_size"],
        metric_keys=["wall_time_s", "max_penetration", "exploded", "max_relative_energy_drift"],
    )
    grouped = defaultdict(list)
    for r in agg:
        grouped[r["mode"]].append(r)
    for mode_rows in grouped.values():
        mode_rows.sort(key=lambda r: to_float(r, "step_size"))

    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))

    for mode, mode_rows in grouped.items():
        dts = [to_float(r, "step_size") for r in mode_rows]
        wall = [to_float(r, "wall_time_s_mean") for r in mode_rows]
        wall_ci = [to_float(r, "wall_time_s_ci95") for r in mode_rows]
        pen = [to_float(r, "max_penetration_mean") for r in mode_rows]
        pen_ci = [to_float(r, "max_penetration_ci95") for r in mode_rows]
        drift = [to_float(r, "max_relative_energy_drift_mean") for r in mode_rows]
        drift_ci = [to_float(r, "max_relative_energy_drift_ci95") for r in mode_rows]
        exploded_rate = [to_float(r, "exploded_mean") for r in mode_rows]

        axes[0].errorbar(dts, wall, yerr=wall_ci, marker="o", capsize=3, label=mode)
        axes[1].errorbar(dts, pen, yerr=pen_ci, marker="o", capsize=3, label=mode)
        axes[2].errorbar(dts, drift, yerr=drift_ci, marker="o", capsize=3, label=mode)

        for i, rate in enumerate(exploded_rate):
            if rate > 0.0:
                axes[1].annotate(f"{rate * 100.0:.0f}%", (dts[i], pen[i]), color="red", fontsize=8, ha="center", va="bottom")

    axes[0].set_title("Wall Time vs dt")
    axes[0].set_xlabel("dt")
    axes[0].set_ylabel("seconds")
    axes[0].grid(True, alpha=0.3)

    axes[1].set_title("Max Penetration vs dt")
    axes[1].set_xlabel("dt")
    axes[1].set_ylabel("penetration")
    axes[1].grid(True, alpha=0.3)

    axes[2].set_title("Max Relative Energy Drift vs dt")
    axes[2].set_xlabel("dt")
    axes[2].set_ylabel("drift")
    axes[2].grid(True, alpha=0.3)

    axes[0].legend()
    axes[1].legend()
    axes[2].legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot chrono_sdf_contact experiment CSV outputs.")
    parser.add_argument("--mesh-csv", required=True, help="Path to sdf_mesh_baseline_compare.csv")
    parser.add_argument("--stability-csv", required=True, help="Path to sdf_sdf_stability_compare.csv")
    parser.add_argument("--out-dir", default=".", help="Output directory for generated PNG files")
    args = parser.parse_args()

    try:
        import matplotlib  # noqa: F401
    except Exception:
        print("matplotlib is required. Install via: pip install matplotlib", file=sys.stderr)
        return 2

    os.makedirs(args.out_dir, exist_ok=True)

    mesh_rows = read_csv(args.mesh_csv)
    stability_rows = read_csv(args.stability_csv)

    mesh_png = os.path.join(args.out_dir, "mesh_compare_summary.png")
    sdf_sdf_png = os.path.join(args.out_dir, "sdf_sdf_stability_summary.png")

    plot_mesh_compare(mesh_rows, mesh_png)
    plot_sdf_sdf_stability(stability_rows, sdf_sdf_png)

    print(f"Generated: {mesh_png}")
    print(f"Generated: {sdf_sdf_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
