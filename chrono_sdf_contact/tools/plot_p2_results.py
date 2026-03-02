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
    m = sum(values) / float(n)
    if n < 2:
        return m, 0.0
    var = sum((x - m) ** 2 for x in values) / float(n - 1)
    std = math.sqrt(max(var, 0.0))
    return m, 1.96 * std / math.sqrt(float(n))


def aggregate(rows, group_keys, metric_keys):
    grouped = defaultdict(list)
    for r in rows:
        grouped[tuple(r.get(k, "") for k in group_keys)].append(r)

    out = []
    for key, grows in grouped.items():
        row = {k: v for k, v in zip(group_keys, key)}
        row["n"] = len(grows)
        for mk in metric_keys:
            vals = [to_float(r, mk, 0.0) for r in grows]
            m, ci = mean_ci95(vals)
            row[f"{mk}_mean"] = m
            row[f"{mk}_ci95"] = ci
        out.append(row)
    return out


def plot_scene_matrix(scene_rows, out_dir):
    import matplotlib.pyplot as plt

    agg = aggregate(
        scene_rows,
        group_keys=["scene", "mode", "step_size"],
        metric_keys=["max_penetration", "wall_time_s", "exploded", "max_relative_energy_drift"],
    )
    scenes = sorted(set(r["scene"] for r in agg))
    if not scenes:
        return

    ncols = 2
    nrows = int(math.ceil(len(scenes) / float(ncols)))
    fig_pen, axes_pen = plt.subplots(nrows, ncols, figsize=(12, 4.2 * nrows))
    fig_wall, axes_wall = plt.subplots(nrows, ncols, figsize=(12, 4.2 * nrows))
    axes_pen = axes_pen.flatten() if hasattr(axes_pen, "flatten") else [axes_pen]
    axes_wall = axes_wall.flatten() if hasattr(axes_wall, "flatten") else [axes_wall]

    mode_order = ["chrono_baseline", "sdf_contact"]

    for i, scene in enumerate(scenes):
        ax_pen = axes_pen[i]
        ax_wall = axes_wall[i]
        scene_rows_agg = [r for r in agg if r["scene"] == scene]
        for mode in mode_order:
            mode_rows = [r for r in scene_rows_agg if r["mode"] == mode]
            if not mode_rows:
                continue
            mode_rows.sort(key=lambda r: to_float(r, "step_size"))
            dts = [to_float(r, "step_size") for r in mode_rows]
            pen = [to_float(r, "max_penetration_mean") for r in mode_rows]
            pen_ci = [to_float(r, "max_penetration_ci95") for r in mode_rows]
            wall = [to_float(r, "wall_time_s_mean") for r in mode_rows]
            wall_ci = [to_float(r, "wall_time_s_ci95") for r in mode_rows]
            unstable = [to_float(r, "exploded_mean") for r in mode_rows]

            ax_pen.errorbar(dts, pen, yerr=pen_ci, marker="o", capsize=3, label=mode)
            ax_wall.errorbar(dts, wall, yerr=wall_ci, marker="o", capsize=3, label=mode)
            for j, u in enumerate(unstable):
                if u > 0:
                    ax_pen.annotate(f"{u * 100.0:.0f}%", (dts[j], pen[j]), color="red", fontsize=8)

        ax_pen.set_title(f"{scene}: penetration vs dt")
        ax_pen.set_xlabel("dt")
        ax_pen.set_ylabel("max_penetration")
        ax_pen.grid(True, alpha=0.3)
        ax_pen.legend()

        ax_wall.set_title(f"{scene}: wall-time vs dt")
        ax_wall.set_xlabel("dt")
        ax_wall.set_ylabel("wall_time_s")
        ax_wall.grid(True, alpha=0.3)
        ax_wall.legend()

    for i in range(len(scenes), len(axes_pen)):
        axes_pen[i].axis("off")
        axes_wall[i].axis("off")

    fig_pen.tight_layout()
    fig_wall.tight_layout()
    pen_path = os.path.join(out_dir, "p2_scene_penetration.png")
    wall_path = os.path.join(out_dir, "p2_scene_walltime.png")
    fig_pen.savefig(pen_path, dpi=160)
    fig_wall.savefig(wall_path, dpi=160)
    plt.close(fig_pen)
    plt.close(fig_wall)
    print(f"Generated: {pen_path}")
    print(f"Generated: {wall_path}")


def plot_scaleup(scale_rows, out_dir):
    import matplotlib.pyplot as plt

    agg = aggregate(
        scale_rows,
        group_keys=["mode", "body_count"],
        metric_keys=["wall_time_s", "avg_contacts", "max_penetration", "max_relative_energy_drift", "exploded"],
    )
    modes = sorted(set(r["mode"] for r in agg))
    if not modes:
        return

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

    for mode in modes:
        mode_rows = [r for r in agg if r["mode"] == mode]
        mode_rows.sort(key=lambda r: to_float(r, "body_count"))
        counts = [to_float(r, "body_count") for r in mode_rows]
        wall = [to_float(r, "wall_time_s_mean") for r in mode_rows]
        wall_ci = [to_float(r, "wall_time_s_ci95") for r in mode_rows]
        contacts = [to_float(r, "avg_contacts_mean") for r in mode_rows]
        contacts_ci = [to_float(r, "avg_contacts_ci95") for r in mode_rows]
        pen = [to_float(r, "max_penetration_mean") for r in mode_rows]
        pen_ci = [to_float(r, "max_penetration_ci95") for r in mode_rows]
        unstable = [to_float(r, "exploded_mean") for r in mode_rows]

        axes[0].errorbar(counts, wall, yerr=wall_ci, marker="o", capsize=3, label=mode)
        axes[1].errorbar(counts, contacts, yerr=contacts_ci, marker="o", capsize=3, label=mode)
        axes[2].errorbar(counts, pen, yerr=pen_ci, marker="o", capsize=3, label=mode)
        for i, u in enumerate(unstable):
            if u > 0:
                axes[2].annotate(f"{u * 100.0:.0f}%", (counts[i], pen[i]), color="red", fontsize=8)

    axes[0].set_title("Scale-up wall time")
    axes[0].set_xlabel("body_count")
    axes[0].set_ylabel("wall_time_s")
    axes[0].grid(True, alpha=0.3)

    axes[1].set_title("Scale-up avg contacts")
    axes[1].set_xlabel("body_count")
    axes[1].set_ylabel("avg_contacts")
    axes[1].grid(True, alpha=0.3)

    axes[2].set_title("Scale-up penetration")
    axes[2].set_xlabel("body_count")
    axes[2].set_ylabel("max_penetration")
    axes[2].grid(True, alpha=0.3)

    for ax in axes:
        ax.legend()

    fig.tight_layout()
    out_path = os.path.join(out_dir, "p2_scaleup_summary.png")
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    print(f"Generated: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Plot P2 benchmark CSV outputs.")
    parser.add_argument("--scene-csv", required=True, help="Path to sdf_p2_scene_matrix.csv")
    parser.add_argument("--scale-csv", required=True, help="Path to sdf_p2_scaleup.csv")
    parser.add_argument("--out-dir", default=".", help="Output directory")
    args = parser.parse_args()

    try:
        import matplotlib  # noqa: F401
    except Exception:
        print("matplotlib is required. Install via: pip install matplotlib", file=sys.stderr)
        return 2

    os.makedirs(args.out_dir, exist_ok=True)
    scene_rows = read_rows(args.scene_csv)
    scale_rows = read_rows(args.scale_csv)
    if not scene_rows:
        print("No rows in scene CSV.", file=sys.stderr)
        return 3
    if not scale_rows:
        print("No rows in scale CSV.", file=sys.stderr)
        return 3

    plot_scene_matrix(scene_rows, args.out_dir)
    plot_scaleup(scale_rows, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
