#!/usr/bin/env python3
"""Plot W3 ablation outputs (kmax/cell, clustering, narrowphase scene exploded rates)."""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from typing import Dict, List, Tuple


def normalize_key(key: str) -> str:
    # Handle UTF-8 BOM and accidental quoted headers from mixed CSV producers.
    return key.replace("\ufeff", "").strip().strip('"').strip("'")


def read_rows(path: str) -> List[dict]:
    rows: List[dict] = []
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            clean: dict = {}
            for k, v in row.items():
                if not k:
                    continue
                nk = normalize_key(k)
                clean[nk] = v.strip() if isinstance(v, str) else v
            rows.append(clean)
    return rows


def to_float(row: dict, key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, default))
    except Exception:
        return default


def plot_kmax_cell(rows: List[dict], out_dir: str) -> str:
    import matplotlib.pyplot as plt

    kmaxs = sorted({int(to_float(r, "kmax")) for r in rows})
    cells = sorted({to_float(r, "cell_size") for r in rows})
    if not kmaxs or not cells:
        raise RuntimeError("kmax/cell summary is empty")

    idx_k = {k: i for i, k in enumerate(kmaxs)}
    idx_c = {c: i for i, c in enumerate(cells)}
    exploded = [[math.nan for _ in kmaxs] for _ in cells]
    overhead = [[math.nan for _ in kmaxs] for _ in cells]

    for r in rows:
        k = int(to_float(r, "kmax"))
        c = to_float(r, "cell_size")
        i = idx_c[c]
        j = idx_k[k]
        exploded[i][j] = to_float(r, "sdf_exploded_rate")
        overhead[i][j] = to_float(r, "overhead_ratio_sdf_over_baseline")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    maps = [
        (exploded, "SDF exploded_rate", "Reds", lambda v: f"{100.0 * v:.0f}%"),
        (overhead, "SDF/Baseline overhead", "Blues", lambda v: f"{v:.2f}x"),
    ]

    for ax, (mat, title, cmap, fmt) in zip(axes, maps):
        im = ax.imshow(mat, origin="lower", aspect="auto", cmap=cmap)
        ax.set_title(title)
        ax.set_xlabel("kmax")
        ax.set_ylabel("cell_size")
        ax.set_xticks(range(len(kmaxs)))
        ax.set_xticklabels([str(k) for k in kmaxs])
        ax.set_yticks(range(len(cells)))
        ax.set_yticklabels([f"{c:.3f}" for c in cells])
        for i in range(len(cells)):
            for j in range(len(kmaxs)):
                val = mat[i][j]
                if math.isnan(val):
                    continue
                ax.text(j, i, fmt(val), ha="center", va="center", fontsize=8, color="black")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.tight_layout()
    out_path = os.path.join(out_dir, "w3_stack3_kmax_cell_heatmaps.png")
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    return out_path


def plot_cluster(rows: List[dict], out_dir: str) -> str:
    import matplotlib.pyplot as plt

    kmaxs = sorted({int(to_float(r, "kmax")) for r in rows})
    if not kmaxs:
        raise RuntimeError("cluster summary is empty")

    modes = ["on", "off"]
    exploded_vals: Dict[Tuple[int, str], float] = {}
    overhead_vals: Dict[Tuple[int, str], float] = {}
    for r in rows:
        k = int(to_float(r, "kmax"))
        mode = r.get("clustering_mode", "")
        exploded_vals[(k, mode)] = to_float(r, "sdf_exploded_rate")
        overhead_vals[(k, mode)] = to_float(r, "overhead_ratio_sdf_over_baseline")

    x = list(range(len(kmaxs)))
    width = 0.35
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))

    for m_idx, mode in enumerate(modes):
        shift = (m_idx - 0.5) * width
        e = [exploded_vals.get((k, mode), math.nan) for k in kmaxs]
        o = [overhead_vals.get((k, mode), math.nan) for k in kmaxs]
        axes[0].bar([v + shift for v in x], e, width=width, label=mode)
        axes[1].bar([v + shift for v in x], o, width=width, label=mode)
        # Keep the "all zero" case visible to avoid appearing as an empty chart.
        axes[0].plot([v + shift for v in x], e, marker="o", linewidth=1.2)

    axes[0].set_title("Stack3 exploded_rate vs clustering")
    axes[0].set_xlabel("kmax")
    axes[0].set_ylabel("exploded_rate")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([str(k) for k in kmaxs])
    axes[0].set_ylim(0.0, 1.0)
    axes[0].grid(True, axis="y", alpha=0.3)
    axes[0].legend(title="cluster")

    axes[1].set_title("Stack3 overhead vs clustering")
    axes[1].set_xlabel("kmax")
    axes[1].set_ylabel("sdf/baseline")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([str(k) for k in kmaxs])
    axes[1].grid(True, axis="y", alpha=0.3)
    axes[1].legend(title="cluster")

    fig.tight_layout()
    out_path = os.path.join(out_dir, "w3_stack3_cluster_ablation.png")
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    return out_path


def plot_np_scene(rows: List[dict], out_dir: str) -> str:
    import matplotlib.pyplot as plt

    scenes = ["box_plane_rest", "box_plane_impact", "nonconvex_t_drop", "stack3_boxes"]
    np_modes = ["default", "no_cull", "no_cache", "no_cache_no_cull"]
    kmaxs = [64, 128]
    cell = "0.05"

    exploded: Dict[Tuple[int, str, str], float] = {}
    for r in rows:
        if r.get("contact_mode") != "sdf_contact":
            continue
        if r.get("clustering_mode") != "on":
            continue
        if r.get("cell_size") != cell:
            continue
        try:
            k = int(float(r.get("kmax", "0")))
        except Exception:
            continue
        np_mode = r.get("narrowphase_mode", "")
        scene = r.get("scene", "")
        exploded[(k, np_mode, scene)] = to_float(r, "exploded_rate")

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=True)
    width = 0.18
    x = list(range(len(scenes)))

    for ax_idx, k in enumerate(kmaxs):
        ax = axes[ax_idx]
        for m_idx, np_mode in enumerate(np_modes):
            shift = (m_idx - 1.5) * width
            y = [exploded.get((k, np_mode, scene), math.nan) for scene in scenes]
            xpos = [v + shift for v in x]
            bars = ax.bar(xpos, y, width=width, label=np_mode)
            # Explicitly mark zero values so "all-zero" scenes do not look like missing bars.
            face = bars[0].get_facecolor() if len(bars) > 0 else None
            zero_x = [xx for xx, yy in zip(xpos, y) if (not math.isnan(yy)) and abs(yy) < 1e-12]
            if zero_x:
                ax.plot(
                    zero_x,
                    [0.0] * len(zero_x),
                    linestyle="None",
                    marker="o",
                    markersize=3.0,
                    color=face,
                    zorder=3,
                )
        ax.set_title(f"SDF exploded_rate by scene (kmax={k})")
        ax.set_xticks(x)
        ax.set_xticklabels(scenes, rotation=15, ha="right")
        ax.set_ylim(0.0, 1.0)
        ax.grid(True, axis="y", alpha=0.3)
        ax.set_ylabel("exploded_rate")
        if ax_idx == 1:
            ax.legend(title="narrowphase", fontsize=8)

    fig.tight_layout()
    out_path = os.path.join(out_dir, "w3_narrowphase_scene_exploded.png")
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--kmax-cell-csv",
        default="build/chrono_sdf_contact_vs/Release/results_w3_stack3_sweep_20seeds_v1/w3_stack3_dt003_summary.csv",
    )
    parser.add_argument(
        "--cluster-csv",
        default="build/chrono_sdf_contact_vs/Release/results_w3_stack3_ablation_20seeds_v1/w3_stack3_dt003_summary.csv",
    )
    parser.add_argument(
        "--np-scene-csv",
        default="build/chrono_sdf_contact_vs/Release/results_w3_np_ablation_20seeds_v1/w3_np_scene_metrics.csv",
    )
    parser.add_argument(
        "--out-dir",
        default="build/chrono_sdf_contact_vs/Release/results_w3_np_ablation_20seeds_v1",
    )
    args = parser.parse_args()

    try:
        import matplotlib  # noqa: F401
    except Exception:
        print("matplotlib is required: pip install matplotlib", file=sys.stderr)
        return 2

    os.makedirs(args.out_dir, exist_ok=True)
    kmax_rows = read_rows(args.kmax_cell_csv)
    cluster_rows = read_rows(args.cluster_csv)
    np_rows = read_rows(args.np_scene_csv)
    if not kmax_rows or not cluster_rows or not np_rows:
        print("One or more input CSV files are empty.", file=sys.stderr)
        return 3

    out1 = plot_kmax_cell(kmax_rows, args.out_dir)
    out2 = plot_cluster(cluster_rows, args.out_dir)
    out3 = plot_np_scene(np_rows, args.out_dir)
    print(f"Generated: {out1}")
    print(f"Generated: {out2}")
    print(f"Generated: {out3}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
