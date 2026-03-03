#!/usr/bin/env python3
"""Analyze and plot COM trajectories for baseline-vs-SDF comparisons."""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


def normalize_key(key: str) -> str:
    return key.replace("\ufeff", "").strip().strip('"').strip("'")


def read_rows(path: str) -> List[dict]:
    rows: List[dict] = []
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        for raw in csv.DictReader(f):
            clean: dict = {}
            for k, v in raw.items():
                if not k:
                    continue
                nk = normalize_key(k)
                clean[nk] = v.strip() if isinstance(v, str) else v
            rows.append(clean)
    return rows


def to_float(v, default=0.0) -> float:
    try:
        return float(v)
    except Exception:
        return default


def to_int(v, default=0) -> int:
    try:
        return int(float(v))
    except Exception:
        return default


def fmt(x: float | None, digits: int = 6) -> str:
    if x is None:
        return "-"
    if not math.isfinite(x):
        return "-"
    return f"{x:.{digits}g}"


def match_pairs(rows_a: Iterable[dict], rows_b: Iterable[dict], key_fields: Sequence[str]) -> List[Tuple[dict, dict]]:
    map_a: Dict[Tuple, dict] = {}
    for r in rows_a:
        key = tuple(r.get(k, "") for k in key_fields)
        map_a[key] = r
    pairs: List[Tuple[dict, dict]] = []
    for r in rows_b:
        key = tuple(r.get(k, "") for k in key_fields)
        if key in map_a:
            pairs.append((map_a[key], r))
    return pairs


def xyz_metrics(pairs: Sequence[Tuple[dict, dict]]) -> dict:
    if not pairs:
        return {
            "n": 0,
            "rmse_x": None,
            "rmse_y": None,
            "rmse_z": None,
            "rmse_l2": None,
            "max_l2": None,
            "final_l2": None,
        }

    dx2 = []
    dy2 = []
    dz2 = []
    l2 = []
    for a, b in pairs:
        dx = to_float(b.get("com_x")) - to_float(a.get("com_x"))
        dy = to_float(b.get("com_y")) - to_float(a.get("com_y"))
        dz = to_float(b.get("com_z")) - to_float(a.get("com_z"))
        d2 = dx * dx + dy * dy + dz * dz
        dx2.append(dx * dx)
        dy2.append(dy * dy)
        dz2.append(dz * dz)
        l2.append(math.sqrt(d2))

    rmse_x = math.sqrt(sum(dx2) / len(dx2))
    rmse_y = math.sqrt(sum(dy2) / len(dy2))
    rmse_z = math.sqrt(sum(dz2) / len(dz2))
    rmse_l2 = math.sqrt(sum(v * v for v in l2) / len(l2))
    max_l2 = max(l2) if l2 else None
    final_l2 = l2[-1] if l2 else None
    return {
        "n": len(pairs),
        "rmse_x": rmse_x,
        "rmse_y": rmse_y,
        "rmse_z": rmse_z,
        "rmse_l2": rmse_l2,
        "max_l2": max_l2,
        "final_l2": final_l2,
    }


def sorted_rows_by_step(rows: Iterable[dict]) -> List[dict]:
    return sorted(rows, key=lambda r: (to_int(r.get("body_id")), to_int(r.get("step"))))


def plot_xyz_all_bodies(
    baseline_rows: List[dict],
    sdf_rows: List[dict],
    out_path: str,
    title: str,
    body_ids: List[int],
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    axes[0].set_title(title)
    colors = plt.cm.tab20.colors

    for bid in body_ids:
        b_rows = sorted((r for r in baseline_rows if to_int(r.get("body_id")) == bid), key=lambda r: to_int(r.get("step")))
        s_rows = sorted((r for r in sdf_rows if to_int(r.get("body_id")) == bid), key=lambda r: to_int(r.get("step")))
        if not b_rows or not s_rows:
            continue
        color = colors[bid % len(colors)]
        tb = [to_float(r.get("time")) for r in b_rows]
        ts = [to_float(r.get("time")) for r in s_rows]
        xb = [to_float(r.get("com_x")) for r in b_rows]
        yb = [to_float(r.get("com_y")) for r in b_rows]
        zb = [to_float(r.get("com_z")) for r in b_rows]
        xs = [to_float(r.get("com_x")) for r in s_rows]
        ys = [to_float(r.get("com_y")) for r in s_rows]
        zs = [to_float(r.get("com_z")) for r in s_rows]
        axes[0].plot(tb, xb, linestyle="--", linewidth=1.0, alpha=0.65, color=color)
        axes[1].plot(tb, yb, linestyle="--", linewidth=1.0, alpha=0.65, color=color)
        axes[2].plot(tb, zb, linestyle="--", linewidth=1.0, alpha=0.65, color=color)
        axes[0].plot(ts, xs, linestyle="-", linewidth=1.0, alpha=0.9, color=color)
        axes[1].plot(ts, ys, linestyle="-", linewidth=1.0, alpha=0.9, color=color)
        axes[2].plot(ts, zs, linestyle="-", linewidth=1.0, alpha=0.9, color=color)

    axes[0].set_ylabel("COM X")
    axes[1].set_ylabel("COM Y")
    axes[2].set_ylabel("COM Z")
    axes[2].set_xlabel("time (s)")
    for ax in axes:
        ax.grid(True, alpha=0.25)

    handles = [
        Line2D([0], [0], color="black", linestyle="--", linewidth=1.2, label="baseline"),
        Line2D([0], [0], color="black", linestyle="-", linewidth=1.2, label="sdf"),
    ]
    axes[0].legend(handles=handles, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


def analyze_mesh(rows: List[dict], out_dir: str, md: List[str]) -> None:
    baseline = sorted_rows_by_step([r for r in rows if r.get("mode") == "chrono_baseline"])
    sdf = sorted_rows_by_step([r for r in rows if r.get("mode") == "sdf_contact"])
    pairs = match_pairs(baseline, sdf, key_fields=("step", "body_id"))
    m = xyz_metrics(pairs)
    fig = os.path.join(out_dir, "m1_mesh_com_xyz_compare.png")
    plot_xyz_all_bodies(baseline, sdf, fig, "M1 mesh baseline vs sdf (all bodies)", body_ids=[0])

    md.append("## M1 Mesh COM Comparison")
    md.append("")
    md.append("| n_samples | rmse_x | rmse_y | rmse_z | rmse_l2 | max_l2 | final_l2 | figure |")
    md.append("| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
    md.append(
        f"| {m['n']} | {fmt(m['rmse_x'])} | {fmt(m['rmse_y'])} | {fmt(m['rmse_z'])} | {fmt(m['rmse_l2'])} | "
        f"{fmt(m['max_l2'])} | {fmt(m['final_l2'])} | `m1_mesh_com_xyz_compare.png` |"
    )
    md.append("")


def analyze_scene(rows: List[dict], out_dir: str, md: List[str]) -> None:
    md.append("## P2 Scene-Matrix COM Comparison")
    md.append("")
    md.append("| scene | dt | n_samples | rmse_x | rmse_y | rmse_z | rmse_l2 | max_l2 | final_l2 | figure |")
    md.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |")

    keys = sorted({(r.get("scene", ""), r.get("step_size", "")) for r in rows})
    for scene, dt in keys:
        baseline = sorted_rows_by_step(
            [r for r in rows if r.get("scene") == scene and r.get("step_size") == dt and r.get("mode") == "chrono_baseline"]
        )
        sdf = sorted_rows_by_step(
            [r for r in rows if r.get("scene") == scene and r.get("step_size") == dt and r.get("mode") == "sdf_contact"]
        )
        if not baseline or not sdf:
            continue
        pairs = match_pairs(baseline, sdf, key_fields=("body_id", "step"))
        m = xyz_metrics(pairs)
        body_ids = sorted({to_int(r.get("body_id")) for r in baseline + sdf})
        dt_tag = dt.replace(".", "p")
        fig_name = f"p2_scene_{scene}_dt_{dt_tag}_com_xyz.png"
        fig = os.path.join(out_dir, fig_name)
        title = f"P2 scene={scene}, dt={dt} baseline vs sdf (all bodies)"
        plot_xyz_all_bodies(baseline, sdf, fig, title, body_ids)
        md.append(
            f"| {scene} | {dt} | {m['n']} | {fmt(m['rmse_x'])} | {fmt(m['rmse_y'])} | {fmt(m['rmse_z'])} | "
            f"{fmt(m['rmse_l2'])} | {fmt(m['max_l2'])} | {fmt(m['final_l2'])} | `{fig_name}` |"
        )
    md.append("")


def analyze_scale(rows: List[dict], out_dir: str, md: List[str]) -> None:
    md.append("## P2 Scale-up COM Comparison")
    md.append("")
    md.append("| body_count | n_samples | rmse_x | rmse_y | rmse_z | rmse_l2 | max_l2 | final_l2 | figure |")
    md.append("| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |")

    counts = sorted({to_int(r.get("body_count")) for r in rows})
    for count in counts:
        baseline = sorted_rows_by_step(
            [r for r in rows if to_int(r.get("body_count")) == count and r.get("mode") == "chrono_baseline"]
        )
        sdf = sorted_rows_by_step([r for r in rows if to_int(r.get("body_count")) == count and r.get("mode") == "sdf_contact"])
        if not baseline or not sdf:
            continue
        pairs = match_pairs(baseline, sdf, key_fields=("body_id", "step"))
        m = xyz_metrics(pairs)
        body_ids = sorted({to_int(r.get("body_id")) for r in baseline + sdf})
        fig_name = f"p2_scale_count_{count}_com_xyz.png"
        fig = os.path.join(out_dir, fig_name)
        title = f"P2 scale-up body_count={count} baseline vs sdf (all bodies)"
        plot_xyz_all_bodies(baseline, sdf, fig, title, body_ids)
        md.append(
            f"| {count} | {m['n']} | {fmt(m['rmse_x'])} | {fmt(m['rmse_y'])} | {fmt(m['rmse_z'])} | {fmt(m['rmse_l2'])} | "
            f"{fmt(m['max_l2'])} | {fmt(m['final_l2'])} | `{fig_name}` |"
        )
    md.append("")


def analyze_m2(rows: List[dict], out_dir: str, md: List[str]) -> None:
    md.append("## M2 SDF-SDF COM Comparison (One-sided vs Bidirectional)")
    md.append("")
    md.append("| dt | compare | n_samples | rmse_x | rmse_y | rmse_z | rmse_l2 | max_l2 | final_l2 | figure |")
    md.append("| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |")

    dts = sorted({r.get("step_size", "") for r in rows}, key=lambda x: to_float(x))
    for dt in dts:
        bid = sorted_rows_by_step([r for r in rows if r.get("step_size") == dt and r.get("mode") == "bidirectional"])
        if not bid:
            continue
        for ref_mode in ("a_to_b_only", "b_to_a_only"):
            ref = sorted_rows_by_step([r for r in rows if r.get("step_size") == dt and r.get("mode") == ref_mode])
            if not ref:
                continue
            pairs = match_pairs(ref, bid, key_fields=("body_id", "step"))
            m = xyz_metrics(pairs)
            fig_name = f"m2_dt_{dt.replace('.', 'p')}_{ref_mode}_vs_bidirectional_com_xyz.png"
            fig = os.path.join(out_dir, fig_name)
            title = f"M2 dt={dt}: {ref_mode} vs bidirectional"
            plot_xyz_all_bodies(ref, bid, fig, title, body_ids=[0])
            md.append(
                f"| {dt} | {ref_mode} vs bidirectional | {m['n']} | {fmt(m['rmse_x'])} | {fmt(m['rmse_y'])} | "
                f"{fmt(m['rmse_z'])} | {fmt(m['rmse_l2'])} | {fmt(m['max_l2'])} | {fmt(m['final_l2'])} | `{fig_name}` |"
            )
    md.append("")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-traj-csv", required=True)
    parser.add_argument("--scene-traj-csv", required=True)
    parser.add_argument("--scale-traj-csv", required=True)
    parser.add_argument("--m2-traj-csv", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--out-summary", default="com_trajectory_diff_summary.md")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    mesh_rows = read_rows(args.mesh_traj_csv)
    scene_rows = read_rows(args.scene_traj_csv)
    scale_rows = read_rows(args.scale_traj_csv)
    m2_rows = read_rows(args.m2_traj_csv)
    if not mesh_rows or not scene_rows or not scale_rows or not m2_rows:
        raise RuntimeError("Trajectory CSV input is empty.")

    md: List[str] = []
    md.append("# COM Trajectory Comparison Report")
    md.append("")
    md.append("本报告比较各数值试验中基线方法与 SDF 方法（以及 M2 单向/双向）的全时程质心轨迹。")
    md.append("差异指标均由同一步号、同一刚体 ID 对齐后计算。")
    md.append("")
    md.append("指标说明：")
    md.append("- `rmse_x/y/z`: 各坐标方向质心差的均方根误差。")
    md.append("- `rmse_l2`: 三维位置差欧氏范数的均方根。")
    md.append("- `max_l2`: 全时程最大三维质心差。")
    md.append("- `final_l2`: 末步三维质心差。")
    md.append("")

    analyze_mesh(mesh_rows, args.out_dir, md)
    analyze_scene(scene_rows, args.out_dir, md)
    analyze_scale(scale_rows, args.out_dir, md)
    analyze_m2(m2_rows, args.out_dir, md)

    md.append("## Interpretation")
    md.append("")
    md.append("- 若 `rmse_l2` 与 `final_l2` 在可接受数量级内且系统未失稳，可视为 SDF 轨迹与基线动力学演化一致。")
    md.append("- 若局部工况出现较大偏差，应结合该工况接触点数量、穿透深度与失稳率共同判断。")
    md.append("- 对 M2，单向与双向轨迹差异反映了方向偏置带来的动力学路径分歧。")
    md.append("")

    out_summary = args.out_summary
    if not os.path.isabs(out_summary):
        out_summary = os.path.join(args.out_dir, out_summary)
    with open(out_summary, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(md))
    print(f"Wrote summary: {out_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
