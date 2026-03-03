#!/usr/bin/env python3
"""Summarize W3 narrowphase ablation CSV directories into merged CSV/Markdown."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Dict, Iterable, List, Tuple


FILENAME_RE = re.compile(
    r"^stack3_seed_\d+_\d+_mode_(?P<cluster>on|off)_(?P<np>[a-z_]+)_k(?P<kmax>\d+)_c(?P<cell>[0-9.]+)\.csv$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--release-dir",
        default="build/chrono_sdf_contact_vs/Release",
        help="Chrono release output directory.",
    )
    parser.add_argument(
        "--output-dir",
        default="build/chrono_sdf_contact_vs/Release/results_w3_np_ablation_20seeds_v1",
        help="Directory where merged CSV/Markdown will be written.",
    )
    parser.add_argument(
        "--mode-dirs",
        nargs="+",
        default=[
            "default=results_w3_np_default_20seeds",
            "no_cache=results_w3_np_nocache_20seeds",
            "no_cull=results_w3_np_nocull_20seeds",
            "no_cache_no_cull=results_w3_np_nocache_nocull_20seeds",
        ],
        help="Mappings in the form mode=subdir.",
    )
    parser.add_argument(
        "--step-size",
        type=float,
        default=0.003,
        help="Only rows with this step_size are included.",
    )
    parser.add_argument(
        "--kmax-filter",
        default="64,128",
        help="Comma-separated kmax values to keep.",
    )
    parser.add_argument(
        "--cell-size-filter",
        default="0.05",
        help="Comma-separated cell_size values to keep (string exact match).",
    )
    return parser.parse_args()


def parse_mode_dirs(entries: Iterable[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for entry in entries:
        if "=" not in entry:
            raise ValueError(f"Invalid --mode-dirs entry: {entry!r}")
        key, value = entry.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key or not value:
            raise ValueError(f"Invalid --mode-dirs entry: {entry!r}")
        out[key] = value
    return out


def parse_token_set(text: str) -> set[str]:
    return {token.strip() for token in text.split(",") if token.strip()}


def fmt_percent(v: float) -> str:
    return f"{100.0 * v:.1f}%"


def main() -> int:
    args = parse_args()
    release_dir = Path(args.release_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mode_dirs = parse_mode_dirs(args.mode_dirs)
    keep_kmax = parse_token_set(args.kmax_filter)
    keep_cell = parse_token_set(args.cell_size_filter)

    # key: (declared_mode, np_mode, cluster_mode, kmax, cell, scene, contact_mode)
    metric_rows: Dict[Tuple[str, str, str, str, str, str, str], Dict[str, List[float]]] = defaultdict(
        lambda: {"exploded": [], "max_penetration": [], "wall_time_s": []}
    )

    for declared_mode, subdir in mode_dirs.items():
        directory = release_dir / subdir
        if not directory.exists():
            raise FileNotFoundError(f"Mode directory not found: {directory}")

        for path in sorted(directory.glob("stack3_seed_*.csv")):
            m = FILENAME_RE.match(path.name)
            if not m:
                continue
            cluster_mode = m.group("cluster")
            np_mode = m.group("np")
            kmax = m.group("kmax")
            cell = m.group("cell")
            if keep_kmax and kmax not in keep_kmax:
                continue
            if keep_cell and cell not in keep_cell:
                continue

            with path.open("r", encoding="utf-8", newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    step = float(row["step_size"])
                    if abs(step - args.step_size) > 1e-12:
                        continue
                    key = (
                        declared_mode,
                        np_mode,
                        cluster_mode,
                        kmax,
                        cell,
                        row["scene"],
                        row["mode"],
                    )
                    metric_rows[key]["exploded"].append(float(row["exploded"]))
                    metric_rows[key]["max_penetration"].append(float(row["max_penetration"]))
                    metric_rows[key]["wall_time_s"].append(float(row["wall_time_s"]))

    merged_csv_path = output_dir / "w3_np_scene_metrics.csv"
    with merged_csv_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "declared_mode",
                "narrowphase_mode",
                "clustering_mode",
                "kmax",
                "cell_size",
                "scene",
                "contact_mode",
                "num_rows",
                "exploded_rate",
                "max_penetration_mean",
                "wall_time_mean",
            ]
        )
        for key in sorted(metric_rows.keys()):
            values = metric_rows[key]
            n = len(values["exploded"])
            writer.writerow(
                list(key)
                + [
                    n,
                    mean(values["exploded"]) if n else "",
                    mean(values["max_penetration"]) if n else "",
                    mean(values["wall_time_s"]) if n else "",
                ]
            )

    # Build compact table: exploded rate (sdf_contact) by scene and narrowphase mode.
    scenes = ["box_plane_rest", "box_plane_impact", "nonconvex_t_drop", "stack3_boxes"]
    table_rows: List[Tuple[str, str, str, str, str, str]] = []
    for kmax in sorted(keep_kmax, key=lambda x: int(x)):
        for declared_mode in sorted(mode_dirs.keys()):
            row = [kmax, declared_mode]
            for scene in scenes:
                matched_keys = [
                    k
                    for k in metric_rows.keys()
                    if k[0] == declared_mode and k[3] == kmax and k[5] == scene and k[6] == "sdf_contact"
                ]
                if not matched_keys:
                    row.append("NA")
                    continue
                # Prefer cluster=on when available, then deterministic by cell size.
                key = sorted(
                    matched_keys,
                    key=lambda x: (x[2] != "on", float(x[4]), x[1]),
                )[0]
                exploded_rate = mean(metric_rows[key]["exploded"])
                row.append(fmt_percent(exploded_rate))
            table_rows.append(tuple(row))

    md_path = output_dir / "w3_np_summary.md"
    with md_path.open("w", encoding="utf-8") as f:
        f.write("# W3 Narrowphase Ablation Summary\n\n")
        f.write(f"- release_dir: `{release_dir}`\n")
        f.write(f"- step_size: `{args.step_size}`\n")
        f.write(f"- kmax_filter: `{','.join(sorted(keep_kmax, key=lambda x: int(x)))}`\n")
        f.write(f"- cell_size_filter: `{','.join(sorted(keep_cell))}`\n")
        f.write("\n## SDF exploded_rate by scene\n\n")
        f.write("| kmax | narrowphase_mode | box_plane_rest | box_plane_impact | nonconvex_t_drop | stack3_boxes |\n")
        f.write("| ---: | --- | ---: | ---: | ---: | ---: |\n")
        for row in table_rows:
            f.write(
                f"| {row[0]} | `{row[1]}` | {row[2]} | {row[3]} | {row[4]} | {row[5]} |\n"
            )
        f.write("\n## Output Files\n\n")
        f.write(f"1. `{merged_csv_path}`\n")
        f.write(f"2. `{md_path}`\n")

    print(f"Wrote {merged_csv_path}")
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
