#!/usr/bin/env python3
import argparse
import csv
import math
import random
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


def to_int(row, key, default=0):
    try:
        return int(float(row.get(key, default)))
    except Exception:
        return default


def to_float_or_none(row, key):
    value = row.get(key, "")
    if value is None or value == "":
        return None
    try:
        return float(value)
    except Exception:
        return None


def mean_ci95(values):
    if not values:
        return None, None
    n = len(values)
    mean = sum(values) / float(n)
    if n < 2:
        return mean, 0.0
    var = sum((x - mean) ** 2 for x in values) / float(n - 1)
    std = math.sqrt(max(var, 0.0))
    ci95 = 1.96 * std / math.sqrt(float(n))
    return mean, ci95


def mean(values):
    if not values:
        return None
    return sum(values) / float(len(values))


def permutation_pvalue(values_a, values_b, samples=2000, seed=12345):
    a = [float(v) for v in values_a if v is not None and math.isfinite(float(v))]
    b = [float(v) for v in values_b if v is not None and math.isfinite(float(v))]
    if len(a) < 2 or len(b) < 2:
        return None
    if samples < 1:
        return None

    obs = abs(mean(a) - mean(b))
    combined = a + b
    n_a = len(a)
    n_total = len(combined)
    if n_total <= 2:
        return None

    rng = random.Random(seed + n_a * 131 + len(b) * 17 + n_total * 19)
    work = combined[:]
    hits = 0
    for _ in range(samples):
        rng.shuffle(work)
        ma = sum(work[:n_a]) / float(n_a)
        mb = sum(work[n_a:]) / float(n_total - n_a)
        if abs(ma - mb) >= obs - 1e-15:
            hits += 1
    return (hits + 1) / float(samples + 1)


def aggregate(rows, group_keys, metric_keys):
    groups = defaultdict(list)
    for r in rows:
        key = tuple(r.get(k, "") for k in group_keys)
        groups[key].append(r)

    out = []
    for key, grows in groups.items():
        row = {k: v for k, v in zip(group_keys, key)}
        row["n"] = len(grows)
        for mk in metric_keys:
            vals = []
            for r in grows:
                fv = to_float_or_none(r, mk)
                if fv is not None:
                    vals.append(fv)
            mean, ci95 = mean_ci95(vals)
            row[f"{mk}_mean"] = mean
            row[f"{mk}_ci95"] = ci95
        out.append(row)
    return out


def fmt(x, digits=6):
    if x is None:
        return "-"
    if isinstance(x, float):
        if math.isinf(x) or math.isnan(x):
            return "-"
        return f"{x:.{digits}g}"
    return str(x)


def fmt_mean_ci(mean, ci, digits=6):
    if mean is None:
        return "-"
    if ci is None or ci <= 0.0:
        return fmt(mean, digits)
    return f"{fmt(mean, digits)} +/- {fmt(ci, digits)}"


def loglog_fit(x_values, y_values):
    pts = []
    for x, y in zip(x_values, y_values):
        if x is None or y is None:
            continue
        if x <= 0.0 or y <= 0.0:
            continue
        if (not math.isfinite(x)) or (not math.isfinite(y)):
            continue
        pts.append((math.log(x), math.log(y)))

    if len(pts) < 2:
        return None

    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    n = len(pts)
    x_mean = sum(xs) / float(n)
    y_mean = sum(ys) / float(n)
    sxx = sum((x - x_mean) ** 2 for x in xs)
    if sxx <= 1e-16:
        return None
    sxy = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = y_mean - slope * x_mean

    y_pred = [intercept + slope * x for x in xs]
    sse = sum((yt - yp) ** 2 for yt, yp in zip(ys, y_pred))
    sst = sum((yt - y_mean) ** 2 for yt in ys)
    if sst <= 1e-16:
        r2 = 1.0 if sse <= 1e-16 else 0.0
    else:
        r2 = 1.0 - sse / sst

    return {
        "n": n,
        "slope": slope,
        "intercept": intercept,
        "prefactor": math.exp(intercept),
        "r2": r2,
    }


def mesh_section(mesh_rows):
    lines = []
    lines.append("## M1 Mesh Baseline Compare")
    lines.append("")

    agg = aggregate(
        mesh_rows,
        group_keys=["mode"],
        metric_keys=[
            "wall_time_s",
            "avg_contacts",
            "max_contacts",
            "max_penetration",
            "final_mesh_height",
            "max_relative_energy_drift",
        ],
    )
    agg.sort(key=lambda r: r.get("mode", ""))

    lines.append(
        "| mode | n | wall_time_s(mean+/-ci95) | avg_contacts(mean+/-ci95) | max_contacts(mean+/-ci95) | "
        "max_penetration(mean+/-ci95) | max_rel_energy_drift(mean+/-ci95) | final_mesh_height(mean+/-ci95) |"
    )
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in agg:
        lines.append(
            "| {mode} | {n} | {wall} | {avgc} | {maxc} | {pen} | {drift} | {h} |".format(
                mode=r.get("mode", "-"),
                n=r.get("n", 0),
                wall=fmt_mean_ci(r.get("wall_time_s_mean"), r.get("wall_time_s_ci95")),
                avgc=fmt_mean_ci(r.get("avg_contacts_mean"), r.get("avg_contacts_ci95")),
                maxc=fmt_mean_ci(r.get("max_contacts_mean"), r.get("max_contacts_ci95")),
                pen=fmt_mean_ci(r.get("max_penetration_mean"), r.get("max_penetration_ci95")),
                drift=fmt_mean_ci(
                    r.get("max_relative_energy_drift_mean"),
                    r.get("max_relative_energy_drift_ci95"),
                ),
                h=fmt_mean_ci(r.get("final_mesh_height_mean"), r.get("final_mesh_height_ci95")),
            )
        )

    by_mode = {r.get("mode", ""): r for r in agg}
    b = by_mode.get("chrono_baseline")
    s = by_mode.get("sdf_contact")
    if b and s:
        b_wall = b.get("wall_time_s_mean")
        s_wall = s.get("wall_time_s_mean")
        b_pen = b.get("max_penetration_mean")
        s_pen = s.get("max_penetration_mean")
        b_avg = b.get("avg_contacts_mean")
        s_avg = s.get("avg_contacts_mean")
        b_drift = b.get("max_relative_energy_drift_mean")
        s_drift = s.get("max_relative_energy_drift_mean")
        speedup = (b_wall / s_wall) if (b_wall is not None and s_wall and s_wall > 0) else None

        lines.append("")
        lines.append("Key metrics:")
        lines.append(f"- `sdf_contact` time ratio vs baseline: `{fmt(speedup)}` (baseline/sdf)")
        lines.append(
            f"- penetration delta mean (`sdf - baseline`): `{fmt((s_pen - b_pen) if (s_pen is not None and b_pen is not None) else None)}`"
        )
        lines.append(
            f"- average contact count delta mean (`sdf - baseline`): `{fmt((s_avg - b_avg) if (s_avg is not None and b_avg is not None) else None)}`"
        )
        lines.append(
            f"- max relative energy drift delta mean (`sdf - baseline`): `{fmt((s_drift - b_drift) if (s_drift is not None and b_drift is not None) else None)}`"
        )
    lines.append("")
    return lines


def ablation_coverage_section(mesh_rows, stability_rows):
    lines = []
    lines.append("## P0 Ablation Coverage")
    lines.append("")

    has_baseline = any(r.get("mode", "") == "chrono_baseline" for r in mesh_rows)
    has_one_sided = any(("a_to_b_only" in r.get("mode", "")) or ("b_to_a_only" in r.get("mode", "")) for r in stability_rows)
    has_bidirectional = any(r.get("mode", "").startswith("bidirectional") for r in stability_rows)

    has_cluster_on = False
    has_cluster_off = False
    has_octree = False
    for r in stability_rows:
        mode = r.get("mode", "")
        cflag = to_int(r, "clustering_enabled", -1)
        oflag = to_int(r, "octree_enabled", -1)
        cluster_on = (cflag == 1) or ("no_cluster" not in mode)
        cluster_off = (cflag == 0) or ("no_cluster" in mode)
        octree_on = (oflag == 1) or ("octree" in mode)
        has_cluster_on = has_cluster_on or cluster_on
        has_cluster_off = has_cluster_off or cluster_off
        has_octree = has_octree or octree_on

    lines.append("| item | covered | evidence |")
    lines.append("| --- | --- | --- |")
    lines.append(f"| baseline | {'yes' if has_baseline else 'no'} | mode=`chrono_baseline` in mesh compare |")
    lines.append(f"| one-sided | {'yes' if has_one_sided else 'no'} | mode contains `a_to_b_only`/`b_to_a_only` |")
    lines.append(f"| bidirectional | {'yes' if has_bidirectional else 'no'} | mode starts with `bidirectional` |")
    lines.append(
        f"| clustering | {'yes' if (has_cluster_on and has_cluster_off) else 'no'} | both clustered and no-cluster modes present |"
    )
    lines.append(f"| octree-enabled | {'yes' if has_octree else 'no'} | mode/flag indicates octree-enabled runs |")
    lines.append("")
    return lines


def stability_section(stability_rows, sig_samples=2000, sig_seed=12345):
    lines = []
    lines.append("## M2 SDF-SDF Stability Compare")
    lines.append("")

    agg = aggregate(
        stability_rows,
        group_keys=["step_size", "mode", "clustering_enabled", "octree_enabled"],
        metric_keys=[
            "wall_time_s",
            "avg_contacts",
            "max_contacts",
            "max_penetration",
            "max_relative_energy_drift",
            "exploded",
        ],
    )
    agg.sort(key=lambda r: (to_float(r, "step_size"), r.get("mode", "")))

    lines.append(
        "| step_size | mode | clustering | octree | n | wall_time_s(mean+/-ci95) | avg_contacts(mean+/-ci95) | "
        "max_contacts(mean+/-ci95) | max_penetration(mean+/-ci95) | max_rel_energy_drift(mean+/-ci95) | instability_rate |"
    )
    lines.append("| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in agg:
        instab = r.get("exploded_mean")
        instab_pct = f"{fmt(100.0 * instab, 4)}%" if instab is not None else "-"
        lines.append(
            "| {dt} | {mode} | {cluster} | {octree} | {n} | {wall} | {avgc} | {maxc} | {pen} | {drift} | {instab} |".format(
                dt=fmt(to_float(r, "step_size")),
                mode=r.get("mode", "-"),
                cluster=fmt(to_int(r, "clustering_enabled", 1)),
                octree=fmt(to_int(r, "octree_enabled", 0)),
                n=r.get("n", 0),
                wall=fmt_mean_ci(r.get("wall_time_s_mean"), r.get("wall_time_s_ci95")),
                avgc=fmt_mean_ci(r.get("avg_contacts_mean"), r.get("avg_contacts_ci95")),
                maxc=fmt_mean_ci(r.get("max_contacts_mean"), r.get("max_contacts_ci95")),
                pen=fmt_mean_ci(r.get("max_penetration_mean"), r.get("max_penetration_ci95")),
                drift=fmt_mean_ci(
                    r.get("max_relative_energy_drift_mean"),
                    r.get("max_relative_energy_drift_ci95"),
                ),
                instab=instab_pct,
            )
        )

    grouped = defaultdict(dict)
    for r in agg:
        dt = to_float(r, "step_size")
        cluster = to_int(r, "clustering_enabled", 1)
        octree = to_int(r, "octree_enabled", 0)
        grouped[(dt, cluster, octree)][r.get("mode", "")] = r

    lines.append("")
    lines.append("Directional bias indicators (using mean values):")
    lines.append("")
    lines.append(
        "| step_size | clustering | octree | penetration_gap(|A->B - B->A|) | avg_contacts_gap(|A->B - B->A|) | "
        "energy_drift_gap(|A->B - B->A|) | bidirectional_instability_rate |"
    )
    lines.append("| ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for key in sorted(grouped.keys()):
        dt, cluster, octree = key
        modes = grouped[key]
        a = modes.get("a_to_b_only")
        b = modes.get("b_to_a_only")
        bi = modes.get("bidirectional")
        if not (a and b and bi):
            a = modes.get("a_to_b_only_no_cluster")
            b = modes.get("b_to_a_only_no_cluster")
            bi = modes.get("bidirectional_no_cluster")
        if a and b:
            pa = a.get("max_penetration_mean")
            pb = b.get("max_penetration_mean")
            ca = a.get("avg_contacts_mean")
            cb = b.get("avg_contacts_mean")
            da = a.get("max_relative_energy_drift_mean")
            db = b.get("max_relative_energy_drift_mean")
            pen_gap = abs(pa - pb) if pa is not None and pb is not None else None
            cnt_gap = abs(ca - cb) if ca is not None and cb is not None else None
            drift_gap = abs(da - db) if da is not None and db is not None else None
        else:
            pen_gap = None
            cnt_gap = None
            drift_gap = None

        bi_instab = bi.get("exploded_mean") if bi else None
        bi_instab_pct = (f"{fmt(100.0 * bi_instab, 4)}%" if bi_instab is not None else "-")
        lines.append(
            f"| {fmt(dt)} | {fmt(cluster)} | {fmt(octree)} | {fmt(pen_gap)} | {fmt(cnt_gap)} | {fmt(drift_gap)} | {bi_instab_pct} |"
        )

    exploded_count = defaultdict(int)
    total_count = defaultdict(int)
    for r in stability_rows:
        mode = r.get("mode", "-")
        total_count[mode] += 1
        exploded_count[mode] += to_int(r, "exploded")

    lines.append("")
    lines.append("Stability summary:")
    for mode in sorted(total_count.keys()):
        total = total_count[mode]
        exploded = exploded_count[mode]
        rate = (100.0 * exploded / float(total)) if total > 0 else 0.0
        lines.append(f"- `{mode}` exploded `{exploded}/{total}` cases (`{fmt(rate, 4)}%`)")

    def rows_for(dt, primary_mode, fallback_mode):
        prim = [
            r
            for r in stability_rows
            if abs(to_float(r, "step_size") - dt) < 1e-12 and r.get("mode", "") == primary_mode
        ]
        if prim:
            return prim
        return [
            r
            for r in stability_rows
            if abs(to_float(r, "step_size") - dt) < 1e-12 and r.get("mode", "") == fallback_mode
        ]

    dt_values = sorted({to_float(r, "step_size") for r in stability_rows})
    lines.append("")
    lines.append(f"Significance checks (permutation test, two-sided, samples={sig_samples}):")
    lines.append("")
    lines.append(
        "| step_size | compare | metric | bidirectional_mean | one_sided_mean | delta(one-bid) | p_value |"
    )
    lines.append("| ---: | --- | --- | ---: | ---: | ---: | ---: |")

    for dt in dt_values:
        bid_rows = rows_for(dt, "bidirectional", "bidirectional_no_cluster")
        if not bid_rows:
            continue
        bid_expl = [to_float(r, "exploded") for r in bid_rows]
        bid_pen = [to_float(r, "max_penetration") for r in bid_rows]
        bid_expl_mean = mean(bid_expl)
        bid_pen_mean = mean(bid_pen)

        for primary, fallback in [("a_to_b_only", "a_to_b_only_no_cluster"), ("b_to_a_only", "b_to_a_only_no_cluster")]:
            one_rows = rows_for(dt, primary, fallback)
            if not one_rows:
                continue

            one_expl = [to_float(r, "exploded") for r in one_rows]
            one_pen = [to_float(r, "max_penetration") for r in one_rows]
            one_expl_mean = mean(one_expl)
            one_pen_mean = mean(one_pen)

            p_expl = permutation_pvalue(one_expl, bid_expl, samples=sig_samples, seed=sig_seed + int(dt * 1e6) + 11)
            p_pen = permutation_pvalue(one_pen, bid_pen, samples=sig_samples, seed=sig_seed + int(dt * 1e6) + 23)

            lines.append(
                "| {dt} | {cmp} | exploded | {bm} | {om} | {delta} | {p} |".format(
                    dt=fmt(dt),
                    cmp=f"bidirectional vs {primary}",
                    bm=fmt(bid_expl_mean),
                    om=fmt(one_expl_mean),
                    delta=fmt((one_expl_mean - bid_expl_mean) if (one_expl_mean is not None and bid_expl_mean is not None) else None),
                    p=fmt(p_expl),
                )
            )
            lines.append(
                "| {dt} | {cmp} | max_penetration | {bm} | {om} | {delta} | {p} |".format(
                    dt=fmt(dt),
                    cmp=f"bidirectional vs {primary}",
                    bm=fmt(bid_pen_mean),
                    om=fmt(one_pen_mean),
                    delta=fmt((one_pen_mean - bid_pen_mean) if (one_pen_mean is not None and bid_pen_mean is not None) else None),
                    p=fmt(p_pen),
                )
            )
    lines.append("")
    return lines


def p2_scene_section(scene_rows, sig_samples=2000, sig_seed=12345):
    lines = []
    lines.append("## P2 Diverse-Scene Benchmark Matrix")
    lines.append("")

    agg = aggregate(
        scene_rows,
        group_keys=["scene", "step_size", "mode"],
        metric_keys=[
            "wall_time_s",
            "avg_contacts",
            "max_contacts",
            "max_penetration",
            "max_relative_energy_drift",
            "exploded",
        ],
    )
    agg.sort(key=lambda r: (r.get("scene", ""), to_float(r, "step_size"), r.get("mode", "")))

    lines.append(
        "| scene | step_size | mode | n | wall_time_s(mean+/-ci95) | avg_contacts(mean+/-ci95) | "
        "max_contacts(mean+/-ci95) | max_penetration(mean+/-ci95) | max_rel_energy_drift(mean+/-ci95) | instability_rate |"
    )
    lines.append("| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in agg:
        instab = r.get("exploded_mean")
        instab_pct = f"{fmt(100.0 * instab, 4)}%" if instab is not None else "-"
        lines.append(
            "| {scene} | {dt} | {mode} | {n} | {wall} | {avgc} | {maxc} | {pen} | {drift} | {instab} |".format(
                scene=r.get("scene", "-"),
                dt=fmt(to_float(r, "step_size")),
                mode=r.get("mode", "-"),
                n=r.get("n", 0),
                wall=fmt_mean_ci(r.get("wall_time_s_mean"), r.get("wall_time_s_ci95")),
                avgc=fmt_mean_ci(r.get("avg_contacts_mean"), r.get("avg_contacts_ci95")),
                maxc=fmt_mean_ci(r.get("max_contacts_mean"), r.get("max_contacts_ci95")),
                pen=fmt_mean_ci(r.get("max_penetration_mean"), r.get("max_penetration_ci95")),
                drift=fmt_mean_ci(
                    r.get("max_relative_energy_drift_mean"),
                    r.get("max_relative_energy_drift_ci95"),
                ),
                instab=instab_pct,
            )
        )

    scene_names = sorted(set(r.get("scene", "") for r in agg))
    lines.append("")
    lines.append(f"- covered scenes: `{len(scene_names)}` -> {', '.join(scene_names)}")

    lines.append("")
    lines.append(f"Significance checks (permutation test, two-sided, samples={sig_samples}):")
    lines.append("")
    lines.append("| scope | metric | chrono_baseline_mean | sdf_contact_mean | delta(sdf-baseline) | p_value |")
    lines.append("| --- | --- | ---: | ---: | ---: | ---: |")

    def metric_vals(rows, mode, metric):
        return [to_float(r, metric) for r in rows if r.get("mode", "") == mode]

    for metric in ["wall_time_s", "max_penetration", "exploded"]:
        b_vals = metric_vals(scene_rows, "chrono_baseline", metric)
        s_vals = metric_vals(scene_rows, "sdf_contact", metric)
        b_mean = mean(b_vals)
        s_mean = mean(s_vals)
        p = permutation_pvalue(b_vals, s_vals, samples=sig_samples, seed=sig_seed + 101 + len(metric))
        lines.append(
            "| all-scenes | {metric} | {bm} | {sm} | {delta} | {p} |".format(
                metric=metric,
                bm=fmt(b_mean),
                sm=fmt(s_mean),
                delta=fmt((s_mean - b_mean) if (s_mean is not None and b_mean is not None) else None),
                p=fmt(p),
            )
        )

    for scene in scene_names:
        scene_subset = [r for r in scene_rows if r.get("scene", "") == scene]
        b_vals = metric_vals(scene_subset, "chrono_baseline", "wall_time_s")
        s_vals = metric_vals(scene_subset, "sdf_contact", "wall_time_s")
        b_mean = mean(b_vals)
        s_mean = mean(s_vals)
        p = permutation_pvalue(b_vals, s_vals, samples=sig_samples, seed=sig_seed + 211 + len(scene))
        lines.append(
            "| {scene} | wall_time_s | {bm} | {sm} | {delta} | {p} |".format(
                scene=scene,
                bm=fmt(b_mean),
                sm=fmt(s_mean),
                delta=fmt((s_mean - b_mean) if (s_mean is not None and b_mean is not None) else None),
                p=fmt(p),
            )
        )

    lines.append("")
    lines.append(f"Exploded significance by scene-step (permutation test, two-sided, samples={sig_samples}):")
    lines.append("")
    lines.append(
        "| scene | step_size | n_baseline | n_sdf | chrono_baseline_exploded_rate | "
        "sdf_contact_exploded_rate | delta(sdf-baseline) | p_value |"
    )
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")

    scene_dt_pairs = sorted(
        {(r.get("scene", ""), to_float(r, "step_size")) for r in scene_rows},
        key=lambda x: (x[0], x[1]),
    )
    for scene, dt in scene_dt_pairs:
        subset = [
            r
            for r in scene_rows
            if r.get("scene", "") == scene and abs(to_float(r, "step_size") - dt) < 1e-12
        ]
        b_vals = metric_vals(subset, "chrono_baseline", "exploded")
        s_vals = metric_vals(subset, "sdf_contact", "exploded")
        b_mean = mean(b_vals)
        s_mean = mean(s_vals)
        scene_code = sum(ord(ch) for ch in scene)
        dt_code = int(round(dt * 1e6))
        p = permutation_pvalue(
            b_vals,
            s_vals,
            samples=sig_samples,
            seed=sig_seed + 503 + scene_code * 17 + dt_code,
        )
        b_rate = (f"{fmt(100.0 * b_mean, 4)}%" if b_mean is not None else "-")
        s_rate = (f"{fmt(100.0 * s_mean, 4)}%" if s_mean is not None else "-")
        delta = (f"{fmt(100.0 * (s_mean - b_mean), 4)}%" if (s_mean is not None and b_mean is not None) else "-")
        lines.append(
            "| {scene} | {dt} | {nb} | {ns} | {br} | {sr} | {delta} | {p} |".format(
                scene=scene,
                dt=fmt(dt),
                nb=len(b_vals),
                ns=len(s_vals),
                br=b_rate,
                sr=s_rate,
                delta=delta,
                p=fmt(p),
            )
        )
    lines.append("")
    return lines


def p2_scale_section(scale_rows, sig_samples=2000, sig_seed=12345):
    lines = []
    lines.append("## P2 Scale-Up Curves")
    lines.append("")

    agg = aggregate(
        scale_rows,
        group_keys=["mode", "body_count"],
        metric_keys=[
            "wall_time_s",
            "avg_contacts",
            "max_contacts",
            "max_penetration",
            "max_relative_energy_drift",
            "exploded",
        ],
    )
    agg.sort(key=lambda r: (r.get("mode", ""), to_float(r, "body_count")))

    lines.append(
        "| mode | body_count | n | wall_time_s(mean+/-ci95) | avg_contacts(mean+/-ci95) | "
        "max_contacts(mean+/-ci95) | max_penetration(mean+/-ci95) | max_rel_energy_drift(mean+/-ci95) | instability_rate |"
    )
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in agg:
        instab = r.get("exploded_mean")
        instab_pct = f"{fmt(100.0 * instab, 4)}%" if instab is not None else "-"
        lines.append(
            "| {mode} | {count} | {n} | {wall} | {avgc} | {maxc} | {pen} | {drift} | {instab} |".format(
                mode=r.get("mode", "-"),
                count=fmt(to_float(r, "body_count")),
                n=r.get("n", 0),
                wall=fmt_mean_ci(r.get("wall_time_s_mean"), r.get("wall_time_s_ci95")),
                avgc=fmt_mean_ci(r.get("avg_contacts_mean"), r.get("avg_contacts_ci95")),
                maxc=fmt_mean_ci(r.get("max_contacts_mean"), r.get("max_contacts_ci95")),
                pen=fmt_mean_ci(r.get("max_penetration_mean"), r.get("max_penetration_ci95")),
                drift=fmt_mean_ci(
                    r.get("max_relative_energy_drift_mean"),
                    r.get("max_relative_energy_drift_ci95"),
                ),
                instab=instab_pct,
            )
        )

    lines.append("")
    lines.append("Scale-up log-log fits (`wall_time_s ~ C * body_count^p`):")
    lines.append("")
    lines.append("| mode | n | slope p | prefactor C | R2 |")
    lines.append("| --- | ---: | ---: | ---: | ---: |")

    modes = sorted(set(r.get("mode", "") for r in agg))
    for mode in modes:
        mode_rows = [r for r in agg if r.get("mode", "") == mode]
        x = [to_float(r, "body_count", None) for r in mode_rows]
        y = [r.get("wall_time_s_mean") for r in mode_rows]
        fit = loglog_fit(x, y)
        if fit:
            lines.append(
                "| {mode} | {n} | {p} | {c} | {r2} |".format(
                    mode=mode,
                    n=fit["n"],
                    p=fmt(fit["slope"]),
                    c=fmt(fit["prefactor"]),
                    r2=fmt(fit["r2"]),
                )
            )
        else:
            lines.append(f"| {mode} | - | - | - | - |")

    lines.append("")
    lines.append(f"Significance checks (permutation test, two-sided, samples={sig_samples}):")
    lines.append("")
    lines.append("| scope | metric | chrono_baseline_mean | sdf_contact_mean | delta(sdf-baseline) | p_value |")
    lines.append("| --- | --- | ---: | ---: | ---: | ---: |")

    def metric_vals(rows, mode, metric):
        return [to_float(r, metric) for r in rows if r.get("mode", "") == mode]

    for metric in ["wall_time_s", "max_penetration", "exploded"]:
        b_vals = metric_vals(scale_rows, "chrono_baseline", metric)
        s_vals = metric_vals(scale_rows, "sdf_contact", metric)
        b_mean = mean(b_vals)
        s_mean = mean(s_vals)
        p = permutation_pvalue(b_vals, s_vals, samples=sig_samples, seed=sig_seed + 307 + len(metric))
        lines.append(
            "| all-body-counts | {metric} | {bm} | {sm} | {delta} | {p} |".format(
                metric=metric,
                bm=fmt(b_mean),
                sm=fmt(s_mean),
                delta=fmt((s_mean - b_mean) if (s_mean is not None and b_mean is not None) else None),
                p=fmt(p),
            )
        )

    counts = sorted({to_int(r, "body_count", 0) for r in scale_rows if to_int(r, "body_count", 0) > 0})
    for count in counts:
        subset = [r for r in scale_rows if to_int(r, "body_count", 0) == count]
        b_vals = metric_vals(subset, "chrono_baseline", "wall_time_s")
        s_vals = metric_vals(subset, "sdf_contact", "wall_time_s")
        b_mean = mean(b_vals)
        s_mean = mean(s_vals)
        p = permutation_pvalue(b_vals, s_vals, samples=sig_samples, seed=sig_seed + 401 + count)
        lines.append(
            "| body_count={count} | wall_time_s | {bm} | {sm} | {delta} | {p} |".format(
                count=count,
                bm=fmt(b_mean),
                sm=fmt(s_mean),
                delta=fmt((s_mean - b_mean) if (s_mean is not None and b_mean is not None) else None),
                p=fmt(p),
            )
        )

    lines.append("")
    return lines


def m3_section(m3_rows):
    lines = []
    lines.append("## M3 Octree Narrow-Band And Error Control")
    lines.append("")
    agg = aggregate(
        m3_rows,
        group_keys=["eps_force"],
        metric_keys=[
            "node_count",
            "leaf_count",
            "min_leaf_size",
            "est_delta_f",
            "query_time_s",
            "rmse_phi",
            "max_abs_phi_err",
            "rmse_force_err",
            "rmse_force_err_linear",
            "force_map_mae",
            "force_map_rel_mae",
            "force_map_r2_fixed",
            "force_map_slope",
            "force_map_bound_violation_rate",
        ],
    )
    rows = sorted(agg, key=lambda r: to_float(r, "eps_force"), reverse=True)
    lines.append(
        "| eps_force | n | node_count(mean+/-ci95) | leaf_count(mean+/-ci95) | min_leaf_size(mean+/-ci95) | "
        "est_delta_f(mean+/-ci95) | query_time_s(mean+/-ci95) | rmse_phi(mean+/-ci95) | max_abs_phi_err(mean+/-ci95) |"
    )
    lines.append("| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in rows:
        lines.append(
            "| {eps} | {n} | {nodes} | {leaves} | {hmin} | {df} | {qt} | {rmse} | {maxe} |".format(
                eps=fmt(to_float(r, "eps_force")),
                n=r.get("n", 0),
                nodes=fmt_mean_ci(r.get("node_count_mean"), r.get("node_count_ci95")),
                leaves=fmt_mean_ci(r.get("leaf_count_mean"), r.get("leaf_count_ci95")),
                hmin=fmt_mean_ci(r.get("min_leaf_size_mean"), r.get("min_leaf_size_ci95")),
                df=fmt_mean_ci(r.get("est_delta_f_mean"), r.get("est_delta_f_ci95")),
                qt=fmt_mean_ci(r.get("query_time_s_mean"), r.get("query_time_s_ci95")),
                rmse=fmt_mean_ci(r.get("rmse_phi_mean"), r.get("rmse_phi_ci95")),
                maxe=fmt_mean_ci(r.get("max_abs_phi_err_mean"), r.get("max_abs_phi_err_ci95")),
            )
        )

    if rows:
        best_rmse = min(rows, key=lambda r: r.get("rmse_phi_mean", float("inf")))
        best_time = min(rows, key=lambda r: r.get("query_time_s_mean", float("inf")))
        lines.append("")
        lines.append("Key metrics:")
        lines.append(
            "- minimum `rmse_phi` mean: eps_force=`{eps}`, value=`{val}`".format(
                eps=fmt(to_float(best_rmse, "eps_force")), val=fmt(best_rmse.get("rmse_phi_mean"))
            )
        )
        lines.append(
            "- minimum `query_time_s` mean: eps_force=`{eps}`, value=`{val}`".format(
                eps=fmt(to_float(best_time, "eps_force")), val=fmt(best_time.get("query_time_s_mean"))
            )
        )

        lines.append("")
        lines.append("P1.1 `log(Delta d)-log(h)` fit:")

        h_vals = [r.get("min_leaf_size_mean") for r in rows]
        rmse_vals = [r.get("rmse_phi_mean") for r in rows]
        max_vals = [r.get("max_abs_phi_err_mean") for r in rows]
        fit_rmse = loglog_fit(h_vals, rmse_vals)
        fit_max = loglog_fit(h_vals, max_vals)

        lines.append("| target Delta d | n | slope p | prefactor C | R2 | model |")
        lines.append("| --- | ---: | ---: | ---: | ---: | --- |")
        if fit_rmse:
            lines.append(
                "| rmse_phi | {n} | {p} | {c} | {r2} | Delta d ~= C * h^p |".format(
                    n=fit_rmse["n"],
                    p=fmt(fit_rmse["slope"]),
                    c=fmt(fit_rmse["prefactor"]),
                    r2=fmt(fit_rmse["r2"]),
                )
            )
        else:
            lines.append("| rmse_phi | - | - | - | - | insufficient points |")

        if fit_max:
            lines.append(
                "| max_abs_phi_err | {n} | {p} | {c} | {r2} | Delta d ~= C * h^p |".format(
                    n=fit_max["n"],
                    p=fmt(fit_max["slope"]),
                    c=fmt(fit_max["prefactor"]),
                    r2=fmt(fit_max["r2"]),
                )
            )
        else:
            lines.append("| max_abs_phi_err | - | - | - | - | insufficient points |")

        lines.append("")
        lines.append("P1.2 `Delta F = k_n * Delta d` mapping quality:")
        lines.append(
            "| eps_force | rmse_force_err(mean+/-ci95) | rmse_force_err_linear(mean+/-ci95) | "
            "rmse_ratio(actual/linear) | map_r2_fixed(mean+/-ci95) | map_slope(mean+/-ci95) | "
            "map_rel_mae(mean+/-ci95) | bound_violation_rate(mean+/-ci95) |"
        )
        lines.append("| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")

        for r in rows:
            rmse_actual = r.get("rmse_force_err_mean")
            rmse_linear = r.get("rmse_force_err_linear_mean")
            rmse_ratio = None
            if rmse_actual is not None and rmse_linear is not None and rmse_linear > 1e-16:
                rmse_ratio = rmse_actual / rmse_linear
            lines.append(
                "| {eps} | {fa} | {fl} | {ratio} | {r2} | {slope} | {relmae} | {viol} |".format(
                    eps=fmt(to_float(r, "eps_force")),
                    fa=fmt_mean_ci(r.get("rmse_force_err_mean"), r.get("rmse_force_err_ci95")),
                    fl=fmt_mean_ci(r.get("rmse_force_err_linear_mean"), r.get("rmse_force_err_linear_ci95")),
                    ratio=fmt(rmse_ratio),
                    r2=fmt_mean_ci(r.get("force_map_r2_fixed_mean"), r.get("force_map_r2_fixed_ci95")),
                    slope=fmt_mean_ci(r.get("force_map_slope_mean"), r.get("force_map_slope_ci95")),
                    relmae=fmt_mean_ci(r.get("force_map_rel_mae_mean"), r.get("force_map_rel_mae_ci95")),
                    viol=fmt_mean_ci(
                        r.get("force_map_bound_violation_rate_mean"),
                        r.get("force_map_bound_violation_rate_ci95"),
                    ),
                )
            )

        avg_rel_mae = mean_ci95([r.get("force_map_rel_mae_mean") for r in rows if r.get("force_map_rel_mae_mean") is not None])
        avg_r2 = mean_ci95([r.get("force_map_r2_fixed_mean") for r in rows if r.get("force_map_r2_fixed_mean") is not None])
        lines.append("")
        lines.append(
            "- aggregate mapping quality: "
            f"`map_r2_fixed ~= {fmt(avg_r2[0])}`; `map_rel_mae ~= {fmt(avg_rel_mae[0])}`"
        )
    lines.append("")
    return lines


def main():
    parser = argparse.ArgumentParser(description="Generate markdown summary from chrono_sdf_contact experiment CSV outputs.")
    parser.add_argument("--mesh-csv", required=True, help="Path to sdf_mesh_baseline_compare.csv")
    parser.add_argument("--stability-csv", required=True, help="Path to sdf_sdf_stability_compare.csv")
    parser.add_argument("--p2-scene-csv", default="", help="Optional path to sdf_p2_scene_matrix.csv")
    parser.add_argument("--p2-scale-csv", default="", help="Optional path to sdf_p2_scaleup.csv")
    parser.add_argument("--m3-csv", default="", help="Optional path to sdf_octree_pareto.csv")
    parser.add_argument("--out", default="results_summary.md", help="Output markdown path")
    parser.add_argument("--title", default="Chrono SDF Contact Experiment Summary", help="Summary title")
    parser.add_argument("--sig-samples", type=int, default=2000, help="Permutation-test sample count for p-value estimation")
    parser.add_argument("--sig-seed", type=int, default=12345, help="Random seed for permutation tests")
    args = parser.parse_args()

    mesh_rows = read_csv(args.mesh_csv)
    stability_rows = read_csv(args.stability_csv)
    p2_scene_rows = read_csv(args.p2_scene_csv) if args.p2_scene_csv else []
    p2_scale_rows = read_csv(args.p2_scale_csv) if args.p2_scale_csv else []
    m3_rows = read_csv(args.m3_csv) if args.m3_csv else []

    lines = [f"# {args.title}", ""]
    lines.extend(ablation_coverage_section(mesh_rows, stability_rows))
    lines.extend(mesh_section(mesh_rows))
    lines.extend(stability_section(stability_rows, sig_samples=args.sig_samples, sig_seed=args.sig_seed))
    if p2_scene_rows:
        lines.extend(p2_scene_section(p2_scene_rows, sig_samples=args.sig_samples, sig_seed=args.sig_seed))
    if p2_scale_rows:
        lines.extend(p2_scale_section(p2_scale_rows, sig_samples=args.sig_samples, sig_seed=args.sig_seed))
    if m3_rows:
        lines.extend(m3_section(m3_rows))

    with open(args.out, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines).rstrip() + "\n")

    print(f"Generated: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
