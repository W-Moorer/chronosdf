#!/usr/bin/env python3
import argparse
import csv
import math
from collections import defaultdict


def read_rows(path):
    rows = []
    with open(path, "r", newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def to_float(row, key, default=None):
    try:
        value = row.get(key, default)
        if value is None or value == "":
            return default
        return float(value)
    except Exception:
        return default


def mean(values):
    if not values:
        return None
    return sum(values) / float(len(values))


def aggregate_by_eps(rows):
    grouped = defaultdict(list)
    for r in rows:
        eps = to_float(r, "eps_force", None)
        if eps is None:
            continue
        grouped[eps].append(r)

    out = []
    for eps, grows in grouped.items():
        item = {
            "eps_force": eps,
            "n": len(grows),
        }
        for k in [
            "min_leaf_size",
            "rmse_phi",
            "max_abs_phi_err",
            "rmse_force_err",
            "rmse_force_err_linear",
            "force_map_r2_fixed",
            "force_map_slope",
            "force_map_rel_mae",
            "force_map_bound_violation_rate",
        ]:
            vals = [to_float(r, k, None) for r in grows]
            vals = [v for v in vals if v is not None and math.isfinite(v)]
            item[k + "_mean"] = mean(vals)
        out.append(item)
    out.sort(key=lambda r: r["eps_force"], reverse=True)
    return out


def loglog_fit(x_values, y_values):
    points = []
    for x, y in zip(x_values, y_values):
        if x is None or y is None or x <= 0.0 or y <= 0.0:
            continue
        if (not math.isfinite(x)) or (not math.isfinite(y)):
            continue
        points.append((math.log(x), math.log(y)))

    if len(points) < 2:
        return None

    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    n = len(points)
    mx = sum(xs) / float(n)
    my = sum(ys) / float(n)
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx <= 1e-16:
        return None
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = my - slope * mx
    y_pred = [intercept + slope * x for x in xs]
    sse = sum((yt - yp) ** 2 for yt, yp in zip(ys, y_pred))
    sst = sum((yt - my) ** 2 for yt in ys)
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


def fmt(x, digits=6):
    if x is None:
        return "-"
    if not isinstance(x, float):
        return str(x)
    if not math.isfinite(x):
        return "-"
    return f"{x:.{digits}g}"


def main():
    parser = argparse.ArgumentParser(description="Validate P1 criteria from M3 CSV.")
    parser.add_argument("--csv", required=True, help="Path to merged sdf_octree_pareto.csv")
    parser.add_argument("--out", default="", help="Optional markdown output path")
    parser.add_argument("--min-fit-r2", type=float, default=0.95, help="Minimum acceptable log-log fit R2")
    parser.add_argument("--min-slope", type=float, default=1.5, help="Minimum acceptable log-log slope")
    parser.add_argument("--max-slope", type=float, default=2.5, help="Maximum acceptable log-log slope")
    parser.add_argument("--max-bound-violation", type=float, default=1e-8, help="Maximum acceptable bound-violation rate")
    parser.add_argument("--min-map-r2", type=float, default=0.4, help="Minimum acceptable force-map R2")
    parser.add_argument("--max-map-rel-mae", type=float, default=0.6, help="Maximum acceptable force-map relative MAE")
    parser.add_argument("--max-rmse-ratio", type=float, default=1.0, help="Maximum acceptable rmse(actual/linear)")
    args = parser.parse_args()

    rows = read_rows(args.csv)
    if not rows:
        print("P1 validation failed: empty CSV input.")
        return 2

    agg = aggregate_by_eps(rows)
    if len(agg) < 2:
        print("P1 validation failed: insufficient eps levels for log-log fit.")
        return 2

    h_vals = [r.get("min_leaf_size_mean") for r in agg]
    rmse_vals = [r.get("rmse_phi_mean") for r in agg]
    max_vals = [r.get("max_abs_phi_err_mean") for r in agg]

    fit_rmse = loglog_fit(h_vals, rmse_vals)
    fit_max = loglog_fit(h_vals, max_vals)
    if not fit_rmse or not fit_max:
        print("P1 validation failed: cannot compute log-log fit.")
        return 2

    map_r2_values = [r.get("force_map_r2_fixed_mean") for r in agg if r.get("force_map_r2_fixed_mean") is not None]
    map_rel_mae_values = [r.get("force_map_rel_mae_mean") for r in agg if r.get("force_map_rel_mae_mean") is not None]
    bound_viol_values = [
        r.get("force_map_bound_violation_rate_mean")
        for r in agg
        if r.get("force_map_bound_violation_rate_mean") is not None
    ]
    rmse_ratio_values = []
    for r in agg:
        actual = r.get("rmse_force_err_mean")
        linear = r.get("rmse_force_err_linear_mean")
        if actual is None or linear is None or linear <= 1e-16:
            continue
        rmse_ratio_values.append(actual / linear)

    map_r2_mean = mean(map_r2_values)
    map_rel_mae_mean = mean(map_rel_mae_values)
    bound_viol_max = max(bound_viol_values) if bound_viol_values else None
    rmse_ratio_max = max(rmse_ratio_values) if rmse_ratio_values else None

    checks = []

    def add_check(name, ok, value, target):
        checks.append({"name": name, "ok": ok, "value": value, "target": target})

    add_check("fit_rmse_r2", fit_rmse["r2"] >= args.min_fit_r2, fit_rmse["r2"], f">= {args.min_fit_r2}")
    add_check("fit_max_r2", fit_max["r2"] >= args.min_fit_r2, fit_max["r2"], f">= {args.min_fit_r2}")
    add_check(
        "fit_rmse_slope",
        (fit_rmse["slope"] >= args.min_slope) and (fit_rmse["slope"] <= args.max_slope),
        fit_rmse["slope"],
        f"in [{args.min_slope}, {args.max_slope}]",
    )
    add_check(
        "fit_max_slope",
        (fit_max["slope"] >= args.min_slope) and (fit_max["slope"] <= args.max_slope),
        fit_max["slope"],
        f"in [{args.min_slope}, {args.max_slope}]",
    )
    add_check(
        "force_map_r2",
        (map_r2_mean is not None) and (map_r2_mean >= args.min_map_r2),
        map_r2_mean,
        f">= {args.min_map_r2}",
    )
    add_check(
        "force_map_rel_mae",
        (map_rel_mae_mean is not None) and (map_rel_mae_mean <= args.max_map_rel_mae),
        map_rel_mae_mean,
        f"<= {args.max_map_rel_mae}",
    )
    add_check(
        "bound_violation_rate",
        (bound_viol_max is not None) and (bound_viol_max <= args.max_bound_violation),
        bound_viol_max,
        f"<= {args.max_bound_violation}",
    )
    add_check(
        "rmse_ratio_actual_linear",
        (rmse_ratio_max is not None) and (rmse_ratio_max <= args.max_rmse_ratio),
        rmse_ratio_max,
        f"<= {args.max_rmse_ratio}",
    )

    passed = all(c["ok"] for c in checks)

    lines = []
    lines.append("# P1 Validation Report")
    lines.append("")
    lines.append("## Fitting Summary")
    lines.append("")
    lines.append(f"- rmse_phi fit: slope=`{fmt(fit_rmse['slope'])}`, R2=`{fmt(fit_rmse['r2'])}`, C=`{fmt(fit_rmse['prefactor'])}`")
    lines.append(f"- max_abs_phi_err fit: slope=`{fmt(fit_max['slope'])}`, R2=`{fmt(fit_max['r2'])}`, C=`{fmt(fit_max['prefactor'])}`")
    lines.append("")
    lines.append("## Force Mapping Summary")
    lines.append("")
    lines.append(f"- mean map_r2_fixed: `{fmt(map_r2_mean)}`")
    lines.append(f"- mean map_rel_mae: `{fmt(map_rel_mae_mean)}`")
    lines.append(f"- max bound_violation_rate: `{fmt(bound_viol_max)}`")
    lines.append(f"- max rmse(actual/linear): `{fmt(rmse_ratio_max)}`")
    lines.append("")
    lines.append("## Checks")
    lines.append("")
    lines.append("| check | value | target | pass |")
    lines.append("| --- | ---: | ---: | --- |")
    for c in checks:
        lines.append(f"| {c['name']} | {fmt(c['value'])} | {c['target']} | {'yes' if c['ok'] else 'no'} |")
    lines.append("")
    lines.append(f"Overall: **{'PASS' if passed else 'FAIL'}**")
    lines.append("")

    report = "\n".join(lines)
    print(report)

    if args.out:
        with open(args.out, "w", encoding="utf-8", newline="\n") as f:
            f.write(report)
        print(f"Generated: {args.out}")

    return 0 if passed else 3


if __name__ == "__main__":
    raise SystemExit(main())
