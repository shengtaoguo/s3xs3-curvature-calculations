#!/usr/bin/env python3
"""Run the paper's exact and numerical curvature checks."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "verification"
ORIGINAL = HERE / "original"
AUDITS = HERE / "audits"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def finite_number(value, label: str):
    require(type(value) in (int, float), f"{label}: expected a number")
    try:
        finite = math.isfinite(value)
    except OverflowError:
        finite = False
    require(finite, f"{label}: nonfinite number")
    return value


def finite_vector(value, length: int, label: str):
    require(isinstance(value, list) and len(value) == length,
            f"{label}: expected a vector of length {length}")
    return [finite_number(item, label) for item in value]


def finite_tree(value) -> None:
    if isinstance(value, dict):
        for item in value.values():
            finite_tree(item)
    elif isinstance(value, list):
        for item in value:
            finite_tree(item)
    elif type(value) in (int, float):
        finite_number(value, "numerical result")


def match_case(values, expected, seen: set, label: str) -> int:
    matches = [i for i, case in enumerate(expected)
               if len(values) == len(case)
               and all(abs(a - b) <= 1e-12 for a, b in zip(values, case))]
    require(len(matches) == 1, f"{label}: unexpected case")
    index = matches[0]
    require(index not in seen, f"{label}: duplicate case")
    seen.add(index)
    return index


def reject_json_constant(value: str):
    raise ValueError(f"nonfinite JSON constant: {value}")


def numerical_diagnostics(result: Path, kind: str) -> None:
    data = json.loads(result.read_text(encoding="utf-8"),
                      parse_constant=reject_json_constant)
    require(isinstance(data, dict), "numerical result must be an object")
    rows = data["records"]
    counts = {"cubic": 3, "transverse": 8, "central": 6}
    require(kind in counts, f"unknown numerical diagnostic: {kind}")
    require(isinstance(rows, list) and len(rows) == counts[kind],
            f"{kind}: expected exactly {counts[kind]} records")
    require(all(isinstance(row, dict) for row in rows), "records must be objects")
    finite_tree(rows)
    seen = set()
    if kind == "cubic":
        cases = [
            (1., 0., 0., 0., 3/5, 4/5, 0., 0.),
            (3/5, 4/5, 0., 0., -7/25, 24/25, 0., 0.),
            (1., 0., 0., 0., 5/13, 12/13, 0., 0.),
        ]
        for row in rows:
            point = finite_vector(row["p"], 4, "p") + finite_vector(row["q"], 4, "q")
            index = match_case(point, cases, seen, kind)
            point = cases[index]
            cos2 = 2 * sum(a*b for a, b in zip(point[:4], point[4:]))**2 - 1
            expected_b = (128/845) * cos2
            expected_cubic = 512/375 + (9728/21125) * cos2
            require(abs(finite_number(row["cos_2theta"], "cos(2 theta)") - cos2) <= 1e-12,
                    "coordinate angle metadata disagrees")
            require(abs(finite_number(row["claimed_B_second"], "claimed quadratic") - expected_b) <= 1e-12,
                    "reported quadratic formula disagrees")
            require(abs(finite_number(row["claimed_B_C_cubic"], "claimed cubic") - expected_cubic) <= 1e-12,
                    "reported cubic formula disagrees")
            observed = finite_vector(row["B_C_K2_K3"]["sectional_coefficients"], 4, "cubic coefficients")
            b = finite_vector(row["B"]["sectional_coefficients"], 4, "B coefficients")
            odd = finite_vector(row["O_K0"]["sectional_coefficients"], 4, "O/K0 coefficients")
            require(max(abs(value) for value in observed[:3]) <= 1e-9,
                    "coordinate lower curvature coefficients disagree")
            require(abs(observed[3] - expected_cubic) <= 1e-8, "coordinate cubic disagrees")
            require(max(abs(value) for value in b[:2]) <= 1e-9,
                    "coordinate B lower coefficients disagree")
            require(abs(b[2] - expected_b) <= 1e-8, "coordinate B quadratic disagrees")
            require(max(abs(value) for value in odd[:3]) <= 1e-9,
                    "coordinate O/K0 cancellation disagrees")
    elif kind == "transverse":
        cases = [(phase, a, b) for phase in (0., .4)
                 for a, b in ((0, 2), (0, 3), (1, 2), (1, 3))]
        for row in rows:
            rotation = finite_vector(row["rotation"], 2, "rotation")
            require(all(type(i) is int for i in rotation), "rotation indices must be integers")
            phase = finite_number(row["phase"], "phase")
            match_case([phase, *rotation], cases, seen, kind)
            step = finite_number(row["step"], "transverse difference step")
            require(math.isclose(step, 1e-5, rel_tol=1e-12, abs_tol=0.),
                    "unexpected transverse difference step")
            values = finite_vector(row["values_at_T0"], 5, "values at T0")
            derivatives = finite_vector(row["transverse_derivatives"], 5, "transverse derivatives")
            require(max(abs(value) for value in values) <= 1e-8,
                    "transverse correction does not vanish at T0")
            require(max(abs(value) for value in derivatives) <= 1e-6,
                    "transverse first-derivative diagnostic disagrees")
    elif kind == "central":
        cases = [(t,) for t in (.01, .005, .0025, -.01, -.005, -.0025)]
        limit = 115712/63375
        for row in rows:
            parameter = finite_number(row["parameter"], "central parameter")
            index = match_case([parameter], cases, seen, kind)
            value = finite_number(row["sectional_curvature"], "central curvature")
            scaled = finite_number(row["curvature_over_t_cubed"], "scaled central curvature")
            claimed = finite_number(row["claimed_limit"], "claimed central limit")
            require(abs(claimed - limit) <= 1e-12, "reported central limit disagrees")
            require(math.isclose(scaled, value / parameter**3, rel_tol=1e-10, abs_tol=1e-10),
                    "central curvature scaling disagrees")
            if abs(cases[index][0]) == .0025:
                require(abs(scaled - limit) <= .02,
                        "finite-parameter central diagnostic disagrees")
    require(len(seen) == counts[kind], f"{kind}: incomplete case coverage")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--suite",
        choices=["smoke", "certificate", "independent", "numerical", "validation", "all"],
        default="all",
    )
    args = parser.parse_args()
    if not __debug__ or sys.flags.optimize:
        parser.error("optimized Python disables certificate assertions; omit -O and -OO")

    try:
        import numpy
        import sympy
    except ImportError as error:
        raise SystemExit(
            "missing dependency; run: python3 -m pip install -r verification/requirements.txt"
        ) from error

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    output = ROOT / "replay-results" / stamp
    output.mkdir(parents=True, exist_ok=False)

    exact_point = (
        "exact-rational-point",
        AUDITS / "audit_user_certificate_engine.py",
        ["--output", str(output / "rational.json")],
        None,
    )
    coordinate = (
        "coordinate-cubic",
        AUDITS / "audit_user_certificate_coordinates.py",
        ["--output", str(output / "coordinate-cubic.json"), "--mode", "cubic"],
        "cubic",
    )
    groups = {
        "certificate": [
            ("original-full-verifier", ORIGINAL / "verify.py", ["--part", "all"], None),
            ("original-supplemental", ORIGINAL / "audit_checks.py", [], None),
            ("original-odd-center", ORIGINAL / "audit_diagonal_odd.py", [], None),
        ],
        "independent": [
            exact_point,
            (
                "contracted-rational-function",
                AUDITS / "audit_user_certificate_engine.py",
                ["--symbolic", "--output", str(output / "symbolic.json")],
                None,
            ),
            (
                "smooth-diagonal-path",
                AUDITS / "audit_user_certificate_center.py",
                ["--output", str(output / "central.json")],
                None,
            ),
            (
                "quadratic-point-audit",
                AUDITS / "audit_user_certificate_quadratic.py",
                ["--kind", "points", "--output", str(output / "quadratic.json")],
                None,
            ),
        ],
        "numerical": [
            coordinate,
            (
                "coordinate-transverse",
                AUDITS / "audit_user_certificate_coordinates.py",
                ["--output", str(output / "coordinate-transverse.json"), "--mode", "transverse"],
                "transverse",
            ),
            (
                "coordinate-central",
                AUDITS / "audit_user_certificate_coordinates.py",
                ["--output", str(output / "coordinate-central.json"), "--mode", "central"],
                "central",
            ),
        ],
    }
    groups["validation"] = [groups["independent"][-1], *groups["numerical"]]
    if args.suite == "smoke":
        tasks = [exact_point, coordinate]
    elif args.suite == "all":
        tasks = [task for group in ("certificate", "independent", "numerical") for task in groups[group]]
    else:
        tasks = groups[args.suite]
    if args.suite in ("all", "validation"):
        tasks.append(("failure-detection-tests", HERE / "test_checks.py",
                      ["--results", str(output)], None))

    environment = os.environ.copy()
    environment.pop("PYTHONOPTIMIZE", None)
    environment["PYTHONPATH"] = os.pathsep.join([str(ORIGINAL), str(AUDITS)])
    environment["PYTHONUNBUFFERED"] = "1"
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment["PYTHONHASHSEED"] = "0"
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        environment[key] = "1"

    summary = {
        "suite": args.suite,
        "python": sys.version,
        "numpy": numpy.__version__,
        "sympy": sympy.__version__,
        "status": "running",
        "checks": [],
    }
    for name, script, extra, diagnostic in tasks:
        print(f"START: {name}", flush=True)
        began = time.monotonic()
        log = output / f"{name}.log"
        with log.open("w", encoding="utf-8") as handle:
            process = subprocess.Popen(
                [sys.executable, str(script), *extra],
                cwd=ROOT,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            assert process.stdout is not None
            for line in process.stdout:
                print(line, end="", flush=True)
                handle.write(line)
                handle.flush()
            returncode = process.wait()
        record = {
            "name": name,
            "exit_code": returncode,
            "elapsed_seconds": time.monotonic() - began,
        }
        summary["checks"].append(record)
        if returncode:
            summary["status"] = "failed"
            record["status"] = "failed"
            (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
            raise SystemExit(f"FAILED: {name}; inspect {log}")
        if diagnostic:
            result = Path(extra[extra.index("--output") + 1])
            try:
                numerical_diagnostics(result, diagnostic)
            except (AssertionError, KeyError, TypeError, ValueError, OSError) as error:
                summary["status"] = "failed"
                record["status"] = "failed"
                record["error"] = str(error)
                (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
                raise SystemExit(f"FAILED: {name}; {error}") from error
            record["diagnostic_tolerances_passed"] = True
        record["status"] = "passed"
        print(f"PASS: {name}", flush=True)
        (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    summary["status"] = "complete"
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"VERIFIED: all selected checks completed; results: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
