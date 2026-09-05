#!/usr/bin/env python3
"""Fail-closed entry point for the S3 x S3 curvature calculations."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "verification"
ORIGINAL = HERE / "original"
AUDITS = HERE / "audits"


def run_integrity() -> None:
    process = subprocess.run(
        [sys.executable, str(HERE / "verify_checksums.py")],
        cwd=ROOT,
        check=False,
    )
    if process.returncode:
        raise SystemExit("FAILED: repository integrity")


def numerical_diagnostics(result: Path, kind: str) -> None:
    rows = json.loads(result.read_text(encoding="utf-8"))["records"]
    if kind == "cubic":
        for row in rows:
            observed = row["B_C_K2_K3"]["sectional_coefficients"]
            if max(abs(value) for value in observed[:3]) > 1e-9:
                raise AssertionError("coordinate lower curvature coefficients disagree")
            if abs(observed[3] - row["claimed_B_C_cubic"]) > 1e-8:
                raise AssertionError("coordinate cubic disagrees")
            if abs(row["B"]["sectional_coefficients"][2] - row["claimed_B_second"]) > 1e-8:
                raise AssertionError("coordinate B quadratic disagrees")
            if abs(row["O_K0"]["sectional_coefficients"][2]) > 1e-9:
                raise AssertionError("coordinate O/K0 cancellation disagrees")
    elif kind == "transverse":
        for row in rows:
            if max(abs(value) for value in row["values_at_T0"]) > 1e-8:
                raise AssertionError("transverse correction does not vanish at T0")
            if max(abs(value) for value in row["transverse_derivatives"]) > 1e-6:
                raise AssertionError("transverse first-derivative diagnostic disagrees")
    elif kind == "central":
        closest = sorted(rows, key=lambda row: abs(row["parameter"]))[:2]
        if any(
            abs(row["curvature_over_t_cubed"] - row["claimed_limit"]) > 0.02
            for row in closest
        ):
            raise AssertionError("finite-parameter central diagnostic disagrees")
    else:
        raise ValueError(f"unknown numerical diagnostic: {kind}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--suite",
        choices=["integrity", "smoke", "certificate", "independent", "numerical", "all"],
        default="all",
    )
    args = parser.parse_args()
    if not __debug__ or sys.flags.optimize:
        parser.error("optimized Python disables certificate assertions; omit -O and -OO")

    run_integrity()
    if args.suite == "integrity":
        print("VERIFIED: repository integrity")
        return 0

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
    if args.suite == "smoke":
        tasks = [exact_point, coordinate]
    elif args.suite == "all":
        tasks = [task for group in ("certificate", "independent", "numerical") for task in groups[group]]
    else:
        tasks = groups[args.suite]

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
            (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
            raise SystemExit(f"FAILED: {name}; inspect {log}")
        if diagnostic:
            result = Path(extra[extra.index("--output") + 1])
            numerical_diagnostics(result, diagnostic)
            record["diagnostic_tolerances_passed"] = True
        print(f"PASS: {name}", flush=True)
        (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    summary["status"] = "complete"
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"VERIFIED: all selected checks completed; results: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
