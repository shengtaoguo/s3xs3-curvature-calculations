"""Check that real successful outputs pass and corrupted copies are rejected."""

from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE / "original"))
sys.path.insert(0, str(HERE / "audits"))
spec = importlib.util.spec_from_file_location("curvature_runner", HERE / "verify.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
import audit_checks
import audit_user_certificate_quadratic as quadratic

RESULTS = None
REJECTIONS = (AssertionError, KeyError, TypeError, ValueError)


class AcceptanceChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.numerical = {
            mode: json.loads((RESULTS / f"coordinate-{mode}.json").read_text(encoding="utf-8"))
            for mode in ("cubic", "transverse", "central")
        }
        cls.quadratic = json.loads((RESULTS / "quadratic.json").read_text(encoding="utf-8"))["results"]

    def validate_numerical(self, data, mode):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "result.json"
            path.write_text(json.dumps(data), encoding="utf-8")
            runner.numerical_diagnostics(path, mode)

    def test_actual_outputs_pass(self):
        for mode, data in self.numerical.items():
            with self.subTest(mode=mode):
                self.validate_numerical(data, mode)
        quadratic.validate_quadratic_results(self.quadratic)

    def test_numerical_case_coverage(self):
        for mode, original in self.numerical.items():
            for fault in ("empty", "missing", "duplicate", "unexpected"):
                with self.subTest(mode=mode, fault=fault):
                    data = deepcopy(original)
                    rows = data["records"]
                    if fault == "empty":
                        rows.clear()
                    elif fault == "missing":
                        rows.pop()
                    elif fault == "duplicate":
                        rows[-1] = deepcopy(rows[0])
                    elif mode == "cubic":
                        rows[0]["p"][0] = 2.
                    elif mode == "transverse":
                        rows[0]["phase"] = .1
                    else:
                        rows[0]["parameter"] = .003
                    with self.assertRaises(REJECTIONS):
                        self.validate_numerical(data, mode)

    def test_nonfinite_and_nonnumeric_values(self):
        for mode, original in self.numerical.items():
            for value in (float("nan"), float("inf"), -float("inf"), "NaN", False, None):
                with self.subTest(mode=mode, value=value):
                    data = deepcopy(original)
                    row = data["records"][0]
                    if mode == "cubic":
                        row["B_C_K2_K3"]["sectional_coefficients"][3] = value
                    elif mode == "transverse":
                        row["transverse_derivatives"][0] = value
                    else:
                        row["curvature_over_t_cubed"] = value
                    with self.assertRaises(REJECTIONS):
                        self.validate_numerical(data, mode)

    def test_vector_lengths(self):
        for mode in ("cubic", "transverse"):
            with self.subTest(mode=mode):
                data = deepcopy(self.numerical[mode])
                row = data["records"][0]
                values = (row["B_C_K2_K3"]["sectional_coefficients"]
                          if mode == "cubic" else row["transverse_derivatives"])
                values.pop()
                with self.assertRaises(REJECTIONS):
                    self.validate_numerical(data, mode)

    def test_wrong_numerical_results(self):
        for mode, original in self.numerical.items():
            with self.subTest(mode=mode):
                data = deepcopy(original)
                if mode == "cubic":
                    data["records"][0]["B_C_K2_K3"]["sectional_coefficients"][3] += 1
                elif mode == "transverse":
                    data["records"][0]["transverse_derivatives"][0] = 1.
                else:
                    row = min(data["records"], key=lambda item: abs(item["parameter"]))
                    row["curvature_over_t_cubed"] = row["claimed_limit"] + 1
                    row["sectional_curvature"] = row["parameter"]**3 * row["curvature_over_t_cubed"]
                with self.assertRaises(REJECTIONS):
                    self.validate_numerical(data, mode)

    def test_reported_reference_is_not_trusted(self):
        data = deepcopy(self.numerical["cubic"])
        data["records"][0]["claimed_B_C_cubic"] += 1
        data["records"][0]["B_C_K2_K3"]["sectional_coefficients"][3] += 1
        with self.assertRaises(REJECTIONS):
            self.validate_numerical(data, "cubic")
        data = deepcopy(self.numerical["central"])
        for row in data["records"]:
            row["claimed_limit"] += 1
            row["curvature_over_t_cubed"] += 1
            row["sectional_curvature"] = row["parameter"]**3 * row["curvature_over_t_cubed"]
        with self.assertRaises(REJECTIONS):
            self.validate_numerical(data, "central")

    def test_quadratic_case_coverage(self):
        for fault in ("empty", "missing", "duplicate"):
            with self.subTest(fault=fault):
                rows = deepcopy(self.quadratic)
                if fault == "empty":
                    rows.clear()
                elif fault == "missing":
                    rows.pop()
                else:
                    rows[-1] = deepcopy(rows[0])
                with self.assertRaises(REJECTIONS):
                    quadratic.validate_quadratic_results(rows)

    def test_quadratic_equalities_and_inequalities(self):
        faults = (
            ("discrepancy", "-1"), ("gram_recurrence_discrepancy", "1"),
            ("Q_B_discrepancy", "1"), ("corrected_odd", "1"),
            ("margin_over_claimed_bound", "-1"), ("rho", "-1"),
            ("hessian_ldlt_pivots", ["1"]*7),
            ("hessian_ldlt_pivots", ["1"]*7 + ["0"]),
            ("hessian_ldlt_pivots", ["1"]*7 + ["-1"]),
            ("first_variations", {"O": "1", "B": "0", "C": "0"}),
        )
        for key, value in faults:
            with self.subTest(field=key, value=value):
                rows = deepcopy(self.quadratic)
                rows[0][key] = value
                with self.assertRaises(REJECTIONS):
                    quadratic.validate_quadratic_results(rows)

    def test_quadratic_missing_geometry(self):
        for index, original in enumerate(self.quadratic):
            for field in ("p", "q", "axis"):
                with self.subTest(case=original["case"], field=field):
                    rows = deepcopy(self.quadratic)
                    del rows[index][field]
                    with self.assertRaises(REJECTIONS):
                        quadratic.validate_quadratic_results(rows)

    def test_quadratic_invalid_quaternion_coordinates(self):
        faults = (
            None, "1,0,0,0", [], [1, 0, 0], [1, 0, 0, 0, 0],
            [0, 0, 0, 0], [2, 0, 0, 0], [1., 0, 0, 0],
            [True, 0, 0, 0], [None, 0, 0, 0], ["1/0", 0, 0, 0],
            ["NaN", 0, 0, 0], ["1e999", 0, 0, 0],
            [float("nan"), 0, 0, 0], [float("inf"), 0, 0, 0],
        )
        for index, original in enumerate(self.quadratic):
            for field in ("p", "q"):
                for value in faults:
                    with self.subTest(case=original["case"], field=field, value=value):
                        rows = deepcopy(self.quadratic)
                        rows[index][field] = deepcopy(value)
                        with self.assertRaises(REJECTIONS):
                            quadratic.validate_quadratic_results(rows)

    def test_quadratic_geometry_must_match_case(self):
        for index, original in enumerate(self.quadratic):
            for field in ("p", "q"):
                with self.subTest(case=original["case"], field=field):
                    rows = deepcopy(self.quadratic)
                    # The antipode is still a unit quaternion, but is not
                    # the geometric input prescribed for this case.
                    rows[index][field] = [str(-Fraction(value)) for value in original[field]]
                    with self.assertRaises(REJECTIONS):
                        quadratic.validate_quadratic_results(rows)

    def test_quadratic_invalid_or_mismatched_axis(self):
        for index, original in enumerate(self.quadratic):
            faults = (None, False, 0., "0", -1, 3, 99, (original["axis"] + 1) % 3)
            for value in faults:
                with self.subTest(case=original["case"], axis=value):
                    rows = deepcopy(self.quadratic)
                    rows[index]["axis"] = value
                    with self.assertRaises(REJECTIONS):
                        quadratic.validate_quadratic_results(rows)

    def test_quadratic_relabelled_clone_is_rejected(self):
        canonical = next(row for row in self.quadratic if row["case"] == "canonical-i")
        rows = [deepcopy(row) if row["case"] == "type-II-first-null"
                else dict(deepcopy(canonical), case=row["case"])
                for row in self.quadratic]
        with self.assertRaises(REJECTIONS):
            quadratic.validate_quadratic_results(rows)

    def test_quadratic_equivalent_rational_geometry_passes(self):
        rows = deepcopy(self.quadratic)
        for row in rows:
            for field in ("p", "q"):
                values = [Fraction(value) for value in row[field]]
                row[field] = [f"{2 * value.numerator}/{2 * value.denominator}" for value in values]
        quadratic.validate_quadratic_results(rows)

    def test_injected_repeated_point_exits_nonzero(self):
        program = """
import sys
import audit_user_certificate_quadratic as audit
original = audit.point_case
def repeated_point(name):
    return original('canonical-i' if name == 'canonical-j' else name)
audit.point_case = repeated_point
sys.argv = ['quadratic-point-fault', '--kind', 'points', '--output', sys.argv[1]]
audit.main()
"""
        with tempfile.TemporaryDirectory() as temporary:
            environment = os.environ.copy()
            environment["PYTHONPATH"] = os.pathsep.join(
                [str(HERE / "original"), str(HERE / "audits")]
            )
            process = subprocess.run(
                [sys.executable, "-c", program, str(Path(temporary) / "wrong.json")],
                cwd=ROOT, env=environment, capture_output=True, text=True, timeout=180,
            )
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("q does not match quadratic case canonical-j", process.stdout + process.stderr)
        self.assertNotIn("all six exact quadratic cases passed", process.stdout)

    def test_injected_quadratic_formula_exits_nonzero(self):
        program = """
import sys
import audit_user_certificate_quadratic as audit
original = audit.claimed_odd
def wrong_formula(*args):
    rho, expected, reversed_j = original(*args)
    return rho, expected + audit.Q(1), reversed_j
audit.claimed_odd = wrong_formula
sys.argv = ['quadratic-fault', '--kind', 'points', '--output', sys.argv[1]]
audit.main()
"""
        with tempfile.TemporaryDirectory() as temporary:
            environment = os.environ.copy()
            environment["PYTHONPATH"] = os.pathsep.join(
                [str(HERE / "original"), str(HERE / "audits")]
            )
            process = subprocess.run(
                [sys.executable, "-c", program, str(Path(temporary) / "wrong.json")],
                cwd=ROOT, env=environment, capture_output=True, text=True, timeout=180,
            )
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("quadratic discrepancy is nonzero", process.stdout + process.stderr)
        self.assertNotIn("all six exact quadratic cases passed", process.stdout)

    def test_formal_K1_rejects_wrong_table(self):
        audit_checks.check_K1_formal()
        changed = audit_checks.boco.copy()
        changed[0, 0] += audit_checks.Q(1)
        with patch.object(audit_checks, "boco", changed):
            with self.assertRaisesRegex(AssertionError, "formal polynomials"):
                audit_checks.check_K1_formal()


def main():
    global RESULTS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    RESULTS = parser.parse_args().results
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(AcceptanceChecks)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
