"""Pole-free exact central-path audit of the user curvature certificate.

Run only on the configured external worker under a resource lease.  The
curvature recurrence is reused from the supplied certificate, but the smooth
Cayley path and coefficient-composition check below are new and independent of
its singular fixed-base optimizer and its diagonal contraction routine.
"""
import argparse
import itertools
import json
import runpy
import time
from pathlib import Path
from types import SimpleNamespace

from exact_engine import *
from exact_cubic import curvature_series, kbbjet, symprod
from exact_diagonal import coef


def array_coefficient(array, degree):
    if degree < 0:
        return zeros(array.shape)
    return ar([coef(x, degree) for x in array.ravel()]).reshape(array.shape)


def contraction(tensor, vectors):
    terms = [[(i, x) for i, x in enumerate(v) if x] for v in vectors]
    value = zero
    for indices in itertools.product(*terms):
        coefficient = tensor[tuple(item[0] for item in indices)]
        for _, scalar in indices:
            coefficient *= scalar
        value += coefficient
    return value


def polynomial_contraction(tensors, vectors, nmax=3):
    result = zeros(nmax + 1)
    for order, tensor in enumerate(tensors):
        for degrees in itertools.product(*(range(len(v)) for v in vectors)):
            degree = order + sum(degrees)
            if degree <= nmax:
                result[degree] += contraction(tensor, [v[d] for v, d in zip(vectors, degrees)])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    target = Path(args.output)
    status = target.with_suffix(".status.json")
    began = time.time()

    def progress(phase, detail):
        data = {"phase": phase, "detail": detail, "elapsed_seconds": time.time() - began}
        temp = status.with_suffix(".tmp")
        temp.write_text(json.dumps(data), encoding="utf-8")
        temp.replace(status)
        print(json.dumps(data), flush=True)

    out = {"scope": "Exact central pole-free path and supplied full-normal diagonal replay", "checks": {}}
    progress("central_cayley_jets", "Composing spatial jets with a smooth Cayley path at p=q=1")
    p = eye(4)[0]
    a = Q(16, 65)
    denominator = 1 + a * a * z * z
    q = ar([(1 - a * a * z * z) / denominator, zero, zero, -2 * a * z / denominator])
    gj = g_jets(q)
    base = SimpleNamespace(g=gj[0], dg=gj[1], ddg=gj[2], ginv=inv(gj[0]))
    bj, cj = killing(p, q, Bco), killing(p, q, Cco)
    hj = tuple(x + y for x, y in zip(bj, cj))
    kj = tuple(x + Q(8, 5) * y for x, y in zip(kbbjet(base, p, q), symprod(base, bj, cj)))
    combined = []
    for degree in range(4):
        progress("central_series_extract", "Extracting exact total metric/spatial-jet order %d" % degree)
        combined.append(tuple(array_coefficient(gj[k], degree)
                              + array_coefficient(hj[k], degree - 1)
                              + array_coefficient(kj[k], degree - 2)
                              for k in range(3)))
    fixed = ExactBase(*combined[0])
    progress("central_curvature_recurrence", "Computing the full curvature series after exact spatial-jet composition")
    Rs, _ = curvature_series(fixed, combined[1:], nmax=3, verbose=True)
    standard = eye(6)
    X = [standard[0], Q(1, 65) * (-96 * standard[1] - 64 * standard[4])]
    Y = [standard[3], Q(1, 65) * (64 * standard[1] + 96 * standard[4])]
    numerator = polynomial_contraction(Rs, [X, Y, Y, X])
    D = fixed.g[0, 0] * fixed.g[3, 3] - fixed.g[0, 3] ** 2
    assert D == Q(15, 8)
    assert all(v == 0 for v in numerator[:3]), [str(v) for v in numerator]
    cubic = numerator[3] / D
    assert cubic == Q(115712, 63375), str(cubic)
    out["checks"]["pole_free_central_curvature_coefficients"] = [str(v / D) for v in numerator]
    out["checks"]["claimed_cubic_difference"] = str(cubic - Q(115712, 63375))

    progress("central_graph_stationarity", "Checking all eight fixed-base graph directions at the smooth optimizer")
    normal = standard[[1, 2, 4, 5]]
    derivatives = []
    for i in range(8):
        xp, yp = [x.copy() for x in X], [y.copy() for y in Y]
        xm, ym = [x.copy() for x in X], [y.copy() for y in Y]
        if i < 4:
            xp[1] += normal[i]
            xm[1] -= normal[i]
        else:
            yp[1] += normal[i - 4]
            ym[1] -= normal[i - 4]
        plus = polynomial_contraction(Rs[:3], [xp, yp, yp, xp], nmax=2)[2]
        minus = polynomial_contraction(Rs[:3], [xm, ym, ym, xm], nmax=2)[2]
        derivative = (plus - minus) / (2 * D)
        assert derivative == 0, (i, str(derivative))
        derivatives.append(str(derivative))
    out["checks"]["eight_graph_stationarity_derivatives"] = derivatives

    progress("supplied_diagonal", "Replaying supplied complete ten-variable/full-normal diagonal routine")
    from verify import diagonal
    diagonal()
    out["checks"]["supplied_full_normal_diagonal"] = "pass"
    progress("supplied_diagonal_odd", "Replaying supplied full-normal diagonal odd-height routine")
    runpy.run_path(str(Path(__file__).resolve().parents[1] / "original" / "audit_diagonal_odd.py"), run_name="__main__")
    out["checks"]["supplied_full_normal_diagonal_odd"] = "pass"
    out["elapsed_seconds"] = time.time() - began
    target.write_text(json.dumps(out, indent=2), encoding="utf-8")
    progress("complete", "Exact pole-free central cubic, graph stationarity, and both full-normal diagonal replays passed")


if __name__ == "__main__":
    main()
