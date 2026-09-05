"""Exact rational point audit of the 2026-09-05 user certificate.

This is an audit, not a positivity certificate.  It uses the supplied ordered
left-frame engine, whose derivatives and convention are inspected separately,
and compares its Gram quadratic formula with the full curvature recurrence.
Run only on the configured external compute resource.
"""
from __future__ import annotations

import argparse
import json
import runpy
import time
from pathlib import Path

import numpy as np
from exact_engine import (
    Q, ar, zeros, eye, eins, inv, adj, qmul, coord_jets, killing,
    Bco, Cco, AP, AQ, ExactBase, g_jets, const_height, zero, one,
)
from exact_cubic import scalar_prod, scalar_times, curvature_series


def addj(*jets):
    return tuple(sum((j[k] for j in jets), zeros(jets[0][k].shape)) for k in range(3))


def k0jet(p, q):
    cs = coord_jets(p, q)
    h = [tuple(x[..., j] for x in cs) for j in range(8)]
    off = np.block([[zeros((3, 3)), eye(3)], [eye(3), zeros((3, 3))]])
    constant = lambda a: (a, zeros((6, 6, 6)), zeros((6, 6, 6, 6)))
    p0q1 = scalar_prod(h[0], h[5])
    p1q0 = scalar_prod(h[1], h[4])
    return addj(
        scalar_times(p0q1, constant(-Q(1, 5) * off)),
        scalar_times(p1q0, constant(-Q(787, 5880) * eye(6))),
    )


def first(base, jet, axis):
    _, r1 = base.linear(jet)
    i, j = axis, 3 + axis
    denom = base.g[i, i] * base.g[j, j] - base.g[i, j] ** 2
    return r1[i, j, j, i] / denom


def tensorjets(p, q):
    cs = coord_jets(p, q)
    return {
        "O": addj(const_height(cs, 0, AP), const_height(cs, 5, AQ)),
        "B": killing(p, q, Bco),
        "C": killing(p, q, Cco),
        "K0": k0jet(p, q),
    }


def claimed_odd(p, axis, costheta, sintheta):
    # Oriented adapted normal basis (u,v,w) is a cyclic permutation of i,j,k.
    vi, wi = (axis + 1) % 3, (axis + 2) % 3
    pv, pw = qmul(p, eye(4)[1 + vi]), qmul(p, eye(4)[1 + wi])
    e, f = ar([pv[0], pw[0]]), ar([pv[1], pw[1]])
    rho = e @ e + f @ f
    inner = e @ f
    jinner = -e[1] * f[0] + e[0] * f[1]
    x = sintheta * sintheta
    den = Q(31752000) * (25 + 3 * x)
    pol = Q(12632468) + Q(1767741) * x + Q(5670) * x * x
    vv = Q(3941932) + Q(686019) * x
    ww = Q(3521524) + Q(636555) * x
    base = pol * rho / (2 * den) + vv * costheta * inner / den
    return rho, base - ww * sintheta * jinner / den, base + ww * sintheta * jinner / den


def ldlt_pivots(a):
    a = a.copy()
    values = []
    for k in range(len(a)):
        pivot = a[k, k]
        values.append(pivot)
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                a[i, j] -= a[i, k] * a[k, j] / pivot
    return values


def clean(v):
    if isinstance(v, dict):
        return {k: clean(x) for k, x in v.items()}
    if isinstance(v, (tuple, list, np.ndarray)):
        return [clean(x) for x in v]
    if isinstance(v, (str, int, float, bool)) or v is None:
        return v
    return str(v)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--kind", default="points")
    args, _ = parser.parse_known_args()
    out = Path(args.output)
    status = out.with_suffix(".status.json")
    start = time.time()
    results = []

    def save(phase, detail):
        record = dict(phase=phase, detail=detail, elapsed_seconds=time.time() - start,
                      completed=len(results), heartbeat=time.time(), expected_finish="unknown")
        tmp = status.with_suffix(".tmp")
        tmp.write_text(json.dumps(record), encoding="utf-8")
        tmp.replace(status)
        print(json.dumps(record), flush=True)

    def publish(rec):
        results.append(clean(rec))
        out.write_text(json.dumps(dict(results=results), indent=2), encoding="utf-8")
        print(json.dumps(clean(rec)), flush=True)

    save("setup", "checking product convention")
    prod = ExactBase(3 * eye(6), zeros((6, 6, 6)), zeros((6, 6, 6, 6)))
    assert prod.R[0, 1, 1, 0] == 3 and prod.R[0, 3, 3, 0] == 0
    if args.kind == "second-audit":
        import verify
        save("second-family", "complete base-and-plane H0 reduction")
        verify.second_family()
        publish(dict(case="complete-second-family", passed=True))
        save("compact-and-recurrence", "compact K1 and independent full-connection cubic")
        runpy.run_path(str(Path(__file__).resolve().parents[1] / "original" / "audit_checks.py"), run_name="__main__")
        publish(dict(case="compact-K1-and-full-connection-cubic", passed=True))
        save("complete", "all selected complete-bundle checks passed")
        return
    cases = [
        ("canonical-i", eye(4)[0], 0, Q(3, 5), Q(4, 5)),
        ("canonical-j", eye(4)[0], 1, Q(3, 5), Q(4, 5)),
        ("left-translate-i", ar([Q(1, 2)] * 4), 0, Q(3, 5), Q(4, 5)),
        ("tilted-i", ar([Q(3, 5), zero, Q(4, 5), zero]), 0, Q(5, 13), Q(12, 13)),
        ("negative-angle-j", eye(4)[0], 1, Q(3, 5), Q(-4, 5)),
    ]
    for name, p, axis, a, b in cases:
        save("quadratic", name)
        r = a * eye(4)[0] + b * eye(4)[1 + axis]
        q = qmul(p, r)
        base = ExactBase(*g_jets(r))
        base.setup_plane(axis)
        i, j = axis, 3 + axis
        assert base.R[i, j, j, i] == 0
        assert base.D == Q(15, 8)
        # Validate curvature tensor identities before extracting the reduction.
        assert all(t == 0 for t in (base.R + base.R.swapaxes(0, 1)).flat)
        assert all(t == 0 for t in (base.R + base.R.swapaxes(2, 3)).flat)
        jets = tensorjets(p, q)
        firsts = {n: first(base, jets[n], axis) for n in ("O", "B", "C")}
        odddata = base.linear_data(jets["O"])
        qo = base.quad(odddata, odddata)
        lk0 = first(base, jets["K0"], axis)
        actual = qo + lk0
        rho, expected, expected_reverse_j = claimed_odd(p, axis, a, b)
        rec = dict(case=name, p=p, q=q, axis=axis, first_variations=firsts,
                   hessian_ldlt_pivots=ldlt_pivots(base.H),
                   Q_O=qo, L_K0=lk0, corrected_odd=actual, rho=rho,
                   claimed_oriented=expected, claimed_reverse_J=expected_reverse_j,
                   discrepancy=actual-expected, discrepancy_reverse_J=actual-expected_reverse_j,
                   margin_over_claimed_bound=actual-rho/Q(250))
        # Full exact recurrence is separate from the Gram-term implementation.
        if name == "canonical-i":
            rs, _ = curvature_series(base, [jets["O"]], nmax=2)
            rr = odddata[1]
            recurrence = (rs[2][i, j, j, i] - (rr @ base.Hinv @ rr)/2) / base.D
            rec["Q_O_full_recurrence"] = recurrence
            rec["gram_recurrence_discrepancy"] = recurrence - qo
        if name == "canonical-i":
            bj = base.linear_data(jets["B"])
            qb = base.quad(bj, bj)
            rec["Q_B"] = qb
            rec["Q_B_claimed"] = Q(128, 845) * (a*a-b*b)
            rec["Q_B_discrepancy"] = qb - rec["Q_B_claimed"]
        publish(rec)
    save("first-null", "type-II rational representative")
    p, q, axis = eye(4)[0], eye(4)[2], 0
    base = ExactBase(*g_jets(q))
    jets = tensorjets(p, q)
    publish(dict(case="type-II-first-null", first_variations={
        n: first(base, jets[n], axis) for n in ("O", "B", "C")
    }))
    save("complete", "exact rational audit complete; identities may pass or fail")


if __name__ == "__main__":
    main()
