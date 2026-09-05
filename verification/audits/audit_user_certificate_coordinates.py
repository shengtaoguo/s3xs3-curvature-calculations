"""Coordinate checks using analytic derivatives and floating-point arithmetic."""
from __future__ import annotations
import argparse
import itertools
import json
import time
from pathlib import Path
import numpy as np
from coordinate_geometry import Jet, block, hstack, vstack, sphere, left_frame, curvature


def adjoint(p):
    w, x, y, z = (p.entry(i) for i in range(4))
    return vstack((
        hstack((w*w+x*x-y*y-z*z, 2*(x*y-w*z), 2*(x*z+w*y))),
        hstack((2*(x*y+w*z), w*w-x*x+y*y-z*z, 2*(y*z-w*x))),
        hstack((2*(x*z-w*y), 2*(y*z+w*x), w*w-x*x-y*y+z*z)),
    ))


def sym(a, b):
    return .5*(a @ b.T+b @ a.T)


def tensors(p0, q0, tp=None, tq=None):
    if tp is None:
        tp = left_frame(Jet(p0[:, None])).v
    if tq is None:
        tq = left_frame(Jet(q0[:, None])).v
    p, jp = sphere(p0, tp, 0)
    q, jq = sphere(q0, tq, 3)
    lp, lq = left_frame(p), left_frame(q)
    co = block(lp.T @ jp, np.zeros((3, 3)), np.zeros((3, 3)), lq.T @ jq)
    rp, rq = adjoint(p), adjoint(q)
    r = rp.T @ rq
    gleft = (block(29*np.eye(3), 9*r+10*np.eye(3),
                   9*r.T+10*np.eye(3), 29*np.eye(3))*(1/30)).inverse()
    basis = [Jet(np.eye(6)[:, i:i+1]) for i in range(6)]
    ah = vstack((rp.T, np.zeros((3, 3))))
    bh = vstack((np.zeros((3, 3)), rq.T))
    ac = Jet(ah.v[:, 0:1], ah.d[:, :, 0:1], ah.dd[:, :, :, 0:1])
    bc = Jet(bh.v[:, 0:1], bh.d[:, :, 0:1], bh.dd[:, :, :, 0:1])
    b = .5*sym(basis[0]+ac, basis[1]+(3/13)*basis[4])
    b -= .5*sym(basis[3]+bc, (3/13)*basis[1]+basis[4])
    c = sym(basis[0], basis[0])
    o = p.entry(0)*block(np.eye(3), -(2/7)*np.eye(3),
                        -(2/7)*np.eye(3), np.zeros((3, 3)))
    o += q.entry(1)*block(np.zeros((3, 3)), -(2/7)*np.eye(3),
                         -(2/7)*np.eye(3), np.eye(3))
    k0 = (-.2*p.entry(0)*q.entry(1))*block(np.zeros((3, 3)), np.eye(3),
                                        np.eye(3), np.zeros((3, 3)))
    k0 -= (787/5880)*p.entry(1)*q.entry(0)*np.eye(6)
    r0 = p.T @ q
    k2 = (-96/845)*(2*r0*r0-1)*gleft
    k3 = (8/5)*(b @ gleft.inverse() @ c+c @ gleft.inverse() @ b)
    def column(mat, j):
        return Jet(mat.v[:, j:j+1], mat.d[:, :, j:j+1], mat.dd[:, :, :, j:j+1])
    def local_ab(offset, right, coords):
        l1, l2, l3 = basis[offset:offset+3]
        r1, r2, r3 = [column(right, i) for i in range(3)]
        uu = l1+r1
        tt = sym(l1, l1)+sym(l1, r1)+sym(r1, r1)
        aa = sym(uu, coords.entry(0)*(l2-r2)+coords.entry(1)*((13/7)*l3-r3))
        aa -= (4/7)*coords.entry(3)*tt
        bb = sym(uu, coords.entry(0)*((13/7)*l3+r3)-coords.entry(1)*(l2+r2))
        bb -= (4/7)*coords.entry(2)*tt
        return aa, bb
    beta_p, _ = local_ab(3, bh, p)
    _, alpha_q = local_ab(0, ah, q)
    _, alpha_p = local_ab(0, ah, p)
    beta_q, _ = local_ab(3, bh, q)
    k1 = (11/325)*(beta_p+alpha_q)+(6/65)*r0*(alpha_p+beta_q)
    v = lp.T @ q
    pv = (v.T @ v)*np.eye(3)-v @ v.T
    h0 = r0*block(np.zeros((3, 3)), pv, pv, np.zeros((3, 3)))
    raw = {"g0": gleft, "B": b, "C": c, "O": o,
           "H0": h0, "K0": k0, "K1": k1, "K2": k2, "K3": k3}
    return {name: co.T @ val @ co for name, val in raw.items()}


def curvature_coefficients(metric, order=3):
    """Coordinate identity R=linear(second g)+Gamma_lower*g^-1*Gamma_lower."""
    n = 6
    zero = Jet(np.zeros((n, n)))
    metric = metric+[zero]*(order+1-len(metric))
    inverses = [np.linalg.inv(metric[0].v)]
    for k in range(1, order+1):
        inverses.append(-inverses[0] @ sum(
            (metric[i].v @ inverses[k-i] for i in range(1, k+1)),
            np.zeros((n, n))))
    lower = [.5*(g.d+g.d.transpose(1, 0, 2)-g.d.transpose(1, 2, 0))
             for g in metric]
    tensors_out = []
    for k in range(order+1):
        dd = metric[k].dd
        rk = np.empty((n, n, n, n))
        for i, j, a, b in itertools.product(range(n), repeat=4):
            rk[i, j, a, b] = .5*(dd[i, a, j, b]+dd[j, b, i, a]
                                     -dd[i, b, j, a]-dd[j, a, i, b])
        for a in range(k+1):
            for b in range(k-a+1):
                inv = inverses[k-a-b]
                rk += np.einsum("jlm,mn,ikn->ijkl", lower[a], inv, lower[b])
                rk -= np.einsum("ilm,mn,jkn->ijkl", lower[a], inv, lower[b])
        tensors_out.append(rk)
    return tensors_out


def curvature_contraction(r, a, b, c, d):
    return np.einsum("ijkl,i,j,k,l", r, a, b, c, d).item()


def reduce_path(metric):
    rr = curvature_coefficients(metric)
    eye = np.eye(6)
    x, y, normals = eye[0], eye[3], eye[[1, 2, 4, 5]]
    hessian = np.zeros((8, 8))
    source = np.zeros(8)
    for a, na in enumerate(normals):
        source[a] = 2*curvature_contraction(rr[1], na, y, y, x)
        source[a+4] = 2*curvature_contraction(rr[1], x, na, y, x)
        for b, nb in enumerate(normals):
            hessian[a, b] = 2*curvature_contraction(rr[0], na, y, y, nb)
            hessian[a+4, b+4] = 2*curvature_contraction(rr[0], x, na, nb, x)
            cross = 2*(curvature_contraction(rr[0], na, nb, y, x)
                       +curvature_contraction(rr[0], na, y, nb, x))
            hessian[a, b+4] = hessian[b+4, a] = cross
    hessian = .5*(hessian+hessian.T)
    eta = -np.linalg.solve(hessian, source)
    xs, ys = [x, eta[:4] @ normals], [y, eta[4:] @ normals]
    nums = np.zeros(4)
    for i, j, k, l in itertools.product(range(2), repeat=4):
        deg = i+j+k+l
        for n in range(deg, 4):
            nums[n] += curvature_contraction(rr[n-deg], xs[i], ys[j], ys[k], xs[l])
    grams = np.zeros((3, 4))
    for a, b in itertools.product(range(2), repeat=2):
        for n in range(a+b, 4):
            if n-a-b >= len(metric):
                continue
            mat = metric[n-a-b].v
            grams[0, n] += xs[a] @ mat @ xs[b]
            grams[1, n] += ys[a] @ mat @ ys[b]
            grams[2, n] += xs[a] @ mat @ ys[b]
    den = np.convolve(grams[0], grams[1])[:4]-np.convolve(grams[2], grams[2])[:4]
    quotient = np.zeros(4)
    for n in range(4):
        quotient[n] = (nums[n]-sum(den[i]*quotient[n-i] for i in range(1, n+1)))/den[0]
    independent_r0 = curvature(metric[0].v, metric[0].d, metric[0].dd)
    return {"sectional_coefficients": quotient.tolist(), "eta": eta.tolist(),
            "hessian_eigenvalues": np.linalg.eigvalsh(hessian).tolist(),
            "background_engine_discrepancy": float(np.max(np.abs(rr[0]-independent_r0))),
            "curvature_symmetry_residual": float(max(np.max(np.abs(rr[0]+rr[0].swapaxes(0, 1))),
                    np.max(np.abs(rr[0]+rr[0].swapaxes(2, 3))),
                    np.max(np.abs(rr[0]-rr[0].transpose(2, 3, 0, 1))))),
            "raw_numerator_coefficients": nums.tolist(), "gram_coefficients": den.tolist()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", default="cubic", choices=["cubic", "transverse", "central"])
    parser.add_argument("--kind", default="cubic")
    args, _ = parser.parse_known_args()
    if args.kind in ["transverse", "central"]:
        args.mode = args.kind
    started = time.time()
    def status(phase, **data):
        value = {"phase": phase, "elapsed_seconds": time.time()-started, **data}
        dest = args.output.with_suffix(".status.json")
        temp = dest.with_suffix(".tmp")
        temp.write_text(json.dumps(value))
        temp.replace(dest)
        print(json.dumps(value), flush=True)
    status("coordinate_jet_audit_started")
    records = []
    if args.mode == "transverse":
        zero = Jet(np.zeros((6, 6)))
        def qval(ts, h, k=zero):
            return reduce_path([ts["g0"], h, k])["sectional_coefficients"][2]
        def corrections(ts):
            oo, bb, cc = [qval(ts, ts[label]) for label in ["O", "B", "C"]]
            return np.array([
                qval(ts, ts["O"]+ts["B"], ts["K1"])-oo-bb,
                qval(ts, ts["O"]+ts["C"])-oo-cc,
                qval(ts, ts["B"], ts["K2"]),
                qval(ts, ts["B"]+ts["C"], ts["K3"])-bb-cc,
                cc,
            ])
        for phase in [0., .4]:
            theta = .9272952180016123
            p0 = np.array([np.cos(phase), np.sin(phase), 0., 0.])
            q0 = np.array([np.cos(phase+theta), np.sin(phase+theta), 0., 0.])
            tp, tq = left_frame(Jet(p0[:, None])).v, left_frame(Jet(q0[:, None])).v
            at_zero = corrections(tensors(p0, q0))
            for a, b in [(0, 2), (0, 3), (1, 2), (1, 3)]:
                vals = []
                step = 1e-5
                for t in [-step, step]:
                    rot = np.eye(4)
                    rot[a, a] = rot[b, b] = np.cos(t)
                    rot[a, b] = -np.sin(t)
                    rot[b, a] = np.sin(t)
                    vals.append(corrections(tensors(rot @ p0, rot @ q0, rot @ tp, rot @ tq)))
                derivative = (vals[1]-vals[0])/(2*step)
                records.append({"phase": phase, "rotation": [a, b],
                                "values_at_T0": at_zero.tolist(),
                                "transverse_derivatives": derivative.tolist(),
                                "step": step})
                status("transverse_rotation_complete", completed=len(records), total=8,
                       maximum_derivative_residual=float(np.max(np.abs(derivative))))
        args.output.write_text(json.dumps({"method": "independent coordinate curvature, central difference only in torus label", "records": records}, indent=2))
        status("complete", completed=8, total=8)
        return
    if args.mode == "central":
        for t in [.01, .005, .0025, -.01, -.005, -.0025]:
            p0 = np.array([1., 0., 0., 0.])
            angle = (-32/65)*t
            q0 = np.array([np.cos(angle), 0., 0., np.sin(angle)])
            ts = tensors(p0, q0)
            g = ts["g0"]+t*(ts["B"]+ts["C"])+t*t*(ts["K2"]+ts["K3"])
            rr = curvature(g.v, g.d, g.dd)
            x, y = np.eye(6)[0].copy(), np.eye(6)[3].copy()
            x += t*np.array([0., -96/65, 0., 0., -64/65, 0.])
            y += t*np.array([0., 64/65, 0., 0., 96/65, 0.])
            val = curvature_contraction(rr, x, y, y, x)/((x @ g.v @ x)*(y @ g.v @ y)-(x @ g.v @ y)**2)
            records.append({"parameter": t, "sectional_curvature": float(val), "curvature_over_t_cubed": float(val/t**3), "claimed_limit": 115712/63375})
            status("central_parameter_complete", completed=len(records), total=6)
        args.output.write_text(json.dumps({"method": "independent coordinate curvature, finite-parameter central path", "records": records}, indent=2))
        status("complete", completed=6, total=6)
        return
    for n, (p0, q0) in enumerate([
        (np.array([1., 0., 0., 0.]), np.array([3/5, 4/5, 0., 0.])),
        (np.array([3/5, 4/5, 0., 0.]), np.array([-7/25, 24/25, 0., 0.])),
        (np.array([1., 0., 0., 0.]), np.array([5/13, 12/13, 0., 0.])),
    ]):
        ts = tensors(p0, q0)
        cos2 = 2*(p0 @ q0)**2-1
        row = {"p": p0.tolist(), "q": q0.tolist(), "cos_2theta": float(cos2)}
        for label, h, k in [("B", ts["B"], Jet(np.zeros((6, 6)))),
                             ("O_K0", ts["O"], ts["K0"]),
                             ("B_C_K2_K3", ts["B"]+ts["C"], ts["K2"]+ts["K3"])]:
            row[label] = reduce_path([ts["g0"], h, k])
        row["claimed_B_second"] = float((128/845)*cos2)
        row["claimed_B_C_cubic"] = float(512/375+(9728/21125)*cos2)
        records.append(row)
        status("point_complete", completed=n+1, total=3,
               observed_cubic=row["B_C_K2_K3"]["sectional_coefficients"][3],
               claimed_cubic=row["claimed_B_C_cubic"])
    args.output.write_text(json.dumps({"method": "independent analytic coordinate jets, float64", "records": records}, indent=2))
    status("complete", completed=len(records), total=3)


if __name__ == "__main__":
    main()
