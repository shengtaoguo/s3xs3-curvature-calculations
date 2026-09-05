"""Exact point and rational-function checks using the original tensor jets."""
import argparse
import json
import itertools
import time
from pathlib import Path
from types import SimpleNamespace

from exact_engine import *
from exact_cubic import curvature_series, kbbjet, moving_cubic, symprod


def addj(*jets):
    return tuple(sum((jet[k] for jet in jets), zeros(jets[0][k].shape)) for k in range(3))


def scalej(scalar, jet):
    return tuple(scalar * x for x in jet)


def zero_count(array):
    flat = np.asarray(array, dtype=object).ravel()
    bad = [str(x) for x in flat if x != 0]
    return {"total_entries": len(flat), "nonzero_entries": len(bad), "first_nonzero": bad[:3]}


def assert_zero_diagnostics(value):
    if isinstance(value,dict):
        if "nonzero_entries" in value:
            assert value["nonzero_entries"] == 0,value
        for child in value.values():
            assert_zero_diagnostics(child)


def riemann_tests(R):
    return {
        "first_pair_antisymmetry": zero_count(R + R.swapaxes(0, 1)),
        "last_pair_antisymmetry": zero_count(R + R.swapaxes(2, 3)),
        "pair_interchange": zero_count(R - R.transpose(2, 3, 0, 1)),
        "first_bianchi": zero_count(R + R.transpose(1, 2, 0, 3) + R.transpose(2, 0, 1, 3)),
    }


def contract_sparse(tensor, vectors):
    """Multilinear contraction, retaining exact zeros before products."""
    nonzero = [[(i,x) for i,x in enumerate(v) if x] for v in vectors]
    answer = zero
    for terms in itertools.product(*nonzero):
        coefficient = tensor[tuple(t[0] for t in terms)]
        if coefficient:
            for _,x in terms:
                coefficient *= x
            answer += coefficient
    return answer


def polynomial_contraction(tensors, vectors, nmax=3):
    answer = zeros(nmax+1)
    for order,tensor in enumerate(tensors):
        for degrees in itertools.product(*(range(len(v)) for v in vectors)):
            degree = order+sum(degrees)
            if degree <= nmax:
                answer[degree] += contract_sparse(tensor,[v[d] for v,d in zip(vectors,degrees)])
    return answer


def contracted_numerator(gjets, perturbation_jets, displacement, nmax=3):
    r"""Independent contracted Levi-Civita formula, no curvature_series.

    For constant-coefficient frame fields X,Y, B_ijk=g(nabla_i e_j,e_k),
    c=[X,Y], one has

      R(X,Y,Y,X) = d_X B(Y,Y,X)-d_Y B(X,Y,X)
                  +|B(X,Y,.)|^2_ginv-<B(X,X,.),B(Y,Y,.)>_ginv
                  -c.B(X,Y,.)-B(c,Y,X).

    Along the flat background torus the three B vectors vanish at order 0,
    so orders through 3 require only the first inverse-metric derivative.
    """
    jets = [gjets]+perturbation_jets
    Bs = [B_lower(j[0],j[1]) for j in jets]
    dBs = [B_lower(j[1],j[2]) for j in jets]
    identity = eye(6)
    X0,Y0 = identity[0],identity[3]
    N = identity[[1,2,4,5]]
    X = [X0,displacement[:4]@N]
    Y = [Y0,displacement[4:]@N]
    numerator = polynomial_contraction(dBs,[X,Y,Y,X],nmax)
    numerator -= polynomial_contraction(dBs,[Y,X,Y,X],nmax)

    def lower_connection(A,B):
        columns = [polynomial_contraction(Bs,[A,B,[e]],nmax) for e in identity]
        return np.stack(columns,axis=1)

    xx,xy,yy = lower_connection(X,X),lower_connection(X,Y),lower_connection(Y,Y)
    assert all(v == 0 for array in (xx[0],xy[0],yy[0]) for v in array)
    inverse = inv(gjets[0])
    first_inverse = -inverse@perturbation_jets[0][0]@inverse
    for i in range(1,nmax+1):
        for j in range(1,nmax+1-i):
            numerator[i+j] += xy[i]@inverse@xy[j]-xx[i]@inverse@yy[j]
            if i+j+1 <= nmax:
                numerator[i+j+1] += xy[i]@first_inverse@xy[j]-xx[i]@first_inverse@yy[j]
    bracket_columns = [polynomial_contraction([C[:,: ,k]],[X,Y],2) for k in range(6)]
    bracket = np.stack(bracket_columns,axis=1)
    for i in range(3):
        for j in range(nmax+1-i):
            numerator[i+j] -= bracket[i]@xy[j]
    numerator -= polynomial_contraction(Bs,[list(bracket),Y,X],nmax)
    return numerator


def symbolic_run(progress, out):
    """Prove generic-angle cubic and graph stationarity over Q(z)."""
    progress("contracted_smoke", "Validating the contracted formula at the independently replayed rational angle")
    p0,q0 = eye(4)[0],ar([Q(3,5),Q(4,5),zero,zero])
    g0j = g_jets(q0)
    base0 = SimpleNamespace(g=g0j[0],dg=g0j[1],ddg=g0j[2],ginv=inv(g0j[0]))
    b0j,c0j = killing(p0,q0,Bco),killing(p0,q0,Cco)
    h0j = addj(b0j,c0j)
    k0j = addj(kbbjet(base0,p0,q0),scalej(Q(8,5),symprod(base0,b0j,c0j)))
    eta0 = ar([Q(-96,65),Q(24,65),Q(-1216,1625),Q(288,1625),Q(1216,1625),Q(288,1625),Q(96,65),Q(24,65)])
    smoke = contracted_numerator(g0j,[h0j,k0j],eta0)/Q(15,8)
    assert all(x==zero for x in smoke[:3]) and smoke[3] == Q(1958912,1584375),[str(x) for x in smoke]
    out["checks"]["independent_contracted_rational_smoke"] = [str(x) for x in smoke]
    cosine = (1-z*z)/(1+z*z)
    sine = 2*z/(1+z*z)
    p = eye(4)[0]
    q = ar([cosine,sine,zero,zero])
    progress("symbolic_jets", "Building exact rational-function metric and tensor jets")
    gjets = g_jets(q)
    base = SimpleNamespace(g=gjets[0],dg=gjets[1],ddg=gjets[2],ginv=inv(gjets[0]))
    bj,cj = killing(p,q,Bco),killing(p,q,Cco)
    hj = addj(bj,cj)
    kj = addj(kbbjet(base,p,q),scalej(Q(8,5),symprod(base,bj,cj)))
    eta = Q(1,65)*ar([-96,32*cosine/sine,-64+24*sine*sine,24*sine*cosine,64-24*sine*sine,24*sine*cosine,96,32*cosine/sine])
    progress("symbolic_cubic", "Contracting third-order curvature directly over Q(z)")
    nums = contracted_numerator(gjets,[hj,kj],eta)
    D = Q(15,8)
    claimed = Q(512,375)+Q(9728,21125)*(cosine*cosine-sine*sine)
    checks = out["checks"]
    checks["contracted_generic_coefficients"] = [str(x/D) for x in nums]
    checks["claimed_generic_cubic"] = str(claimed)
    checks["generic_cubic_difference"] = str(nums[3]/D-claimed)
    assert all(x == zero for x in nums[:3]), checks["contracted_generic_coefficients"]
    assert nums[3]/D == claimed, checks["generic_cubic_difference"]
    checks["quadratic_graph_stationarity"] = []
    for a in range(8):
        progress("symbolic_stationarity", "Checking exact graph derivative %d of 8"%(a+1))
        e = eye(8)[a]
        plus = contracted_numerator(gjets,[hj,kj],eta+e,nmax=2)[2]
        minus = contracted_numerator(gjets,[hj,kj],eta-e,nmax=2)[2]
        derivative = (plus-minus)/(2*D)
        checks["quadratic_graph_stationarity"].append(str(derivative))
        assert derivative == zero, (a,str(derivative))
    checks["scope_warning"] = "Generic rational-function plane calculation only; singular moving-base geometry and global assembly not established by this check."


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--symbolic", action="store_true")
    parser.add_argument("--verify-first", action="store_true")
    args = parser.parse_args()
    target = Path(args.output)
    status = target.with_suffix(".status.json")
    began = time.time()

    def progress(phase, detail):
        message = {"phase": phase, "detail": detail, "elapsed_seconds": time.time() - began}
        temporary = status.with_suffix(".tmp")
        temporary.write_text(json.dumps(message), encoding="utf-8")
        temporary.replace(status)
        print(json.dumps(message), flush=True)

    out = {"scope": "Exact rational point p=1, q=(3/5,4/5,0,0); not global verification", "checks": {}}
    if args.verify_first:
        import verify
        out["scope"] = "Supplied verify.first_family() replay over Q(z), including O margin and mixed jets"
        def reported_print(*messages, **kwargs):
            message = " ".join(str(value) for value in messages)
            progress("verify_first", message)
        verify.print = reported_print
        progress("verify_first", "Entering the full supplied first-family verification")
        verify.first_family()
        out["checks"]["supplied_first_family_assertions"] = "all passed"
        out["elapsed_seconds"] = time.time()-began
        target.write_text(json.dumps(out,indent=2),encoding="utf-8")
        progress("complete", "Supplied first-family verification completed with all assertions passing")
        return
    if args.symbolic:
        out["scope"] = "Exact rational-function generic-angle check p=1, q=((1-z^2)/(1+z^2),2z/(1+z^2),0,0)"
        symbolic_run(progress,out)
        out["elapsed_seconds"] = time.time()-began
        target.write_text(json.dumps(out,indent=2),encoding="utf-8")
        progress("complete", "Rational-function cubic and eight graph-stationarity checks completed")
        return
    progress("base", "Computing exact background connection and curvature")
    p = eye(4)[0]
    q = ar([Q(3,5), Q(4,5), zero, zero])
    gj = g_jets(q)
    base = ExactBase(*gj)
    base.setup_plane()
    checks = out["checks"]
    checks["metric_cometric_formula"] = zero_count(inv(gj[0]) - Q(1,30)*np.block([[29*I3,10*I3+9*adj(q)],[10*I3+9*adj(q).T,29*I3]]))
    checks["metric_second_derivative_commutator"] = zero_count(gj[2]-gj[2].swapaxes(0,1)-eins('abk,kij->abij',C,gj[1]))
    checks["background_riemann"] = riemann_tests(base.R)
    checks["background_plane_numerator"] = str(base.R[0,3,3,0])
    checks["background_plane_gram"] = str(base.D)

    progress("linear", "Checking tensor jets, linearized curvature, and reduced quadratic coefficients")
    bj = killing(p,q,Bco)
    cj = killing(p,q,Cco)
    hj = addj(bj,cj)
    kj = addj(kbbjet(base,p,q),scalej(Q(8,5),symprod(base,bj,cj)))
    checks["B_second_derivative_commutator"] = zero_count(bj[2]-bj[2].swapaxes(0,1)-eins('abk,kij->abij',C,bj[1]))
    bd = base.linear_data(bj)
    hd = base.linear_data(hj)
    kd = base.linear_data(kj)
    checks["B_linearized_riemann"] = riemann_tests(bd[3])
    checks["B_first_curvature"] = str(bd[3][0,3,3,0]/base.D)
    checks["B_reduced_second"] = str(base.quad(bd,bd))
    checks["B_reduced_second_claim"] = str(Q(128,845)*Q(-7,25))
    checks["B_reduced_second_difference"] = str(base.quad(bd,bd)-Q(128,845)*Q(-7,25))
    checks["BC_corrected_second"] = str(base.quad(hd,hd)+kd[3][0,3,3,0]/base.D)

    progress("cubic", "Computing full third-order connection recurrence and moving plane coefficient")
    coefficients, eta, Rs, Gs = moving_cubic(base,[hj,kj])
    checks["moving_plane_coefficients"] = [str(x) for x in coefficients]
    checks["claimed_cubic"] = str(Q(512,375)+Q(9728,21125)*Q(-7,25))
    checks["cubic_difference"] = str(coefficients[3]-Q(512,375)-Q(9728,21125)*Q(-7,25))
    checks["minimizing_displacement"] = [str(x) for x in eta]
    checks["cubic_riemann"] = riemann_tests(Rs[3])
    contracted = contracted_numerator(gj,[hj,kj],eta)/base.D
    checks["independently_contracted_coefficients"] = [str(x) for x in contracted]
    checks["contracted_versus_full_connection_difference"] = [str(x-y) for x,y in zip(contracted,coefficients)]
    assert_zero_diagnostics(checks)
    assert base.R[0,3,3,0] == zero and base.D == Q(15,8)
    assert bd[3][0,3,3,0] == zero
    assert base.quad(bd,bd) == Q(128,845)*Q(-7,25)
    assert base.quad(hd,hd)+kd[3][0,3,3,0]/base.D == zero
    assert all(x == zero for x in coefficients[:3])
    assert coefficients[3] == Q(512,375)+Q(9728,21125)*Q(-7,25)
    assert all(x == y for x,y in zip(contracted,coefficients))
    out["elapsed_seconds"] = time.time()-began
    target.write_text(json.dumps(out,indent=2),encoding="utf-8")
    progress("complete", "Exact point replay complete; inspect all equalities, not only process status")


if __name__ == "__main__":
    main()
