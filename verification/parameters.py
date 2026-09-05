"""Reconstruct the coefficient choices described in the construction section."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import json
import os
from pathlib import Path
import sys
import time

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "original"))

import numpy as np
import sympy as sp
from exact_engine import Q, z, zero, eye, zeros, ar, inv, ExactBase, g_jets, killing, coord_jets, const_height
from exact_sparse import Linear, cubic_sparse
from exact_cubic import symprod, scalar_prod, scalar_times, r0jet
from exact_cross_checks import L_on_flat
from exact_h0 import h0jet, at0, deriv0
from exact_diagonal import base_hessian


def expr(value):
    return sp.cancel(value.as_expr())


def field_number(value):
    value = sp.Rational(value)
    return Q(int(value.p), int(value.q))


def round_down_one_digit(value):
    value = sp.Rational(value)
    if value <= 0:
        raise ValueError("a positive bound is required")
    unit = sp.Integer(1)
    while unit > value:
        unit /= 10
    while 10*unit <= value:
        unit *= 10
    return sp.floor(value/unit)*unit


def h0_coefficient(first_time, second_time):
    relative = ar([2*z/(1+z*z), zero, (1-z*z)/(1+z*z), zero])
    base = ExactBase(*g_jets(relative, lam=1/(1+6*first_time), tau=second_time))
    base.X, base.Y = eye(6)[0], eye(6)[3]
    base.N = eye(6)[[1,2,4,5]]
    X, Y = base.X, base.Y
    gram = (X@base.g@X)*(Y@base.g@Y)-(X@base.g@Y)**2
    h8, gradient = base_hessian(base)
    h9 = zeros((9,9))
    h9[:8,:8] = ar([[at0(value) for value in row] for row in h8])
    h9[:8,8] = h9[8,:8] = ar([deriv0(value) for value in gradient])
    h9[8,8] = deriv0(base.R[0,3,3,0],2)
    variation = Linear(base,h0jet(relative))
    shape, source = variation.data()
    source9 = np.r_[ar([at0(value) for value in source]),deriv0(variation.component(0,3,3,0))]
    shape0 = ar([[at0(value) for value in row] for row in shape])
    g0 = ar([[at0(value) for value in row] for row in base.g])
    fixed = shape0[1]@g0@shape0[1] - shape0[0]@g0@shape0[2]
    indices = [0,1,2,3,4,5,6,8]
    h_inverse = inv(h9[np.ix_(indices,indices)])
    right = source9[indices]
    solution = zeros(9)
    solution[indices] = -h_inverse@right
    assert all(value == 0 for value in h9@solution + source9)
    return sp.Rational(str((fixed-right@h_inverse@right/2)/at0(gram)))


def in_sin_squared(value, zz, denominator=None):
    """Recover a degree-two numerator and linear denominator in sin(theta)^2."""
    xx = 4*zz**2/(1+zz**2)**2
    coefficients = sp.symbols("p0:3")
    slope = sp.Symbol("denominator_slope") if denominator is None else denominator
    target = sum(coefficient*xx**i for i,coefficient in enumerate(coefficients))
    numerator = sp.together(value*(1+slope*xx)-target).as_numer_denom()[0]
    unknowns = (*coefficients, slope) if denominator is None else coefficients
    solutions = sp.solve(sp.Poly(numerator,zz).all_coeffs(), unknowns, dict=True)
    assert len(solutions) == 1 and all(key in solutions[0] for key in unknowns)
    answer = [sp.Rational(solutions[0][key]) for key in coefficients]
    return answer, sp.Rational(solutions[0][slope] if denominator is None else denominator)


def sampled_height_search(family, zz, diagonal, off_diagonal, gram):
    theta = np.linspace(0,np.pi/2,513)
    nodes = np.tan(theta/2)
    arrays = [np.stack([np.broadcast_to(np.asarray(sp.lambdify(zz,value,"numpy")(nodes),dtype=float), nodes.shape)
                        for value in row]) for row in family]
    assert all(np.isfinite(row).all() for row in arrays)
    cosine = np.cos(theta)

    def score(a,b):
        powers = np.array([1.,a,a*a])
        self_term, real_term, imaginary_term = [powers@row for row in arrays]
        v = b - diagonal*a - off_diagonal*(a*a+.25)
        real_term -= v*cosine/gram
        values = self_term-.5*np.sqrt(real_term**2+imaginary_term**2)
        if not np.isfinite(values).all():
            raise ValueError("the sampled height score is nonfinite")
        return float(np.min(values))

    best = (-np.inf,0.,0.)
    evaluations = 0
    for a in np.linspace(-1,1,81):
        for b in np.linspace(-1,1,81):
            value = score(a,b)
            evaluations += 1
            if value > best[0] + 1e-12:
                best = (value,float(a),float(b))
    step = .025
    while step > 1e-6:
        previous = best
        for da, db in ((1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)):
            a,b = previous[1]+step*da, previous[2]+step*db
            if max(abs(a),abs(b)) <= 1:
                value = score(a,b)
                evaluations += 1
                if value > best[0] + 1e-12:
                    best = (value,a,b)
        if previous == best:
            step /= 2
    assert best[0] > 0
    rational_trials = []
    for bound in range(2,33):
        a = Fraction(best[1]).limit_denominator(bound)
        b = Fraction(best[2]).limit_denominator(bound)
        value = score(float(a),float(b))
        rational_trials.append({"bound":bound,"a":str(a),"u":str(b),"sampled_margin":value})
        if value >= .75*best[0]:
            return sp.Rational(str(a)), sp.Rational(str(b)), {
                "angular_nodes":513, "box":[-1,1], "grid_points_per_axis":81,
                "pattern_step_tolerance":1e-6, "retained_fraction":.75,
                "candidate":{"a":best[1],"u":best[2],"sampled_margin":best[0]},
                "evaluations":evaluations, "rational_trials":rational_trials,
                "proof_status":"sampled search only",
            }
    raise AssertionError("no bounded-denominator approximation retained the requested margin")


def mixed_coefficients(base, p, q, a, b_data, height_data, zz):
    """Solve the transverse-jet equations in a six-coefficient tensor family."""
    cs = coord_jets(p,q)
    I, Z = eye(3), zeros((3,3))
    ap = np.block([[I,field_number(a)*I],[field_number(a)*I,Z]])
    selected = [tuple(x+field_number(a)*y for x,y in zip(pair[0],pair[1]))
                for pair in height_data]
    p3 = Linear(base,const_height(cs,3,ap)).data()
    variations = [selected[0],p3,selected[1],selected[2]]
    source = [expr(-2*base.quad(data,b_data)) for data in variations]
    ct, st = (1-zz**2)/(1+zz**2), 2*zz/(1+zz**2)
    unknowns = sp.symbols("d0 d1 d2 e0 e1 e2")
    rows = []
    axes = ((0,1,0,-1),(-1,0,1,0),(1,0,1,0),(0,1,0,1))

    def local(x0,x1,dx2,dx3,axis):
        uj,uk,wj,wk = axis
        aa = (2*(x0*(uj-wj)-x1*wk),2*x1*uk,-3*dx3)
        bb = (2*(x0*wk-x1*(uj+wj)),2*x0*uk,-3*dx2)
        return aa,bb

    for index,axis in enumerate(axes):
        dp2,dp3 = ((1,0),(0,1),(0,0),(0,0))[index]
        dq2,dq3 = ((ct,0),(0,ct),(st,0),(0,st))[index]
        ap,bp = local(1,0,dp2,dp3,axis)
        aq,bq = local(ct,st,dq2,dq3,axis)
        derivative = sum(unknowns[i]*(ap[i]+bq[i]) for i in range(3))
        derivative += ct*sum(unknowns[3+i]*(bp[i]+aq[i]) for i in range(3))
        numerator = sp.together(derivative/(2*expr(base.D))+source[index]).as_numer_denom()[0]
        rows.extend(sp.Poly(numerator,zz).all_coeffs())
    affine = sp.linsolve(rows,unknowns)
    shape = -2*a
    solutions = sp.solve(rows+[unknowns[2]-shape*unknowns[0],unknowns[5]-shape*unknowns[3]],
                         unknowns,dict=True)
    assert len(solutions) == 1 and all(key in solutions[0] for key in unknowns)
    solution = solutions[0]
    assert all(sp.expand(row.subs(solution)) == 0 for row in rows)
    return {"affine_solution":str(affine), "shape_choice":str(shape),
            "coefficients":{str(key):str(value) for key,value in solution.items()},
            "coframe_ratio":str(1+sp.Rational(3,2)*shape)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    if not __debug__:
        parser.error("optimized Python disables certificate assertions; omit -O and -OO")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output is None:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
        args.output = HERE.parent / "replay-results" / f"parameters-{stamp}.json"
    args.output.parent.mkdir(parents=True,exist_ok=True)
    began = time.monotonic()
    result = {
        "scope":"constructive reconstruction; sampled search followed by exact checks",
        "environment":{"python":sys.version, "numpy":np.__version__, "sympy":sp.__version__},
    }

    def report(phase,**data):
        result.update(data)
        result["phase"] = phase
        result["elapsed_seconds"] = time.monotonic()-began
        args.output.write_text(json.dumps(result,indent=2),encoding="utf-8")
        print(json.dumps({"phase":phase,"elapsed_seconds":result["elapsed_seconds"]}),flush=True)

    tt = sp.Symbol("t")
    first_time = sp.solve(sp.Eq(1/(sp.Rational(1,3)+2*tt),1),tt)[0]
    screen = []
    for index in range(1,5):
        second_time = sp.Rational(index,10)
        report("screen_background", background_screen=screen)
        h0 = h0_coefficient(field_number(first_time),field_number(second_time))
        screen.append({"second_time":str(second_time),"H0_coefficient":str(h0)})
        if h0 > 0:
            break
    else:
        raise AssertionError("no background passes the stated screening grid")
    gamma_bound = round_down_one_digit(h0)
    report("background_selected",first_time=str(first_time),second_time=str(second_time),
           background_screen=screen,gamma_bound=str(gamma_bound))

    ct,st = (1-z*z)/(1+z*z),2*z/(1+z*z)
    p,q = eye(4)[0],ar([ct,st,zero,zero])
    base = ExactBase(*g_jets(q,lam=1/(1+6*field_number(first_time)),tau=field_number(second_time)))
    base.setup_plane()
    zz = sp.Symbol("z")
    cos2 = expr(ct*ct-st*st)
    wdiag = sp.Rational(1,3)+first_time+second_time
    woff = first_time+second_time
    conformal_factor = 4*(wdiag-woff)
    b0 = {(0,1):Q(1,2),(1,3):Q(1,2),(6,7):Q(-1,2),(7,9):Q(-1,2)}
    b1 = {(0,7):Q(1,2),(1,6):Q(-1,2),(1,9):Q(-1,2),(3,7):Q(1,2)}
    jb0,jb1 = killing(p,q,b0),killing(p,q,b1)
    db0,db1 = Linear(base,jb0).data(),Linear(base,jb1).data()
    qb = [base.quad(db0,db0),2*base.quad(db0,db1),base.quad(db1,db1)]
    nu,kappa = sp.symbols("nu kappa")
    qexpr = sum(expr(value)*nu**i for i,value in enumerate(qb))
    numerator = sp.together(qexpr+conformal_factor*kappa*cos2).as_numer_denom()[0]
    equations = sp.Poly(numerator,zz).all_coeffs()
    solutions = sp.solve(equations,(nu,kappa),dict=True)
    assert len(solutions) == 1 and set(solutions[0]) == {nu,kappa}
    nu0,kappa0 = sp.Rational(solutions[0][nu]),sp.Rational(solutions[0][kappa])
    kernel_coeffs,kernel_slope = in_sin_squared(-expr(qb[2]),zz)
    assert kernel_slope >= 0 and all(value <= 0 for value in kernel_coeffs[1:])
    assert sum(kernel_coeffs) > 0
    affine_cos = sp.cancel((qexpr+(-expr(qb[2]))*(nu-nu0)**2)/cos2)
    assert zz not in affine_cos.free_symbols
    report("B_and_K2_recovered",nu=str(nu0),K2_coefficient=str(kappa0),
           B_cosine_coefficient=str(affine_cos),
           B_square_kernel={"numerator":[str(v) for v in kernel_coeffs],"denominator_slope":str(kernel_slope)})

    jb = tuple(x+field_number(nu0)*y for x,y in zip(jb0,jb1))
    db = Linear(base,jb).data()
    k3 = []
    for index in (1,2):
        jd = killing(p,q,{(0,index):Q(2)})
        numerator = 2*base.quad(db,Linear(base,jd).data())
        denominator = L_on_flat(symprod(base,jb,jd),base)
        k3.append(expr(-numerator/denominator))
    assert k3[0] == k3[1] and k3[0].is_Rational
    report("K3_recovered",K3_coefficient=str(k3[0]))

    I,Z = eye(3),zeros((3,3))
    ap0,aq0,off = np.block([[I,Z],[Z,Z]]),np.block([[Z,Z],[Z,I]]),np.block([[Z,I],[I,Z]])
    cs = coord_jets(p,q)
    height_data = [[Linear(base,const_height(cs,index,matrix)).data() for matrix in (diagonal,off)]
                   for index,diagonal in ((2,ap0),(6,aq0),(7,aq0))]
    pp,qq,qr = height_data
    self_coeffs = [base.quad(pp[0],pp[0]),2*base.quad(pp[0],pp[1]),base.quad(pp[1],pp[1])]

    def cross(pair):
        return [2*base.quad(pp[0],pair[0]),
                2*(base.quad(pp[0],pair[1])+base.quad(pp[1],pair[0])),
                2*base.quad(pp[1],pair[1])]

    family = [[expr(value) for value in row] for row in (self_coeffs,cross(qq),cross(qr))]
    for row,sign in zip(family,(1,-1,1)):
        assert all(sp.cancel(value.subs(zz,1/zz)-sign*value) == 0 for value in row)
    for row,sign in zip(family,(1,1,-1)):
        assert all(sp.cancel(value.subs(zz,-zz)-sign*value) == 0 for value in row)
    a,u,search = sampled_height_search(family,zz,float(wdiag),float(woff),float(expr(base.D)))
    v = sp.factor(u-wdiag*a-woff*(a*a+sp.Rational(1,4)))
    report("height_selected",height_a=str(a),K0_cross=str(u),K0_diagonal=str(v),search=search,
           height_family=[[str(value) for value in row] for row in family])

    selected = [sp.cancel(sum(value*a**i for i,value in enumerate(row))) for row in family]
    selected[1] -= v*expr(ct)/expr(base.D)
    s_coeff, slope = in_sin_squared(selected[0],zz)
    c_coeff, _ = in_sin_squared(sp.cancel(selected[1]/expr(ct)),zz,slope)
    j_coeff, _ = in_sin_squared(sp.cancel(selected[2]/expr(st)),zz,slope)
    assert slope >= 0 and all(value >= 0 for value in s_coeff)
    assert all(value >= 0 for value in c_coeff) and all(value <= 0 for value in j_coeff)
    lower = (s_coeff[0]-max(sum(c_coeff),-sum(j_coeff))/2)/(1+slope)
    assert lower > 0
    delta = round_down_one_digit(lower)
    report("height_certified",exact_height_lower_bound=str(lower),delta_bound=str(delta))

    mixed = mixed_coefficients(base,p,q,a,db,height_data,zz)
    report("mixed_correction_recovered",mixed=mixed)

    cj = killing(p,q,{(0,0):Q(1)})
    sq = scalar_prod(r0jet(p,q),r0jet(p,q))
    phi = (2*sq[0]-1,2*sq[1],2*sq[2])
    k2 = scalar_times(tuple(field_number(kappa0)*entry for entry in phi),(base.g,base.dg,base.ddg))
    bc = symprod(base,jb,cj)
    correction = tuple(x+field_number(k3[0])*y for x,y in zip(k2,bc))
    second,third,_ = cubic_sparse(base,tuple(x+y for x,y in zip(jb,cj)),correction)
    assert second == 0
    mean, harmonic = sp.symbols("mean harmonic")
    numerator = sp.together(expr(third)-mean-harmonic*cos2).as_numer_denom()[0]
    cubic_solution = sp.solve(sp.Poly(numerator,zz).all_coeffs(),(mean,harmonic),dict=True)
    assert len(cubic_solution) == 1 and cubic_solution[0][mean] > 0
    report("complete",cubic_mean=str(cubic_solution[0][mean]),
           cubic_harmonic=str(cubic_solution[0][harmonic]),
           epsilon_times_Cstar=str(delta/4),cstar_over_Mstar=str(2/gamma_bound),
           proof_status="exact coefficient identities and rational margin; sampled search is not used as proof")
    print(f"VERIFIED: parameter reconstruction; results: {args.output}",flush=True)


if __name__ == "__main__":
    main()
