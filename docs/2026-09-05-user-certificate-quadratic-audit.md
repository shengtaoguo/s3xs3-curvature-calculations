# Quadratic audit of the supplied September 5 certificate

Target: `EXT-POSCURV-S3XS3`.  This note records an independent audit of
selected assertions in the user-supplied `PROOF.md`, not an acceptance of its
conclusion.  The local source copy is in
`references/user-certificate-2026-09-05/`.

## First-nullity: geometric verification

Let a surviving flat torus be parametrized by
\((p(s),q(t))=(p_*e^{su},q_*e^{tu})\), with \(|u|=1\).  Its coordinate
vector fields are \(X=E_u\) and \(Y=F_u\).  Their left Maurer--Cartan
components are constant, and the right components are constant as well:
\[
\operatorname{Ad}_{p_*e^{su}}u=\operatorname{Ad}_{p_*}u,
\qquad
\operatorname{Ad}_{q_*e^{tu}}u=\operatorname{Ad}_{q_*}u.
\]
Consequently, the restrictions of the displayed tensors \(B\) and \(C\)
to the tangent bundle of this torus have constant coefficients.  For \(O\),
the tangent coefficients are
\[
 O(X,X)=p_0(s),\quad O(Y,Y)=q_1(t),\quad
 O(X,Y)=-\frac27\bigl(p_0(s)+q_1(t)\bigr).
\]
Thus all three linearized Gaussian curvatures vanish.  Since the background
torus is totally geodesic, the ambient and intrinsic first variations agree.

For \(H_0\), a plus zero plane has \(v\parallel u\), hence
\((|v|^2I-vv^T)u=0\).  At a minus zero plane, \(r_0=0\).  Therefore
\(H_0|_{TT}=0\) on every surviving flat torus, and its first curvature
variation vanishes there too.  This verifies global first-nullity on both
zero families, independently of the coefficient tables.

## Exact computational audit method

The audit script is `scripts/audit_user_certificate_quadratic.py`.
Substantive computation runs only on `windows`, with a one-core/one-GiB
reservation and stable Agent ID `s3xs3-certificate-quadratic`.

The first pass uses the supplied `exact_engine.py` after checking its
ordered left-frame derivative rules against
\(E_a r=-ar\), \(F_a r=ra\), and the stated Koszul convention.  Its
curvature tensor recovers product curvature \(1/3\), mixed zero curvature,
and the required antisymmetries.  For the first representative, the fixed
quadratic Gram expression is additionally compared with the separate full
connection recurrence in `exact_cubic.py`.  This is independent of the
claimed coefficient table, but it is not a wholly independent implementation
of differential geometry; the primary agent is supplying that separate
coordinate-jet check.

## An all-height coefficient test

By diagonal \(SO(4)\)-covariance, the global odd identity can be tested at
\(p=1\), \(q=\cos\theta+i\sin\theta\), with arbitrary orthogonal unit
height vectors \(e,f\in\mathbb R^4\).  The transformed tensors are
\[
 O_{e,f}=\langle e,p\rangle A_P+\langle f,q\rangle A_Q,
\]
\[
 K_{0,e,f}=-\frac15\langle e,p\rangle\langle f,q\rangle
 \begin{pmatrix}0&I\\I&0\end{pmatrix}
 -\frac{787}{5880}\langle f,p\rangle\langle e,q\rangle I_6.
\]
For this torus, the linear correction can be derived without a curvature
engine:
\[
 \mathcal L K_{0,e,f}=
 \frac{-\frac15\langle e,pi\rangle\langle f,qi\rangle
       -\frac{787}{5880}\langle f,p\rangle\langle e,q\rangle}{15/8}.
\]
Thus the full claim is a check of an eight-dimensional quadratic form, with
two four-dimensional self blocks and one cross block.  The cross block is
the raw \(2\mathcal Q(A_P,A_Q)\) block plus the explicit bilinear expression
above.  This does not require sampling orthogonal height pairs, and avoids
assuming that a few positive \(\rho>0\) examples imply the global estimate.

## Point audit results recorded so far

At \(p=1\), \(q=3/5+(4/5)j\), \(u=j\), one has \(\rho=1\) and the
exact computation gives
\[
 \mathcal Q(O)+\mathcal LK_0=\frac{8779429}{1090260000},
\]
equal to the supplied formula.  Its margin over \(1/250\) is
\(4418389/1090260000>0\).

At \(p=(1+i+j+k)/2\), \(q=(-1+7i+7j-k)/10\), \(u=i\), one again has
\(\rho=1\), and
\[
 \mathcal Q(O)=\frac{88951703}{7631820000},\qquad
 \mathcal LK_0=-\frac{389}{220500},
\]
\[
 \mathcal Q(O)+\mathcal LK_0=\frac{528415001}{53422740000},
\]
again exactly equal to the supplied formula with the oriented normal
quarter-turn taking \(pj\) to \(pk\).  This also checks the sign of its
antisymmetric cross term.

At the canonical residual torus point \(p=1\), \(q=3/5+(4/5)i\),
the odd corrected coefficient is zero and
\[
 \mathcal Q(B)=-\frac{896}{21125}
 =\frac{128}{845}\cos(2\theta).
\]
All eight exact LDL pivots of the regular fixed-base curvature Hessian are
positive at the representatives above.  These point checks support, but do
not by themselves prove, the global rational-function identity.

## Current scope boundary

The full rational-function height matrix is assigned to a separate agent;
this agent's duplicate run `CERT-QUAD-002` was deliberately stopped after
that division of labor, with all its processes verified stopped and the
lease released.  This was a canceled duplicate, not a failed identity.

The first pass `CERT-QUAD-001` completed all five rational plus-family
tests and the type-II first-nullity test.  Their evidence is in
`evidence/R20260905-PC-S3XS3-CERT-QUAD-001/`.  The complete-bundle run
`CERT-QUAD-003` has additionally passed the exact full base-and-plane
second-family reduction
\[
 \mathcal Q(H_0)=\frac{282877}{2278125}>\frac1{10}.
\]
Its normal stationarity residual is exactly zero and the fixed quadratic
numerator is \(1/6\).

The same run completed the supplemental `audit_checks.py` assertions:

- the compact displayed \(K_1\) equals the coefficient-table tensor,
  including all first and second ordered frame derivatives, at the
  symbolic relative-angle test supplied by that audit;
- the curvature convention gives product sectional curvature \(1/3\);
- the full connection recurrence, not the sparse cubic expression,
  reproduces the claimed \(B+C\) cubic at
  \(\cos\theta=3/5\), \(\sin\theta=4/5\).

All selected complete-bundle checks passed in approximately 79 seconds.
Evidence is in `evidence/R20260905-PC-S3XS3-CERT-QUAD-003/`; the worker and
its children stopped and the resource lease was released.

No conclusion of globally positive sectional curvature follows from these
checks alone; the remaining mixed, cubic, central, and geometric assembly
obligations must also be verified.

## Why four common-phase checks suffice for the mixed jet

The supplied driver checks the \(O\)-\(B\) transverse derivative at four
common phases.  The following geometric argument justifies its comment
that the expression is linear in the common phase.  An infinitesimal
\(SO(4)\) rotation can be written
\(\dot p=ap+pb\), with fixed imaginary quaternions \(a,b\).  Differentiating
the right and left tangent axes of the rotated circle gives, on \(T_0\),
\[
 \dot u=[i,b],\qquad \dot w=[a,i].
\]
Both are common-phase independent.  The scalar coordinates in \(K_1\)
have weight one, so its differentiated torus restriction has weight one.

There is also a conservative elementary interpolation proof that needs
only an unsimplified degree bound, recorded here to make the finite-phase
verification particularly transparent.

At a residual-torus point, write
\(p=e^{\varphi i}\), \(q=e^{(\varphi+\theta)i}\), keeping \(\theta\)
fixed.  Coordinate heights have common weights \(\pm1\).  The derivatives
of the left and right tangent axes under a fixed ambient \(SO(4)\)
rotation are quadratic expressions in \(p\) and \(pi\), and hence have
weights \(0,\pm2\).  The differentiated torus restriction of \(K_1\)
therefore has weights at most three and is odd under
\(\varphi\mapsto\varphi+\pi\).  The same bound applies to the full
mixed identity (its stated target has weight one).  It lies in
\[
 \operatorname{span}_{\mathbb R}
 \{\cos\varphi,\sin\varphi,\cos3\varphi,\sin3\varphi\}.
\]
Vanishing at \(\varphi=0,\pi/2\) reduces a member of this space to
\[
 A(\cos\varphi-\cos3\varphi)
 +B(\sin\varphi+\sin3\varphi).
\]
At \((\cos\varphi,\sin\varphi)=(3/5,4/5)\) and \((4/5,3/5)\), its
two values are
\[
 \frac1{125}
 \begin{pmatrix}192&144\\144&192\end{pmatrix}
 \begin{pmatrix}A\\B\end{pmatrix}.
\]
The determinant is positive.  Hence the four exact phase checks are
unisolvent and certify the identity at every common phase, for every
relative angle covered by the rational-function computation.  This
argument supplies a rigorous bridge from the finite common-phase checks
to the continuous phase variable without numerical interpolation.
