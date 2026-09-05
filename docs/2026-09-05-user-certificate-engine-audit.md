# Exact engine and generic cubic audit of the supplied certificate

Date: 2026-09-05. Target: `EXT-POSCURV-S3XS3`.

Status: the generic-angle cubic identity and its eight graph-stationarity
equations have passed an exact rational-function verification using a separately
derived contracted-curvature formula. This note does **not** certify the global
positive-curvature construction: the singular base-and-plane geometry, the
remaining coefficient identities, and the geometric assembly have separate
proof obligations.

## 1. Scope and inputs

The immutable supplied files are in
`references/user-certificate-2026-09-05/`. The notation \(g_0,B,C,K_2,K_3\)
in this note refers to the explicit formulas in that bundle's `PROOF.md`.

The audit implementation is `scripts/audit_user_certificate_engine.py`; its
native-Windows launcher is `scripts/run_user_certificate_engine_windows.ps1`.
All substantive computations ran externally on `windows`, under stable agent
ID `s3xs3-certificate-engine`, one logical processor and a 1 GiB aggregate
reservation at a time. No numerical computation ran on the Mac.

## 2. Static audit of the frame conventions

The supplied engine uses left-invariant frame fields \(E_i,F_i\), with
\([E_i,E_j]=2\epsilon_{ijk}E_k\) and the corresponding relation on the
second factor. Its derivative array convention is
\[
 (dg)_{aij}=e_a(g_{ij}),\qquad
 (ddg)_{abij}=e_a e_b(g_{ij}).
\]
For \(R=\operatorname{Ad}_{p^{-1}q}\), differentiation gives
\[
 E_iR=-2J_iR,\quad F_iR=2RJ_i,\quad
 E_lE_iR=4J_iJ_lR.
\]
The reversed order in the last matrix product is therefore intentional, not a
commutator bug. The other second derivative blocks in `g_jets` agree with the
same convention.

The lower Christoffel formula, curvature convention
\(R(e_i,e_j)e_k=\nabla_i\nabla_j e_k-\nabla_j\nabla_i e_k
-\nabla_{[e_i,e_j]}e_k\), and lowering of the last index agree with
\(\sec(X\wedge Y)=R(X,Y,Y,X)/|X\wedge Y|^2\). The coefficient dictionary
`Bco` agrees with the factor \(1/2\) in the displayed definition of \(B\)
and the convention \(\lambda\odot\mu=(\lambda\otimes\mu+
\mu\otimes\lambda)/2\).

## 3. Exact rational-point replay

Run `R20260905-PC-S3XS3-CERT-ENG-001` evaluated
\[
 p=1,\qquad q=\frac35+\frac45i.
\]
It used the supplied full connection recurrence, without importing
`exact_sparse`. The exact results were
\[
 \mathcal Q(B)=-\frac{896}{21125},\qquad
 \mathcal Q(B+C)+\mathcal L(K_2+K_3)=0,
\]
and
\[
 \mathcal R(B+C,K_2+K_3)=\frac{1958912}{1584375}
 =\frac{512}{375}-\frac{7}{25}\frac{9728}{21125}.
\]
The four moving-plane coefficients through degree three were exactly
\(0,0,0,1958912/1584375\).

The following additional tests had zero exact residual in every entry:

- the displayed cometric identity;
- the frame-commutator identities for the metric and \(B\) two-jets;
- first-pair and last-pair skew symmetry, pair interchange, and the first
  Bianchi identity for background curvature, the \(B\) linearized curvature,
  and the cubic curvature tensor.

Each curvature symmetry test covered all \(6^4=1296\) entries. The run
completed in 39.1 seconds and its lease was released after process exit.
The detailed evidence is
`evidence/R20260905-PC-S3XS3-CERT-ENG-001/result.json`.

## 4. An independent contracted-curvature identity

The rational-function replay did not call either supplied cubic evaluator.
For frame fields \(X,Y\) with coefficients constant near the evaluated point,
write
\[
 B(X,Y,Z)=g(\nabla_XY,Z),\qquad c=[X,Y].
\]
The defining curvature formula and metric compatibility give
\[
\begin{aligned}
 R(X,Y,Y,X)={}&X\bigl(B(Y,Y,X)\bigr)
             -Y\bigl(B(X,Y,X)\bigr)\\
 &+B(X,Y,\cdot)^Tg^{-1}B(X,Y,\cdot)
  -B(X,X,\cdot)^Tg^{-1}B(Y,Y,\cdot)\\
 &-B(X,Y,c)-B(c,Y,X).
\end{aligned}
\]
For example, the first connection-product term is obtained by combining
\(-X(g)(X,\nabla_YY)\) with \(B(X,\nabla_YY,X)\); metric compatibility
turns their sum into \(-B(X,X,\nabla_YY)\). The second product similarly
becomes \(B(Y,X,\nabla_XY)\), and torsion freeness supplies
\(\nabla_YX=\nabla_XY-c\).

On a background flat torus, \(B(X,X,\cdot)\), \(B(X,Y,\cdot)\), and
\(B(Y,Y,\cdot)\) vanish at order zero. Consequently their Gram terms through
order three require only
\[
 g_\tau^{-1}=g_0^{-1}-\tau g_0^{-1}h g_0^{-1}+O(\tau^2).
\]
The implementation expands the displayed scalar curvature expression directly,
using exact sparse multilinear contractions. This avoids constructing a full
third-order Riemann tensor and provides an algebraically distinct check on the
supplied cubic routines. It still shares the supplied metric and tensor two-jet
definitions; independence of those definitions is a separate issue.

## 5. Exact rational-function result

Run `R20260905-PC-S3XS3-CERT-ENG-002` used
\[
 p=1,\quad q=\cos\theta+i\sin\theta,\qquad
 \cos\theta=\frac{1-z^2}{1+z^2},\quad
 \sin\theta=\frac{2z}{1+z^2}
\]
over \(\mathbb Q(z)\), and the explicit first graph displacement in equation
(26) of the supplied text. Before the symbolic calculation, the contracted
formula independently reproduced the exact rational-point cubic above.

The direct contraction then gave
\[
 [\tau^0]\sec=[\tau^1]\sec=[\tau^2]\sec=0,
\]
and
\[
 [\tau^3]\sec
 =\frac{115712z^4-2048z^2+115712}{63375(1+z^2)^2}
 =\frac{512}{375}+\frac{9728}{21125}\cos(2\theta).
\]
Every equality here is an identity of rational functions, not interpolation
or numerical reconstruction.

For each of the eight graph-coordinate unit vectors \(e_a\), the run also
computed the second coefficient at \(\eta+e_a\) and \(\eta-e_a\). Since this
coefficient is quadratic in the first displacement, their half-difference is
its exact derivative at \(\eta\). All eight derivatives vanished identically
over \(\mathbb Q(z)\). This verifies stationarity, not merely evaluation on
one favorable moving plane.

The run completed in 106.8 seconds and its lease was released. Evidence:
`evidence/R20260905-PC-S3XS3-CERT-ENG-002/result.json`.

This check applies at generic angles where the fixed-base Hessian is
invertible. The rational expression has a smooth finite limit at the central
angles, but that observation alone does not replace the required moving-base
normal calculation there.

## 6. Common-phase sampling in the supplied mixed-jet check

The supplied `verify.first_family()` tests the \(O\)-\(B\) correction at four
exact common phases. Those checks require a finite-dimensional justification.
On the residual torus set
\(p=(\cos\varphi,\sin\varphi,0,0)\), \(q=p r\). In `bo_Lgrad`, the
undifferentiated left and right axes are constant \(i\); their transverse
variations are quadratic in \(\cos\varphi,\sin\varphi\). The height factors
and their variations are linear. Thus the difference of the proposed two
sides lies in
\[
 \operatorname{span}\{\cos\varphi,\sin\varphi,
                    \cos3\varphi,\sin3\varphi\}.
\]
The phases \((1,0),(0,1),(3/5,4/5),(4/5,3/5)\) are unisolvent for this
space. Indeed, vanishing at the first two phases leaves a combination
\(c(\cos3\varphi-\cos\varphi)+d(\sin3\varphi+\sin\varphi)\).
Vanishing at the last two phases respectively gives \(d=4c/3\) and
\(d=3c/4\), forcing \(c=d=0\). Hence exact tests at these four phases can
establish the full common-phase identity once the stated degree bound is
checked. They should not be described as unexplained representative samples.

## 7. Complete supplied first-family replay

Run `R20260905-PC-S3XS3-CERT-ENG-003` completed in 81.6 seconds with every
assertion in `verify.first_family()` passing over \(\mathbb Q(z)\). Its scope
includes:

- all background curvature symmetries and the frame-commutator identity;
- the \(H_0\) connection and first-curvature-gradient cancellation on the
  generic first-family torus;
- the entire two self-coefficient matrices and corrected cross-coefficient
  matrix for \(\mathcal Q(O)+\mathcal L K_0\), and the exact lower bound
  \(\rho/250\);
- the \(B\)-\(B\) coefficient and \(K_2\) cancellation;
- the two \(B\)-\(C\) transverse derivative identities and the \(8/5\)
  factor in \(K_3\);
- the four \(O\)-\(B\) transverse derivative identities and the \(K_1\)
  cancellation, with the common-phase justification in Section 6;
- agreement of the supplied sparse and full linearized-curvature evaluators;
- the generic cubic coefficient and its explicit minimizing displacement.

The terminal output is
`evidence/R20260905-PC-S3XS3-CERT-ENG-003/stdout.log`, and the terminal result is
`evidence/R20260905-PC-S3XS3-CERT-ENG-003/result.json`. Its resource lease was
released after process exit.

The independent runs initially reported residuals explicitly, rather than
raising on every nonzero residual. Every recorded residual was exactly zero.
The current implementation has subsequently been made fail-closed: it asserts
the vanishing of all lower coefficients, the cubic identity, each of the eight
stationarity equations, and all rational-point diagnostics before reporting
completion. The first strict-assert replay, `R20260905-PC-S3XS3-CERT-ENG-004`,
was interrupted by a registry-lock collision in the monitor, not an algebraic
failure. All recorded processes were verified to have exited; the last worker
status was preserved and the lease released. The private launcher now retries
such lock collisions. Its replay is `R20260905-PC-S3XS3-CERT-ENG-005`.
That strict-assert replay completed successfully in 96.8 seconds: all three
lower coefficients vanished, the cubic rational function agreed identically,
and each of the eight graph-stationarity assertions passed. Its lease was
released after all worker processes exited. The final evidence is
`evidence/R20260905-PC-S3XS3-CERT-ENG-005/result.json`.

The frozen audit script has SHA-256
`8901061d497a9a93f1a7ef60301f17317b5d62e59249b162bb5d6310ed83fb76`.
The externally inspected dependency versions were Python 3.11.9, NumPy 2.4.6,
and SymPy 1.14.0 (native 64-bit Windows). The mathematical audit script is
portable; the separate PowerShell resource launchers are machine-specific and
are not required for a collaborator's certificate check.

## 8. Verdict and remaining scope

No algebraic error was found in the supplied first-family identities. The
strongest independently checked statement here is the exact generic-angle
cubic identity together with all eight stationarity equations. The full
first-family coefficient package also reproduces exactly, and the finite
common-phase checks have the unisolvence justification above.

This leaves no identified gap **within this audit's first-family algebraic
scope**. This audit does not by itself certify the clean full normal Hessian,
the singular moving-base calculations, the second zero-plane component, the
global symmetry arguments for the cubic mean, or the final perturbation
assembly. Those are distinct checks required before promoting the construction
to a proof of positive curvature on all planes.
