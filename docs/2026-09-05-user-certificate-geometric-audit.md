# Geometric audit of the September 5 user certificate

Target: EXT-POSCURV-S3XS3.

Source under review: references/user-certificate-2026-09-05/PROOF.md.
This note audits the geometric implications of the displayed formulas. The
exact quadratic and cubic identities are separate algebraic obligations.

## 1. Background and clean zero locus

The background cometric has the stated deformation order. Adding the
right-diagonal Cheeger term with parameter \(1/3\) to the cometric of
\(3I_6\) gives

\[
\frac13\begin{pmatrix}2I&I\\I&2I\end{pmatrix},
\]

whose inverse is \(\left(\begin{smallmatrix}2I&-I\\-I&2I\end{smallmatrix}\right)\).
The left-diagonal deformation adds
\(\frac3{10}\left(\begin{smallmatrix}I&R\\R^T&I\end{smallmatrix}\right)\)
to the cometric, giving exactly the displayed \(g_0^{-1}\).
Since \(10I+9R\) has singular values at most \(19\), the cometric
eigenvalues lie in \([1/3,8/5]\).

For the clean-Hessian assertion, the product curvature Hessian has rank four
at a mixed plane: the independent components are the first-factor component
of the second spanning vector and the second-factor component of the first.
On the four-dimensional family of product-flat planes, the first Cheeger
bracket has rank two and zero condition \([u]=[w]\). Thus the intermediate
zero locus has codimension six and positive normal Hessian.

On this locus the next bracket is a nonzero multiple of \(u\times Ru\).
At \(r=(a,bu)\), differentiation in \(v\mapsto v+t w\), \(w\perp u\), gives

\[
D(u\times Ru)[w]=2(aI+bJ_u)w.
\]

This is invertible on \(u^\perp\), since \(a^2+b^2=1\). At a minus zero
\(r=(0,v)\), \(v\perp u\), varying \(r_0\) and \(v\cdot u\) gives the
independent vectors \(2v\) and \(2u\times v\). These give two additional
normal directions, hence total rank eight. The smooth Cheeger plane
reparametrizations preserve rank. At the zero conditions, they preserve the
actual plane \(P_u\), since the added action field is a multiple of
\((u,\pm u)\).

The plus and minus families are embedded and disjoint: the point and plane
recover \(L=\operatorname{span}(p,pu)\). In particular there is no unaccounted
identification at \(p=\pm q\). The same rank-two section proves smoothness
there.

The fixed-set proof of total geodesicity is valid. The map equal to \(+I\)
on \(L\) and \(-I\) on \(L^\perp\) has determinant \(+1\); its simultaneous
action fixes a plus torus. Composing with the second-factor antipodal map
fixes a minus torus. Both are background isometries. Restriction to the
\(u\)-axis invariant block gives the two stated constant induced metric
matrices and determinants.

## 2. First-nullity and the \(H_0\) assertion

On a surviving torus the left axes are constant \(u\). The right axes are
constant \(w,w\) on a plus torus or \(w,-w\) on a minus torus. Hence \(B,C\)
restrict to constant tensors. The restriction of \(O\) is

\[
h_{ss}=p_0(s),\qquad h_{tt}=q_1(t),\qquad
h_{st}=-\frac27\bigl(p_0(s)+q_1(t)\bigr),
\]

which is first-null by the linearized Gaussian-curvature formula.
The restriction of \(H_0\) is zero: on plus tori \(P(v)u=0\), and on minus
tori \(r_0=0\).

Every plus torus retains its induced metric and total geodesicity along
\(g_0+\tau H_0\). Thus its full tangential connection variation vanishes.
Codazzi gives zero first plane-curvature gradient. At a regular plus point,
fixed-base plane variables form a full normal complement, so the full first
normal gradient vanishes; continuity covers the center. This establishes
the polarized radical identity \(\mathcal Q(H_0,h)=0\) on the plus locus.

There is no symmetry obstruction to opening the minus family with \(H_0\).
It is diagonal \(SO(4)\)-invariant but odd under either individual antipodal
map. The latter symmetry is essential to the background fixed-set proof for
minus tori and is lost under this perturbation.

The fixed quadratic numerator \(1/6\) at a minus zero has an independent
calculation. At \(p=1,q=j,u=i\),

\[
H_0=0,\qquad dH_0=dr_0\otimes
\begin{pmatrix}0&P(j)\\P(j)&0\end{pmatrix},\qquad
dr_0=\alpha_2-\beta_2.
\]

Therefore

\[
D_{H_0}(E_i,E_i)=D_{H_0}(F_i,F_i)=0,\qquad
D_{H_0}(E_i,F_i)=-\tfrac12\operatorname{grad}_{g_0}r_0.
\]

The cometric gives \(|dr_0|^2=2/3\), proving the numerator \(1/6\).
The full Schur subtraction is a separate algebraic obligation.
One minus representative suffices: the ordered orthonormal frame
\((p,pu,q,qu)\) always has ambient orientation \(-1\), so diagonal \(SO(4)\)
is transitive on these frames. There is no omitted second orientation orbit.

## 3. Independent \(O,K_0\) check on the residual torus

Put \(p=e^{si},q=e^{ti}\), \(P=\cos s\), \(Q=\sin t\), \(c=-2/7\). The
inverse induced background metric is

\[
G_T^{-1}=\frac1{30}\begin{pmatrix}29&19\\19&29\end{pmatrix}.
\]

The lower connection variations are

\[
D_{O,ss}^{\flat}=P'(1/2,c),\qquad
D_{O,tt}^{\flat}=Q'(c,1/2),\qquad D_{O,st}^{\flat}=0.
\]

The pairing \((1/2,c)G_T^{-1}(c,1/2)^T\) equals \(-389/5880\).
Thus the second curvature numerator is
\(-389\sin s\cos t/5880\). The off-diagonal and diagonal parts of \(K_0\)
contribute respectively \(\sin s\cos t/5\) and
\(-787\sin s\cos t/5880\), cancelling it exactly.
This verifies \(\mathcal Q(O)+\mathcal LK_0=0\) on \(T_0\), but does not
verify the positive transverse margin away from it.

## 4. Normal-minimum gauge and pole removal

Let
\[
F_\tau=F_0+\tau F_1+\tau^2F_2+\tau^3F_3+O(\tau^4)
\]
have a clean minimum manifold \(Z\), positive normal Hessian, and
\(F_1|_Z=0\). A tubular choice supplies first optimizer \(n(z)\) and reduced
coefficients \(q,r\). Another first optimizer differs by a tangent vector
\(t\in T_zZ\). Its footpoint is \(z_\tau=z+\tau t+O(\tau^2)\), so the
reduced value along it is

\[
\tau^2q(z)+\tau^3\bigl(r(z)+dq_z[t]\bigr)+O(\tau^4).
\]

Changing the second displacement contributes only \(O(\tau^4)\), since the
first displacement makes the normal gradient vanish through order one.
Thus \(q\) is intrinsic, and \(r(z)\) is intrinsic where \(dq_z=0\).
This concerns normalized sectional curvature: numerator calculations must
retain the Gram factor unless their lower coefficients have vanished.

The stated pole removal works at every common phase. Hold \(p=e^{si}\)
fixed, vary \(u=i+t k\), and set
\(q=p(\cos\theta+\sin\theta\,u)\). Then

\[
q^{-1}\dot q=\sin^2\theta\,j+\sin\theta\cos\theta\,k,
\]

and the plane derivatives are \((E_k,F_k)\). Subtracting
\(32\cot\theta/65\) times this tangent vector removes both cotangent
entries and gives exactly \(w_\theta,U_\theta,V_\theta\) in the proof.
No common-phase-dependent adjoint factor is missing.
The resulting smooth expression and the gauge lemma extend the generic
cubic to \(\theta=0,\pi\).

Only \(B\) contributes to the first optimizer on \(T_0\). When its amplitude
is zero, the remaining formal path preserves the reflection fixing \(T_0\),
which stays totally geodesic. Its normalized sectional curvature is
stationary under all fixed-base plane variations, so the first normal
source vanishes at regular points and by continuity at the center.

## 5. Cubic mean and parity

Assume the stated corrected quadratic value and first-jet identities.
The reduced cubic is homogeneous of degree three in
\((\lambda,\delta,\mu,c)\), because \(h\) is linear and \(k\) quadratic in
these amplitudes. Reflection by \(i\) fixes \(T_0\) pointwise and sends
\(\delta\) to \(-\delta\), so the cubic is even in \(\delta\).

When \(\delta=0\), total geodesicity identifies the reduced cubic with the
intrinsic cubic. The intrinsic coefficients of orders zero, one, and two
vanish. Gauss--Bonnet then makes its mean zero with respect to the original
constant area form. This removes every \(\delta\)-free monomial from the
mean, not merely \(\mu^3\).

The remaining monomials are \(\delta^2\lambda,\delta^2\mu,\delta^2c\).
The simultaneous antipodal map reverses \(O,K_1\), preserves all other
relevant tensors and torus area, and kills the first mean. To isolate the
last term set \(\lambda=0\); then the first-factor antipodal map reverses
\(H_0\) and preserves \(B,C,K_2,K_3\), killing its mean. Hence the mean in
the \(B+C,K_2+K_3\) calculation is exactly the coefficient needed for the
full family. If the value \(512/375\) is verified, the full mean is
\((512/375)\varepsilon^3\); no mixed-amplitude mean has been omitted.

An equivariant tubular choice is unnecessary: the gauge lemma gives cubic
invariance under changes of choice on \(T_0\).

## 6. Constants, extension, and uniform positivity

If the corrected cross-functions vanish with their first jets on \(T_0\),
Taylor's theorem and \(\rho\asymp d(T_0,\cdot)^2\) give finite \(C_*\).
Compactness away from a tube supplies the rest of the bound. The choices of
\(\varepsilon,c\) then give the stated plus margin and minus positivity.
There is no circularity: the constants use only known smooth jets at \(g_0\).

The torus area density is constant, so the coordinate Fourier mean and the
normalized Riemannian mean agree. The Poisson denominator satisfies

\[
29(a^2+b^2)+38ab=10(a^2+b^2)+19(a+b)^2>0
\]

away from the zero character. A smooth right-hand side gives a smooth
mean-zero solution. The extension is smooth since its cutoff vanishes near
either undefined normalized projection. The conformal variation is
\(-\tfrac12\Delta_Tf\), with the stated sign.

A compact clean zero locus permits uniform normal minimization and a
uniform fourth-order remainder. If \(q\ge0\), \(q^{-1}(0)=T_0\), and the
corrected cubic is positive there, it stays positive on a fixed neighborhood.
There the cubic dominates the fourth-order error. On the compact remainder
of \(Z\), a positive quadratic minimum dominates all higher terms. Uniform
strict convexity in the normal fibers gives positivity on the tube in the
whole plane bundle; openness preserves positivity on its compact complement.
All tensors are fixed before the parameter limit, so small constants or
large coefficients do not invalidate the argument.

## 7. Exact central replay

Run R20260905-PC-S3XS3-CERT-GEO-001 completed on the Windows worker in
148.34 seconds, using one core and a 1 GiB reservation. Its entire process
tree exited and its lease was released.

The new calculation uses \(p=1\) and the smooth Cayley path

\[
q(t)=\frac{1-a^2t^2}{1+a^2t^2}
     -\frac{2at}{1+a^2t^2}k,\qquad a=\frac{16}{65},
\]

with spanning vectors

\[
X(t)=E_i+\frac{t}{65}(-96E_j-64F_j),\qquad
Y(t)=F_i+\frac{t}{65}(64E_j+96F_j).
\]

Full spatial jets of the metric and tensors were composed with this moving
base before their coefficients through order three were passed to the full
Levi-Civita curvature recurrence. The resulting exact sectional-curvature
coefficients were

\[
(0,0,0,115712/63375).
\]

All eight quadratic graph-stationarity derivatives were exactly zero.
This check reuses the certificate's curvature recurrence, but neither its
singular fixed-base optimizer nor its diagonal sparse contraction.
Both supplied full-normal central routines were then replayed: the diagonal
optimizer/cubic check passed, and the odd-height diagonal quadratic check
returned \(3158117/198450000\), as claimed.

Evidence is in evidence/R20260905-PC-S3XS3-CERT-GEO-001/result.json and
stdout.log. The new reproducible driver is
scripts/audit_user_certificate_center.py.

## Geometric verdict

No missing geometric implication has been found in this assembly. This is
conditional on correct algebraic certificates for the global quadratic
identity, transverse corrections, full minus Schur value, and cubic
identity. The central pole removal and central cubic now additionally have
the exact pole-free replay just recorded. This note does not by itself
replace the separate algebraic audits of the remaining global identities.
