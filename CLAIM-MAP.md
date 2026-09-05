# Paper-to-Code Verification Map

This map binds each computer-assisted statement in the 23-page manuscript
snapshot to its executable check. Formula numbers refer to the PDF in
`artifacts/manuscript/paper.pdf`. The `PASS` marker is printed only after the
corresponding assertions have succeeded.

| Paper statement | Exact executable location | Accepted evidence | What is and is not checked |
| --- | --- | --- | --- |
| Proposition 4.2, (4.9), Appendix A.4: type-II coefficient `282877/2278125` | `verification/original/verify.py`, function `second_family`, marker `PASS: full base-and-plane second-family coefficient`; `verification/audits/audit_user_certificate_quadratic.py --kind points` | `artifacts/accepted-results/R20260905-PC-S3XS3-CERT-QUAD-003/` | Includes the base direction absent from a fixed-base calculation and the Schur subtraction. Transitivity and the global zero-locus geometry are proved in the paper. |
| Proposition 5.1, (5.2)--(5.3), (A.5): three coefficient matrices and the `1/250` lower bound | `verification/original/verify.py`, function `first_family`, markers `PASS: exact uniform first-family quadratic positivity margin` and `PASS: background identities and first-family H0 cancellation`; `verification/audits/audit_user_certificate_engine.py --symbolic` | `artifacts/accepted-results/R20260905-PC-S3XS3-CERT-ENG-003/` and `artifacts/accepted-results/R20260905-PC-S3XS3-CERT-ENG-005/` | Checks every displayed rational-function matrix entry and the final rational inequality. The SO(4) reduction from arbitrary tori is mathematical. |
| Proposition 5.4, (5.10)--(5.13), (A.9)--(A.11): mixed transverse jets | `verification/original/verify.py`, function `first_family`, markers `PASS: B-C transverse-jet compatibility` and `PASS: B-O transverse-jet compatibility`; `verification/original/audit_checks.py` | `artifacts/accepted-results/R20260905-PC-S3XS3-CERT-QUAD-003/` | Checks the four transverse generators and exact common phases for the algebraic `K1`. The paper proves that the phases determine the relevant trigonometric function and separately constructs the integral `K1`. |
| Proposition 6.1, (6.1)--(6.3): generic cubic and eight stationarity equations | `verification/original/verify.py`, function `first_family`, marker `PASS: cubic coefficient 512/375 + (9728/21125) cos(2 theta)`; `verification/audits/audit_user_certificate_engine.py --symbolic` | `artifacts/accepted-results/R20260905-PC-S3XS3-CERT-ENG-002/` and `artifacts/accepted-results/R20260905-PC-S3XS3-CERT-ENG-005/` | Exact over `Q(z)`. The contracted calculation does not call either original cubic evaluator, but it shares the original metric and tensor two-jets. |
| Appendix A.4: full normal-space calculation at `p=q=1` | `verification/original/verify.py`, function `diagonal`; `verification/original/audit_diagonal_odd.py` | `artifacts/accepted-results/R20260905-PC-S3XS3-CERT-GEO-001/` and `artifacts/accepted-results/R20260905-PC-S3XS3-CERT-QUAD-003/` | Checks the ten-variable stationarity system, rank-eight normal Hessian, and central quadratic identities. The identification of the displayed kernel with tangent directions is explained geometrically in the paper. |
| Lemma 6.2 and Appendix A.4: smooth diagonal curve and cubic `115712/63375` | `verification/audits/audit_user_certificate_center.py` | `artifacts/accepted-results/R20260905-PC-S3XS3-CERT-GEO-001/` | Reuses the curvature recurrence but not the singular fixed-base optimizer or the original diagonal contraction. It composes the spatial jets with a smooth moving-base curve. |
| Cross-check of selected coefficients in independent sphere coordinates | `verification/audits/audit_user_certificate_coordinates.py` | `artifacts/accepted-results/R20260905-PC-S3XS3-CERT-ROOT-002/` and `artifacts/accepted-results/R20260905-PC-S3XS3-CERT-ROOT-003/` | Analytically differentiates separate sphere-coordinate formulas and imports no original curvature routine. It is floating-point corroboration, not an exact certificate. |

The suite names in `verification/verify.py` select these checks as follows:

- `certificate`: the original full verifier and its two supplemental audits;
- `independent`: rational-point replay, exact contracted cubic, smooth diagonal
  path, and exact quadratic point audit;
- `numerical`: the separately differentiated coordinate diagnostics;
- `all`: all three groups, after the checksum inventory has passed.

The global implication from these identities to positive curvature on every
two-plane is not an executable assertion in this repository. It uses the
normal-Hessian, smoothness, parity, Poisson, and compactness arguments in the
paper.
