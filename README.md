# Curvature Calculations for S³ × S³

This repository contains the computations accompanying *A Metric with Positive
Sectional Curvature on S³ × S³*. The construction perturbs a nonnegatively
curved metric. The scripts verify its quadratic and cubic curvature identities,
including the correction obtained by minimizing over nearby base points and
two-planes.

## Run

Use Python 3.11 and run from the repository root:

```sh
python3 -m pip install -r verification/requirements.txt
python3 verification/verify.py
```

Success prints `VERIFIED:`; failures exit with a nonzero status.
Results are written to `replay-results/`.

## Files

`original/` contains the main exact calculations used in the paper.
`audits/` contains supplementary recomputations of selected coefficients and
numerical checks in sphere coordinates. The entry point runs both groups.

```text
verification/
├── verify.py                                    # runs the checks
├── requirements.txt                             # Python dependencies
├── original/                                    # main exact calculations
│   ├── verify.py                                # main exact assertions
│   ├── exact_engine.py                          # metric, tensor jets, and curvature
│   ├── exact_sparse.py                          # linear variation and contractions
│   ├── exact_cubic.py                           # third-order curvature expansion
│   ├── exact_cross_checks.py                    # mixed terms
│   ├── exact_diagonal.py                        # full normal calculation at the diagonal
│   ├── exact_h0.py                              # type-II perturbation
│   ├── audit_checks.py                          # mixed correction and cubic checks
│   └── audit_diagonal_odd.py                    # diagonal height-tensor check
└── audits/                                      # supplementary checks
    ├── audit_user_certificate_engine.py         # contracted-curvature check
    ├── audit_user_certificate_center.py         # smooth diagonal path
    ├── audit_user_certificate_quadratic.py      # additional quadratic checks
    ├── audit_user_certificate_coordinates.py    # numerical coordinate checks
    └── coordinate_geometry.py                   # coordinate derivatives and curvature
```

## Paper Claims and Scripts

Links below open the relevant verification code.

| Paper result | Main calculation | Additional check |
| --- | --- | --- |
| Proposition 4.2 | [Type-II coefficient](verification/original/verify.py#L122) | — |
| Proposition 5.1 | [Quadratic matrices and lower bound](verification/original/verify.py#L23) | [Quadratic check](verification/audits/audit_user_certificate_quadratic.py#L97) |
| Proposition 5.4 | [Mixed transverse derivatives](verification/original/verify.py#L23) | [Mixed-term check](verification/original/audit_checks.py) |
| Proposition 6.1 | [Cubic coefficient and stationarity](verification/original/verify.py#L23) | [Contracted calculation](verification/audits/audit_user_certificate_engine.py#L117) |
| Lemma 6.2 | [Smooth diagonal path](verification/audits/audit_user_certificate_center.py#L48) | — |
| Appendix A.4 | [Full normal calculation at the diagonal](verification/original/verify.py#L143) | [Height-tensor check](verification/original/audit_diagonal_odd.py) |

These calculations use exact arithmetic over the rationals or ℚ(z).
The contracted calculation shares the main metric and tensor two-jets.
[Separate coordinate checks](verification/audits/audit_user_certificate_coordinates.py)
use floating-point arithmetic. The global geometric argument and the integral
mixed correction are proved in the paper.
