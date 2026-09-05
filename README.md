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

```text
verification/
├── verify.py                              # runs the checks
├── requirements.txt                       # Python dependencies
├── original/
│   ├── verify.py                          # main exact assertions
│   ├── exact_engine.py                    # metric, tensor jets, and curvature
│   ├── exact_sparse.py                    # linear variation and contractions
│   ├── exact_cubic.py                     # third-order curvature expansion
│   ├── exact_cross_checks.py              # mixed terms
│   ├── exact_diagonal.py                  # full normal calculation at the diagonal
│   ├── exact_h0.py                        # type-II perturbation
│   ├── audit_checks.py                    # mixed correction and cubic checks
│   └── audit_diagonal_odd.py              # diagonal height-tensor check
└── audits/
    ├── audit_user_certificate_engine.py   # contracted-curvature check
    ├── audit_user_certificate_center.py   # smooth diagonal path
    ├── audit_user_certificate_quadratic.py # additional quadratic checks
    ├── audit_user_certificate_coordinates.py # numerical coordinate checks
    └── coordinate_geometry.py            # coordinate derivatives and curvature
```

## Paper Claims and Scripts

| Paper claim | Script and function |
| --- | --- |
| Proposition 4.2: type-II coefficient `282877/2278125` | [verify.py](verification/original/verify.py), `second_family` |
| Proposition 5.1: quadratic matrices and lower bound | [verify.py](verification/original/verify.py), `first_family`; [quadratic audit](verification/audits/audit_user_certificate_quadratic.py) |
| Proposition 5.4: mixed transverse derivatives | [verify.py](verification/original/verify.py), `first_family`; [mixed checks](verification/original/audit_checks.py) |
| Proposition 6.1: generic cubic and stationarity equations | [verify.py](verification/original/verify.py), `first_family`; [contracted audit](verification/audits/audit_user_certificate_engine.py), `symbolic_run` |
| Appendix A.4: full normal calculation at the diagonal | [verify.py](verification/original/verify.py), `diagonal`; [diagonal height check](verification/original/audit_diagonal_odd.py) |
| Lemma 6.2: smooth diagonal path | [diagonal path audit](verification/audits/audit_user_certificate_center.py) |
| Additional numerical checks | [coordinate audit](verification/audits/audit_user_certificate_coordinates.py) |

The exact scripts verify finite identities over the rationals or `Q(z)`.
The contracted audit shares the original metric and tensor two-jets;
coordinate checks use floating-point arithmetic. The global geometric
argument and the integral mixed correction are proved in the paper.
