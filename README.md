# Curvature Calculations for S³ × S³

## Run

Use Python 3.11 and run from the repository root:

```sh
python3 -m pip install -r verification/requirements.txt
python3 verification/verify.py
```

Success prints `VERIFIED:`; failures exit with a nonzero status.
Results are written to `replay-results/`.

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
