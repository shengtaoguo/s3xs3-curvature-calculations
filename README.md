# Curvature calculations for S3 x S3

Computational companion to *A Metric with Positive Sectional Curvature on
S3 x S3*, September 5, 2026. The version cited by the paper is
`v0.1.0-paper-2026-09-05`.

## Run

Tested with Python 3.11.9, NumPy 2.4.1, and SymPy 1.14.0.
From the repository root:

```sh
python3 -m pip install -r verification/requirements.txt
python3 verification/verify.py
```

The command runs the exact checks and the additional numerical checks
sequentially. It prints `PASS:` for completed checks, ends with
`VERIFIED:`, and exits with status zero on success. Failed assertions or
numerical comparisons produce a nonzero exit status. Do not use Python
`-O` or `-OO`, which disable assertions.

Results and logs are generated in `replay-results/`; they are not committed.
For a shorter run, choose a suite:

```sh
python3 verification/verify.py --suite certificate  # original exact checks
python3 verification/verify.py --suite independent  # additional exact checks
python3 verification/verify.py --suite numerical    # coordinate diagnostics
python3 verification/verify.py --suite smoke        # quick point checks
```

## Files and paper references

`verification/original/` contains the original exact curvature engine,
its verifier, and the two supplemental checks.
`verification/audits/` contains the additional exact checks and the
coordinate implementation. Each file is needed by a listed check or its
imports; there are no packaging or release scripts.

| Paper calculation | Code, relative to `verification/` |
| --- | --- |
| Proposition 4.2: type-II coefficient `282877/2278125` | `original/verify.py:second_family` |
| Proposition 5.1: quadratic matrices and lower bound | `original/verify.py:first_family` |
| Proposition 5.4: mixed transverse derivatives | `original/verify.py:first_family`, `original/audit_checks.py` |
| Proposition 6.1: generic cubic and stationarity equations | `original/verify.py:first_family`, `audits/audit_user_certificate_engine.py --symbolic` |
| Appendix A.4: full normal calculation at the diagonal | `original/verify.py:diagonal`, `original/audit_diagonal_odd.py` |
| Lemma 6.2: smooth diagonal path | `audits/audit_user_certificate_center.py` |
| Additional quadratic and coordinate checks | `audits/audit_user_certificate_quadratic.py`, `audits/audit_user_certificate_coordinates.py` |

Use the top-level verification command above; it sets the required import
paths and runs the assertions.

## Scope

The exact routines check the displayed finite algebraic identities over
the rationals or `Q(z)`. The paper supplies the geometric arguments for
the zero locus, normal Hessian, smooth corrections, averaging, and uniform
positivity, including why the finite common-phase checks suffice.

The contracted-curvature audit shares the original metric and tensor
two-jets. The separate coordinate implementation uses floating-point
arithmetic and provides corroborating checks. The integral mixed correction
is proved in the paper; the corresponding code checks the algebraic
implementation. These programs are not an end-to-end formalization of the
geometric theorem.
