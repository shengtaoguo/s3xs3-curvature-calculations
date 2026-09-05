# Certificate for Positive Sectional Curvature on S3 x S3

This repository is the versioned computational companion to
*A Metric with Positive Sectional Curvature on S3 x S3*.
It contains the exact symbolic calculations used for the finite curvature
identities in the paper, accepted replay outputs, and additional audits of
the most sensitive normal-minimization calculations.

The repository does **not** formalize or independently prove the whole
geometric theorem. The paper proves the zero-locus geometry, smooth global
extensions, parity and averaging arguments, and the final compactness step.
The programs verify the finite algebraic identities identified in
[CLAIM-MAP.md](CLAIM-MAP.md).

## Paper snapshot

The exact manuscript snapshot checked by this repository is
[artifacts/manuscript/paper.pdf](artifacts/manuscript/paper.pdf). Its title,
page count, SHA-256 digest, and source checkpoint are recorded in
[repository-manifest.json](repository-manifest.json). The release tag
`v0.1.0-paper-2026-09-05` fixes the repository version cited by the paper.

## Verification

The tested environment is Python 3.11.9 with NumPy 2.4.6 and SymPy 1.14.0.
No GPU, network service, API key, server account, TeX installation, or
optimizer is needed. From the repository root, create a virtual environment
if desired, install the pinned dependencies, and run

```sh
python3 -m pip install -r verification/requirements.txt
python3 verification/verify.py
```

The default command runs all exact checks and the separately labelled
floating-point diagnostics. It is single-process and normally takes about
5--15 minutes, depending on the machine. A successful run prints a line
beginning `VERIFIED:` and exits with status `0`. Do not use Python's `-O`
or `-OO` options: the certificate programs deliberately use assertions.

For a quick standard-library-only integrity check, run

```sh
python3 verification/verify_checksums.py
```

Available bounded suites are

```sh
python3 verification/verify.py --suite integrity
python3 verification/verify.py --suite smoke
python3 verification/verify.py --suite certificate
python3 verification/verify.py --suite independent
python3 verification/verify.py --suite numerical
python3 verification/verify.py --suite all
```

Each mathematical run writes new logs beneath `replay-results/`; that
directory is ignored by Git, so accepted evidence is never overwritten.
A failed assertion, nonzero child exit, missing dependency, or checksum
mismatch makes the entry point fail.

## Repository layout

```text
s3xs3-positive-curvature-certificate/
├── artifacts/
│   ├── manuscript/
│   │   └── paper.pdf                 # exact paper snapshot
│   └── accepted-results/             # recorded JSON outputs and logs
├── docs/                              # scope and dependency audits
├── verification/
│   ├── verify.py                      # only mathematical entry point
│   ├── verify_checksums.py            # complete inventory check
│   ├── requirements.txt               # pinned Python dependencies
│   ├── verification-report.json       # accepted-run summary
│   ├── original/                       # original exact certificate engine
│   └── audits/                         # additional exact and numerical audits
├── CLAIM-MAP.md                        # paper claim to executable check
├── repository-manifest.json            # claim, artifact, and version metadata
└── SHA256SUMS                          # hashes of all distributed files
```

Only `verification/verify.py` is intended as the mathematical entry point.
Its exact flow is

```text
verify.py
├── original/verify.py
├── original/{audit_checks.py,audit_diagonal_odd.py}
├── audits/audit_user_certificate_engine.py
├── audits/audit_user_certificate_center.py
├── audits/audit_user_certificate_quadratic.py
└── audits/audit_user_certificate_coordinates.py
    └── audits/{so4_matrix_twojet.py,so4_fd_probe.py}
```

The original exact engine and the contracted-curvature audit share the
metric and tensor two-jets. The coordinate audit differentiates the sphere
charts and tensors separately but uses floating-point arithmetic. These
dependency boundaries are stated in the claim map and the audit notes.

## What the checks establish

The exact suites check, among other identities:

- the type-I quadratic coefficient matrices and the uniform lower bound;
- the type-II full base-and-plane Schur coefficient
  `282877/2278125`;
- the mixed transverse-jet cancellations, including the algebraic `K1`;
- the generic cubic coefficient
  `512/375 + (9728/21125) cos(2 theta)` and all eight stationarity equations;
- the full-normal diagonal calculation and cubic coefficient
  `115712/63375` along a smooth moving-base curve.

The finite common-phase checks are not presented as sampling evidence:
the paper proves why evaluation at those phases is injective on the relevant
finite-dimensional trigonometric space. Likewise, the rational angular
parameterization is supplemented by a direct calculation at the diagonal.

The code does not establish the classification of all zero planes, the
positive normal Hessian globally, the smooth integral construction of the
geometric mixed correction, the parity reduction of the cubic average, the
Poisson extension, or the final uniform positivity argument. Those are
mathematical parts of the paper.

## Provenance and limitations

`verification/original/` preserves the exact scripts supplied with the
construction. `verification/audits/` contains later checks, including a
contracted-curvature calculation that avoids the original cubic evaluator,
a smooth diagonal path, and a separately differentiated coordinate audit.
The accepted outputs are included for inspection; reruns create separate
results and do not silently replace them.

This is a certificate companion, not a search repository. It does not include
the exploratory parameter search, unrelated failed candidates, machine
launchers, credentials, or the larger research workspace. The rational data
are fixed inputs, and no numerical minimum is used as proof of positive
sectional curvature.
