# Accepted Replay Results

These directories preserve outputs from the internal September 5, 2026
review. They are evidence of completed runs, not inputs trusted by a fresh
verification. Running `verification/verify.py` creates a separate timestamped
directory under `replay-results/`.

| Run | Role |
| --- | --- |
| `ENG-001` | Exact rational-point replay. |
| `ENG-002` | Exact generic-angle contracted-curvature calculation over `Q(z)`. |
| `ENG-003` | Original exact type-I verifier replay. |
| `ENG-005` | Fail-closed repeat of the generic exact calculation and stationarity equations. |
| `QUAD-001` | Preliminary exact quadratic audit. |
| `QUAD-003` | Full type-II and algebraic-`K1` exact audits. |
| `GEO-001` | Smooth diagonal path and full-normal diagonal audit. |
| `ROOT-002` | Separately differentiated floating-point coordinate checks. |
| `ROOT-003` | Additional floating-point central diagnostics. |

The precise scope and shared dependencies are recorded in `CLAIM-MAP.md`,
`verification/verification-report.json`, and the three files in `docs/`.
