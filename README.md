# Computational supplement: Projective covers and block reconstruction in a Virasoro tensor category

Version v14, 4 November 2026. Author: Yueqiang Zhao.

This package contains the exact-arithmetic programs and finite data accompanying the manuscript. The full mathematical arguments and all-index recoupling proof are in the manuscript PDF. No installation of third-party Python packages is required.

## Reproduce the calculations

From this package directory, using Python 3 without the `-O` option:

```sh
python3 code/replay.py --output replay_fresh
```

Choose an absent or empty output directory. The driver copies the programs into that directory, regenerates the multiplication data, runs eight checks, and saves stdout, stderr, return codes and source hashes. It leaves the packaged input files unchanged. The recorded run used Python 3.12.14; its environment and exact results are included under `verification/replay_v14/`. Absolute paths in that original run identify its execution environment and are not paths that the recipient must create.

`source/constants.tex` is machine-readable input for comparison against the printed nonunit constants. It is an excerpt of data from the manuscript, not a standalone paper or a LaTeX driver.

## Recorded exact results

| Quantity | Result |
|---|---|
| Homogeneous quotient dimensions, degrees 0 through 5 | 5, 20, 33, 20, 5, 0 |
| Independent scalar quadratic relations | 53 |
| Basis vectors | 83 |
| Directly reduced ordered basis products | 6889 |
| Nonzero basis products | 480 |
| Nonzero multiplet constants | 270 = 105 unit + 165 nonunit |
| Printed nonunit constants compared with regenerated data | 165, all equal |
| Centre dimension; ordered basis products; socle dimension | 9; 81; 5 |
| Bounded recoupling diagnostics | 623 channels, 3805 input monomials, indices at most 10 |

The finite quotient calculation is exhaustive: once degree five vanishes, all higher degrees vanish. `verify_quadratic_table.py` independently encodes the quadratic relations and checks the 83-vector basis and all ordered products, without importing the potential-based construction modules. See `code/README.md` for the role of each program.

The bounded recoupling tests do not replace the general proof. The formal projection identities do not establish that the specified matrices realize the VOA action. These calculations do not certify the external triplet-cover results, categorical descent, or every proof in the manuscript. Generative-AI assistance is disclosed in the manuscript's AI Declaration; execution of these programs is not formal proof-assistant verification or independent human peer review.

`manifest.json` lists every packaged file other than itself with its SHA-256 digest. All original Python programs and the packaged finite multiplication data are unchanged from the supplied v13.1 source. The recorded run was performed afresh for this revision.
