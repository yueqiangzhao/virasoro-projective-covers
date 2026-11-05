# Reproducible exact calculations

Use Python 3 with its standard library; no third-party Python packages are needed. Do not use Python's `-O` option, because the checkers use assertions.

From the package root:

```sh
python3 code/replay.py --output replay_fresh
```

The output directory must be absent or empty. The driver copies the code into `replay_fresh/working`, regenerates the multiplication table there, and runs eight checks sequentially. It never overwrites the source table. Each program's stdout, stderr, elapsed time and return status are retained, with source SHA-256 values and the actual Python version.

## Checks and their scope

- `candidate_certify.py` and `mixed_explicit_relations.py`: the homogeneous quotient through degree five, the 53 independent scalar relations, endpoint dimensions, and equivariance checks. Degree-five vanishing implies vanishing in all higher degrees.
- `finite_multiplication.py`: the 83-vector basis, 480 nonzero basis products, 270 nonzero multiplet constants, associativity and equivariance.
- `verify_quadratic_table.py`: a separate direct encoding of the printed quadratic relations. It imports none of the potential-based modules; it proves basis independence and checks all 6889 ordered basis products by exact homogeneous ideal reduction.
- `check_constants.py`: compares all 165 printed nonunit constants in `source/constants.tex` against the freshly generated table (not a stale packaged copy). The remaining 105 constants are unit channels.
- `center_algebra.py`: the complete 9-dimensional centre, its 81 ordered basis products, and its 5-dimensional socle.
- `actual_projection_check.py`: formal identities in the specified matrix-block pattern, not a proof that the pattern is the VOA action.
- `transport_check.py`: bounded diagnostics, 623 channels and 3805 monomials with indices at most 10. The unbounded formula is proved in the manuscript.

For a standalone comparison after generating a different table:

```sh
python3 code/verify_quadratic_table.py --data replay_fresh/working/finite_multiplication.json
python3 code/check_constants.py --data replay_fresh/working/finite_multiplication.json
```

All arithmetic is exact (`int` and `fractions.Fraction`). These tests do not certify the external triplet-cover theorems, the general categorical descent, the actual VOA comparison, or the arbitrary-index recoupling proof. They are not formal verification or independent human peer review.
