#!/usr/bin/env python3
"""The entire centre of the fixed 83-dimensional quadratic candidate Q.

Read the existing exact multiplication table; do not rerun its construction.
All coefficients are rational Fractions. Table entries are chronological
products, which reverse conventional path multiplication; centres are unchanged.
This is a calculation about Q, not about an unidentified VOA endomorphism ring.
"""
from fractions import Fraction as F
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import platform
import sys

BASE = Path(__file__).resolve().parent
INPUT = BASE / 'finite_multiplication.json'
D = json.loads(INPUT.read_text())
B = D['basis']
assert len(B) == 83
assert hashlib.sha256((BASE/'finite_multiplication.py').read_bytes()).hexdigest() == D['source_sha256']
T = {(m['first'], m['second']): {int(k): F(v) for k, v in m['terms'].items()}
     for m in D['multiplication']}
names = {b['name']: i for i, b in enumerate(B)}

def add(*terms):
    z = defaultdict(F)
    for c, v in terms:
        for i, a in v.items():
            z[i] += F(c)*a
    return {i: a for i, a in z.items() if a}

def mul(x, y):
    z = defaultdict(F)
    for i, a in x.items():
        for j, b in y.items():
            for k, c in T.get((i, j), {}).items():
                z[k] += a*b*c
    return {k: c for k, c in z.items() if c}

def rref(rows):
    piv = {}
    for row in rows:
        r = {j: F(c) for j, c in row.items() if c}
        while r:
            j = min(r)
            if j in piv:
                c = r[j]
                r = add((1, r), (-c, piv[j]))
            else:
                c = r[j]
                piv[j] = {k: v/c for k, v in r.items()}
                break
    return piv

def nullspace(rows, n):
    piv = rref(rows)
    ans = []
    for t in range(n):
        if t in piv:
            continue
        v = {t: F(1)}
        for j in sorted(piv, reverse=True):
            c = -sum((a*v.get(k, F(0)) for k, a in piv[j].items() if k != j), F(0))
            if c:
                v[j] = c
        ans.append(v)
    return ans

def columns_kernel(cols):
    rows = defaultdict(dict)
    for j, col in enumerate(cols):
        for i, c in col.items():
            rows[i][j] = c
    return nullspace(rows.values(), len(cols))

def v(name):
    return {names[name+':0']: F(1)}

generators = [i for i, b in enumerate(B) if b['degree'] <= 1]
comm_rows = []
for j in generators:
    rows = defaultdict(dict)
    for i in range(len(B)):
        c = add((1, T.get((i, j), {})), (-1, T.get((j, i), {})))
        for k, a in c.items():
            rows[k][i] = a
    comm_rows.extend(rows.values())
centre_basis = nullspace(comm_rows, 83)
assert len(centre_basis) == 9
centre_degrees = {}
for d in range(5):
    ids = [i for i, b in enumerate(B) if b['degree'] == d]
    rows = [{j: r.get(i, F(0)) for j, i in enumerate(ids) if r.get(i)} for r in comm_rows]
    centre_degrees[str(d)] = len(nullspace(rows, len(ids)))
assert centre_degrees == {'0': 1, '1': 0, '2': 3, '3': 0, '4': 5}

def N(a, b, c):
    s = a+b+c
    return add((F(a+c, 2), v('t')), (b, v('x_A')), (-s, v('y_A')),
               (c, v('x_B')), (a, v('y_B')), (a, v('x_C')),
               (-s, v('y_C')), (c, v('x_D')), (b, v('y_D')))

X, Y, Z = N(1, 0, 0), N(0, 1, 0), N(0, 0, 1)
for n in (X, Y, Z):
    assert all(mul(n, {j: F(1)}) == mul({j: F(1)}, n) for j in range(83))

mon_names = ['X^2', 'Y^2', 'Z^2', 'XY', 'XZ', 'YZ']
mon = [mul(X, X), mul(Y, Y), mul(Z, Z), mul(X, Y), mul(X, Z), mul(Y, Z)]
expected_kernel = {0: F(1), 1: F(1), 2: F(-1), 3: F(-2)}
kernel = columns_kernel(mon)
assert len(kernel) == 1
assert not add(*[(c, mon[j]) for j, c in expected_kernel.items()])
socles = [names['s_'+s+':0'] for s in 'TABCD']
matrix = [[m.get(i, F(0)) for m in mon] for i in socles]
assert len(rref([{j: c for j, c in enumerate(row) if c} for row in matrix])) == 5

u, vv, w = add((1, X), (-1, Y), (-1, Z)), add((1, X), (-1, Y), (1, Z)), Y
assert not mul(u, vv)
quad = [mul(u, u), mul(u, w), mul(vv, vv), mul(vv, w), mul(w, w)]
assert len(columns_kernel(quad)) == 0
unit = add(*[(1, v('e_'+s)) for s in 'TABCD'])
cb = [unit, u, vv, w]+quad
assert len(columns_kernel(cb)) == 0
for a in (u, vv, w):
    for b in (u, vv, w):
        for c in (u, vv, w):
            assert not mul(mul(a, b), c)

# Check all 81 basis products against the monomial presentation uv=0, m^3=0.
exponents = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
             (2, 0, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2)]
for i, a in enumerate(exponents):
    for j, b in enumerate(exponents):
        e = tuple(x+y for x, y in zip(a, b))
        expected = {} if sum(e) >= 3 or (e[0] and e[1]) else cb[exponents.index(e)]
        assert mul(cb[i], cb[j]) == expected

socle_rows = []
for q in (u, vv, w):
    rows = defaultdict(dict)
    for j, c in enumerate(cb):
        for k, a in mul(c, q).items():
            rows[k][j] = a
    socle_rows.extend(rows.values())
assert len(nullspace(socle_rows, 9)) == 5

def rank_mul(n):
    cols = [mul(n, c) for c in cb]
    return 9-len(columns_kernel(cols))

sample_ranks = {}
for label, n in [('u', u), ('v', vv), ('w', w), ('u+v+w', add((1, u), (1, vv), (1, w))), ('u^2', quad[0])]:
    sample_ranks[label] = [rank_mul(n), rank_mul(mul(n, n))]
assert sample_ranks == {'u': [3, 1], 'v': [3, 1], 'w': [4, 1], 'u+v+w': [4, 1], 'u^2': [1, 0]}

report = {
    'scope': 'Full centre of the fixed explicit quadratic candidate Q; no VOA comparison',
    'input': str(INPUT.name),
    'input_sha256': hashlib.sha256(INPUT.read_bytes()).hexdigest(),
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'program_bridge': 'Read the certified complete Q basis/table; commutation with all idempotents and arrows is necessary and sufficient for centrality',
    'centre_dimension': len(centre_basis),
    'centre_dimension_by_path_degree': centre_degrees,
    'commutator_constraint_rank': len(rref(comm_rows)),
    'degree_two_centre_generators': ['N(1,0,0)', 'N(0,1,0)', 'N(0,0,1)'],
    'product_matrix_rows': list('TABCD'),
    'product_matrix_columns': mon_names,
    'product_matrix': [[str(c) for c in row] for row in matrix],
    'quadratic_kernel': 'X^2+Y^2-Z^2-2XY',
    'presentation': 'C[u,v,w]/(uv,(u,v,w)^3), u=X-Y-Z, v=X-Y+Z, w=Y',
    'basis_product_checks': 81,
    'cubic_generator_checks': 27,
    'centrality_checks_all_basis': 3*83,
    'socle_dimension': 5,
    'sample_ranks_m_n_m_n_squared': sample_ranks,
    'general_rank_boundary': 'All-parameter rank and zero-product claims are proved algebraically in center_algebra.tex, not inferred from these samples',
    'arithmetic': 'Exact Fraction over Q; one serial process; no added dependencies',
    'runtime': {'python': sys.version, 'executable': sys.executable, 'platform': platform.platform()},
}
print(json.dumps(report, indent=2))
