#!/usr/bin/env python3
"""Verify the printed quadratic presentation directly against the finite table.

This checker uses no Frobenius potential and imports none of the calculation
modules. It independently encodes the 53 scalar relations of the manuscript,
forms every homogeneous prefix-relation-suffix multiple through degree five,
and uses exact sparse Gaussian elimination over Q. It checks independence of
the supplied 83 basis representatives and every one of their 6889 ordered
products, including the zero products. The arrows and proposed answers are
read from JSON; their correctness is tested, not assumed.

This verifies a finite presented algebra, not its identification with a VOA
endomorphism algebra. That comparison remains a theorem of the manuscript.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
import hashlib
import json
from pathlib import Path
import sys

F = Fraction
PathWord = tuple[int, ...]
Vector = dict[int, Fraction]


def reduce_vector(row: Vector, pivots: dict[int, Vector]) -> Vector:
    """Return the unique remainder in the free-column complement."""
    r = {j: F(c) for j, c in row.items() if c}
    for j in sorted(pivots):
        c = r.get(j, F(0))
        if c:
            for k, v in pivots[j].items():
                z = r.get(k, F(0)) - c*v
                if z:
                    r[k] = z
                else:
                    r.pop(k, None)
    return r


def insert(row: Vector, pivots: dict[int, Vector]) -> None:
    r = dict(row)
    while r:
        j = min(r)
        if j not in pivots:
            c = r[j]
            pivots[j] = {k: v/c for k, v in r.items() if v}
            return
        c = r[j]
        for k, v in pivots[j].items():
            z = r.get(k, F(0)) - c*v
            if z:
                r[k] = z
            else:
                r.pop(k, None)


def verify(data_path: Path) -> dict:
    data = json.loads(data_path.read_text())
    vertices = tuple('TABCD')
    arrows = [tuple(x) for x in data['arrows']]
    # Independently prescribe the quiver, but accept the JSON's numbering.
    wanted = set()
    for i, j, weights in [('T','A',(0,)), ('T','B',(0,)),
                          ('A','C',(1,-1)), ('A','D',(1,-1)),
                          ('B','C',(1,-1)), ('B','D',(1,-1))]:
        for s,t in [(i,j),(j,i)]:
            wanted.update((s,t,w) for w in weights)
    assert len(arrows) == 20 and set(arrows) == wanted
    number = {a:i for i,a in enumerate(arrows)}
    outgoing = {s:tuple(i for i,a in enumerate(arrows) if a[0] == s) for s in vertices}

    @lru_cache(None)
    def paths(s: str, d: int) -> tuple[PathWord, ...]:
        if d == 0:
            return ((),)
        return tuple((i,)+tail for i in outgoing[s]
                     for tail in paths(arrows[i][1],d-1))

    def end(s: str, p: PathWord) -> str:
        return arrows[p[-1]][1] if p else s

    def word(vs: str, ws: tuple[int, ...]) -> PathWord:
        return tuple(number[(vs[k],vs[k+1],w)] for k,w in enumerate(ws))

    relations: list[tuple[str,str,dict[PathWord,Fraction]]] = []
    def relation(*terms: tuple[int,str,tuple[int,...]]) -> None:
        r: dict[PathWord,Fraction] = {}
        endpoints = {(v[0],v[-1]) for _,v,_ in terms}
        assert len(endpoints) == 1
        for c,vs,ws in terms:
            p = word(vs,ws)
            r[p] = r.get(p,F(0))+c
        s,t = next(iter(endpoints))
        relations.append((s,t,{p:c for p,c in r.items() if c}))

    # Symmetric returns (three per directed V-edge).
    for i,j in [('A','C'),('A','D'),('B','C'),('B','D'),
                ('C','A'),('C','B'),('D','A'),('D','B')]:
        relation((1,i+j+i,(1,1)))
        relation((1,i+j+i,(-1,-1)))
        relation((1,i+j+i,(1,-1)),(1,i+j+i,(-1,1)))
    relation((1,'TAT',(0,0)),(1,'TBT',(0,0)))
    for i in 'AB':
        relation((1,i+'T'+i,(0,0)),(1,i+'C'+i,(1,-1)),
                 (1,i+'D'+i,(1,-1)))
    for s,m1,m2,t in [('A','D','C','B'),('B','D','C','A'),
                       ('C','B','A','D'),('D','B','A','C')]:
        for u,v in product((1,-1),repeat=2):
            relation((1,s+m1+t,(u,v)),(-1,s+m2+t,(v,u)))
    for s,t in [('A','B'),('B','A')]:
        relation((1,s+'T'+t,(0,0)),(-1,s+'C'+t,(1,-1)),
                 (1,s+'C'+t,(-1,1)))
    for w in (1,-1):
        relation((1,'TAC',(0,w)),(-1,'TBC',(0,w)))
        relation((1,'CAT',(w,0)),(-1,'CBT',(w,0)))
        relation((1,'TAD',(0,w)),(1,'TBD',(0,w)))
        relation((1,'DAT',(w,0)),(1,'DBT',(w,0)))
    assert len(relations) == 53

    # Each reducer represents the full degree-d homogeneous ideal component.
    reducers = {}
    path_counts=[]; quotient_dimensions=[]; relation_ranks=[]
    for d in range(6):
        count = rank = 0
        for s in vertices:
            for t in vertices:
                pp = tuple(p for p in paths(s,d) if end(s,p) == t)
                ids = {p:j for j,p in enumerate(pp)}
                piv: dict[int,Vector] = {}
                if d >= 2:
                    for left_length in range(d-1):
                        for pre in paths(s,left_length):
                            u = end(s,pre)
                            for rs,rt,rel in relations:
                                if rs != u:
                                    continue
                                for post in paths(rt,d-left_length-2):
                                    if end(rt,post) == t:
                                        insert({ids[pre+p+post]:c for p,c in rel.items()},piv)
                reducers[d,s,t] = (ids,piv)
                count += len(pp); rank += len(piv)
        path_counts.append(count); relation_ranks.append(rank)
        quotient_dimensions.append(count-rank)
    assert path_counts == [5,20,86,360,1548,6480]
    assert quotient_dimensions == [5,20,33,20,5,0]
    assert relation_ranks[2] == 53

    basis = data['basis']
    assert len(basis) == 83
    polys=[]; groups=defaultdict(list)
    for i,b in enumerate(basis):
        d,s,t = b['degree'],b['source'],b['target']
        poly = {tuple(r['arrows']):F(r['coefficient']) for r in b['representative']}
        ids,piv= reducers[d,s,t]
        assert all(len(p)==d and p in ids for p in poly)
        reduced = reduce_vector({ids[p]:c for p,c in poly.items()},piv)
        polys.append(poly)
        groups[d,s,t].append(reduced)
    for key,(ids,piv) in reducers.items():
        if key[0] == 5:
            continue
        independent: dict[int,Vector] = {}
        for r in groups[key]:
            insert(r,independent)
        assert len(independent) == len(groups[key]) == len(ids)-len(piv), key

    table={}
    for x in data['multiplication']:
        key=(x['first'],x['second'])
        assert key not in table
        table[key]={int(k):F(v) for k,v in x['terms'].items() if F(v)}
    products_checked=nonzero=0
    for i,a in enumerate(basis):
        for j,b in enumerate(basis):
            expected=table.get((i,j),{})
            products_checked+=1
            if a['target'] != b['source'] or a['degree']+b['degree']>4:
                assert not expected,(i,j,'zero product has nonzero table entry')
                continue
            key=(a['degree']+b['degree'],a['source'],b['target'])
            ids,piv=reducers[key]
            raw: dict[PathWord,Fraction]={}
            for p,c in polys[i].items():
                for q,d in polys[j].items():
                    raw[p+q]=raw.get(p+q,F(0))+c*d
            for k,c in expected.items():
                target=basis[k]
                assert (target['degree'],target['source'],target['target'])==key
                for p,d in polys[k].items():
                    raw[p]=raw.get(p,F(0))-c*d
            remainder=reduce_vector({ids[p]:c for p,c in raw.items() if c},piv)
            assert not remainder,(i,j,remainder)
            nonzero+=bool(expected)
    assert products_checked==6889 and nonzero==480
    return {'status':'PASS', 'arithmetic':'exact fractions.Fraction',
            'method':'direct homogeneous reduction; no potential or imported calculation code',
            'path_counts':path_counts,'relation_ranks':relation_ranks,
            'quotient_dimensions':quotient_dimensions,'relation_count':53,
            'independent_basis_vectors':83,'ordered_products_checked':products_checked,
            'nonzero_products':nonzero,
            'table_sha256':hashlib.sha256(data_path.read_bytes()).hexdigest(),
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'python':sys.version,
            'scope':'The fixed quadratic algebra and complete table; not the VOA comparison or an independent human review.'}


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=Path(__file__).with_name('finite_multiplication.json'))
    args=parser.parse_args()
    print(json.dumps(verify(args.data),indent=2))

if __name__=='__main__':
    main()
