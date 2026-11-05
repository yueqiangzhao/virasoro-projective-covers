#!/usr/bin/env python3
"""Exact finite checks for the conditional SL2 transport formula.

Program objects are coefficient dictionaries in the monomial basis
v[n,i] = x**(n-i) y**i of Sym^n(C^2), with rational coefficients.
No VOA or projective-module identification is asserted by this computation.
Only Python's standard library is used; one serial process, no random tests.
"""

from fractions import Fraction
from math import comb, factorial
from functools import lru_cache
import json
import platform
import sys


def fall(n, a):
    if n < 0 or a < 0 or a > n:
        return 0
    return factorial(n) // factorial(n - a)


def admissible(j, k, r):
    return min(j, k, r) >= 0 and abs(k-r) <= j <= k+r and (k+r-j) % 2 == 0


@lru_cache(None)
def embedding(j, k, r, m):
    """I^r_{j,k}(v[j,m]), returned as ((k-index,r-index), coefficient)."""
    if not admissible(j, k, r) or not 0 <= m <= j:
        return ()
    d = (k+r-j)//2
    out = {}
    prefactor = Fraction(factorial(j-m), factorial(j))
    for a in range(d+1):
        for b in range(m+1):
            coeff = prefactor * (-1)**a * comb(d,a) * comb(m,b)
            coeff *= fall(k-a,b) * fall(r-d+a,m-b)
            idx = (a+b, d-a+m-b)
            if coeff:
                out[idx] = out.get(idx, Fraction(0)) + coeff
    return tuple((i,c) for i,c in out.items() if c)


@lru_cache(None)
def transvectant(s, r, t, u, v):
    """Coefficient of v[t,u+v-h] in T^t_{s,r}(v[s,u] tensor v[r,v])."""
    if not admissible(t, s, r) or not (0 <= u <= s and 0 <= v <= r):
        return 0
    h = (s+r-t)//2
    if not 0 <= u+v-h <= t:
        return 0
    return sum((-1)**z * comb(h,z) * fall(s-u,h-z) * fall(u,z)
               * fall(r-v,z) * fall(v,h-z) for z in range(h+1))


def recoupling(j, k, ell, r, s, t):
    if not all((admissible(j,k,r), admissible(k,ell,s),
                admissible(t,s,r), admissible(j,ell,t))):
        return Fraction(0)
    d = (k+r-j)//2
    e = (ell+s-k)//2
    return sum((Fraction((-1)**a * comb(d,a) * fall(s-e,a), fall(k,a))
                * transvectant(s,r,t,e+a,d-a) for a in range(d+1)), Fraction(0))


def composition(j, k, ell, r, s, t, m):
    """(id tensor T_sr^t) (I_kl^s tensor id) I_jk^r on v[j,m]."""
    if not all((admissible(j,k,r), admissible(k,ell,s), admissible(t,s,r))):
        return {}
    h = (s+r-t)//2
    out = {}
    for (ki,ri), c in embedding(j,k,r,m):
        for (li,si), b in embedding(k,ell,s,ki):
            a = transvectant(s,r,t,si,ri)
            idx = (li,si+ri-h)
            if a:
                out[idx] = out.get(idx,Fraction(0)) + c*b*a
    return {i:c for i,c in out.items() if c}


def main():
    bound = 10
    channels = vectors = 0
    failures = []
    for j in range(bound+1):
        for k in range(bound+1):
            for r in range(3):
                if not admissible(j,k,r):
                    continue
                for ell in range(bound+1):
                    for s in range(3):
                        if not admissible(k,ell,s):
                            continue
                        for t in range(abs(s-r),s+r+1,2):
                            channels += 1
                            c = recoupling(j,k,ell,r,s,t)
                            for m in range(j+1):
                                vectors += 1
                                lhs = composition(j,k,ell,r,s,t,m)
                                rhs = {idx:c*a for idx,a in embedding(j,ell,t,m) if c*a}
                                if lhs != rhs:
                                    failures.append([j,k,ell,r,s,t,m,str(lhs),str(rhs)])
    return_cases = []
    for j in range(bound+1):
        plus0 = recoupling(j,j+1,j,1,1,0)
        plus2 = recoupling(j,j+1,j,1,1,2)
        assert plus0 == Fraction(j+2,j+1)
        assert plus2 == Fraction(j,j+1)
        minus0 = recoupling(j,j-1,j,1,1,0)
        minus2 = recoupling(j,j-1,j,1,1,2)
        assert minus0 == (-1 if j else 0)
        assert minus2 == (1 if j else 0)
        return_cases.append({'j':j, 'via_j_plus_1':[str(plus0),str(plus2)],
                             'via_j_minus_1':[str(minus0),str(minus2)]})
    report = {
        'claim':'Polynomial-coordinate recoupling equality on every monomial in each tested channel',
        'evidence_class':'Exact finite rational computation; not a uniform proof or VOA identification',
        'program_to_mathematics':'Dictionary keys are monomial indices in Sym^ell(C^2) tensor Sym^t(C^2); values are Fraction coefficients',
        'algorithm':'Expand both normalized embeddings, apply derivative transvectant, compare all coefficients with the closed finite sum',
        'coverage':{'j_k_ell':[0,bound], 'r_s':[0,2], 't':'all Clebsch-Gordan channels 0..r+s, including output-inadmissible channels'},
        'tested_channels':channels, 'tested_input_monomials':vectors,
        'failures':failures, 'return_channel_examples':return_cases,
        'runtime':{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
                   'arithmetic':'fractions.Fraction and arbitrary-precision integers','parallelism':1},
    }
    print(json.dumps(report,ensure_ascii=False,indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
