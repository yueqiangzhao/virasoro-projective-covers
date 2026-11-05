"""Exact formal checks for Nakano's four-lattice projection maps.

No VOA modes or screening integrals are approximated or simulated.  Each
operator-valued block is an independent formal symbol.  The mathematical
bridge to the source formula is proved in actual_maps_continuation.tex.
"""
from fractions import Fraction as Q
from itertools import product

V = tuple(product(range(2), repeat=2))
name = {(0,0):'A',(1,1):'B',(1,0):'C',(0,1):'D'}
x = lambda i: (1-i[0],i[1])
y = lambda i: (i[0],1-i[1])
bar = lambda i: (1-i[0],1-i[1])

def k(i,v):
    if v in (i,bar(i)): return 0
    return -1 if v==x(bar(i)) else 1

# A block is (operator label, exponent of d).  Absent blocks are zero.
def action(i):
    s=bar(i)
    M={(v,v):(f'Y_{name[v]}',0) for v in V}
    for m in (x(s),y(s)):
        M[m,s]=(f'D_{name[m]}_{name[s]}',1)
        M[i,m]=(f'D_{name[i]}_{name[m]}',0)
    M[i,s]=(f'H_{name[i]}',0)
    return M

def arrow(j,i):
    assert j in (x(i),y(i))
    return {bar(i):1,j:0}  # diagonal scalars d^exponent

checked=0
for i in V:
    for j in (x(i),y(i)):
        f=arrow(j,i)
        lhs={(v,u):(op,p+f[v]) for (v,u),(op,p) in action(i).items() if v in f}
        rhs={(v,u):(op,p+f[u]) for (v,u),(op,p) in action(j).items() if u in f}
        assert lhs==rhs, (name[i],name[j],lhs,rhs)
        weights={k(j,v)-k(i,v) for v in f}
        assert weights==({-1} if j==x(i) else {1})
        checked+=1
    for v in V:
        # minus lattices have odd Fock coset parity, plus lattices even.
        assert (sum(v)+k(i,v)-sum(i))%2==0
    for (v,u),(op,p) in action(i).items():
        if op.startswith('D_'):
            screening_shift = 1 if v==x(u) else -1
            assert k(i,v)-k(i,u)+screening_shift==0
        elif op.startswith('H_'):
            assert k(i,v)==k(i,u)

def compose_diag(g,f):
    return {v:g[v]+f[v] for v in g.keys() & f.keys()}

for i in V:
    via_x=compose_diag(arrow(bar(i),x(i)),arrow(x(i),i))
    via_y=compose_diag(arrow(bar(i),y(i)),arrow(y(i),i))
    assert via_x==via_y=={bar(i):1}
    for j in (x(i),y(i)):
        assert compose_diag(arrow(i,j),arrow(j,i))=={}

# Two path spaces, ordered C_{++},C_{+-},C_{-+},C_{--},D_{++},...
labels=list(product((1,-1),repeat=2))*2

def act(v,raising):
    out=[Q(0)]*8
    for block in range(2):
        for a,(s,t) in enumerate(product((1,-1),repeat=2)):
            q=v[4*block+a]
            if not q: continue
            for pos,z in enumerate((s,t)):
                if z==(-1 if raising else 1):
                    zz=[s,t]; zz[pos]=-z
                    b=list(product((1,-1),repeat=2)).index(tuple(zz))
                    out[4*block+b]+=q
    return out

def rank(rows):
    rows=[list(map(Q,r)) for r in rows]; r=0
    for c in range(len(rows[0])):
        pivot=next((j for j in range(r,len(rows)) if rows[j][c]),None)
        if pivot is None: continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        q=rows[r][c]; rows[r]=[v/q for v in rows[r]]
        for j in range(len(rows)):
            if j!=r:
                q=rows[j][c]; rows[j]=[a-q*b for a,b in zip(rows[j],rows[r])]
        r+=1
        if r==len(rows): break
    return r

r0=[Q(0)]*8; r0[1]=1; r0[6]=-1
orbit=[r0,act(r0,True),act(r0,False),act(act(r0,False),True)]
assert rank(orbit)==4
commutation=[]
for a,(s,t) in enumerate(product((1,-1),repeat=2)):
    b=list(product((1,-1),repeat=2)).index((t,s))
    row=[Q(0)]*8; row[a]=1; row[4+b]=-1
    commutation.append(row)
assert rank(orbit+commutation)==4
print('8 full formal block intertwining identities: PASS')
print('4 mixed diagonal compositions = d*pi_bar(i): PASS')
print('8 reverse same-axis compositions = 0: PASS')
print('Cartan semilinearity, central parity and arrow weights: PASS')
print('SL2 orbit of mixed seed has rank 4 and equals tensor-flip relations: PASS')
print('Scope: exact formal/rational identities; no analytic screening certification.')
