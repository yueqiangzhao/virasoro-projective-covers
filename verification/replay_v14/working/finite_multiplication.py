"""An explicit equivariant 83-element basis and complete multiplication table.

The input is the fixed quadratic candidate, not an unknown End algebra.
All coordinates are exact Fractions, obtained by its degree-four pairing.
The table is then checked for associativity and sl2 equivariance.
"""
from contextlib import redirect_stdout
from io import StringIO
with redirect_stdout(StringIO()):
    import mixed_candidate as q
from fractions import Fraction as F
from math import factorial,comb
from pathlib import Path
import json,hashlib,sys

def add(*terms):
    result={}
    for coeff,poly in terms:
        for p,x in poly.items():
            result[p]=result.get(p,F(0))+coeff*x
            if not result[p]:del result[p]
    return result
def concat(P,Q):
    result={}
    for p,x in P.items():
        for t,y in Q.items():
            result[p+t]=result.get(p+t,F(0))+x*y
    return {p:x for p,x in result.items() if x}
def path(vs,ws):
    return {tuple(q.idx[(vs[i],vs[i+1],ws[i])] for i in range(len(ws))):F(1)}
def trace(P):return sum((x*q.phi(p) for p,x in P.items()),F(0))
def inverse(M):
    n=len(M); a=[list(row)+[F(int(i==j)) for j in range(n)] for i,row in enumerate(M)]
    for j in range(n):
        k=next(k for k in range(j,n) if a[k][j])
        a[j],a[k]=a[k],a[j];c=a[j][j];a[j]=[x/c for x in a[j]]
        for i in range(n):
            if i!=j:
                c=a[i][j]
                if c:a[i]=[x-c*y for x,y in zip(a[i],a[j])]
    return [row[n:] for row in a]

mults=[];basis=[]
def multiplet(name,degree,s,t,r,polys):
    assert len(polys)==r+1
    mi=len(mults); ids=[]
    for m,P in enumerate(polys):
        ids.append(len(basis))
        basis.append({'name':name+':'+str(m),'degree':degree,'source':s,'target':t,'r':r,'m':m,'multiplet':mi,'poly':P})
    mults.append({'name':name,'degree':degree,'source':s,'target':t,'r':r,'basis_ids':ids})

for s in q.V:multiplet('e_'+s,0,s,s,0,[{():F(1)}])
arrow_mults=[]
for s in q.V:
    for t in q.V:
        aa=[a for a in q.out[s] if q.a[a][1]==t]
        if not aa:continue
        r=0 if len(aa)==1 else 1
        arrow_mults.append((s,t,r,[{(a,):F(1)} for a in aa]))
        multiplet('a_'+s+t,1,s,t,r,arrow_mults[-1][3])

tloop=path('TAT',[0,0]);x={};y={}
for s in 'AB':x[s]=path(s+'C'+s,[1,-1]);y[s]=path(s+'D'+s,[1,-1])
for s in 'CD':x[s]=path(s+'A'+s,[1,-1]);y[s]=path(s+'B'+s,[1,-1])
multiplet('t',2,'T','T',0,[tloop])
for s in 'ABCD':
    multiplet('x_'+s,2,s,s,0,[x[s]]);multiplet('y_'+s,2,s,s,0,[y[s]])
for s,u,t in [('A','C','B'),('B','C','A'),('C','A','D'),('D','A','C')]:
    pp=path(s+u+t,[1,-1]);pm=path(s+u+t,[-1,1])
    multiplet('b0_'+s+t,2,s,t,0,[add((1,pp),(-1,pm))])
    multiplet('b2_'+s+t,2,s,t,2,[path(s+u+t,[1,1]),add((F(1,2),pp),(F(1,2),pm)),path(s+u+t,[-1,-1])])
for s,t in [('T','C'),('C','T'),('T','D'),('D','T')]:
    multiplet('b1_'+s+t,2,s,t,1,[path(s+'A'+t,[0,w] if s=='T' else [w,0]) for w in [1,-1]])

Z={'T':tloop,'A':add((1,x['A']),(-3,y['A'])),'B':add((1,x['B']),(1,y['B'])),
   'C':add((1,x['C']),(-3,y['C'])),'D':add((1,x['D']),(1,y['D']))}
for s,t,r,pp in arrow_mults:multiplet('za_'+s+t,3,s,t,r,[concat(Z[s],P) for P in pp])
socle_traces={}
for s in q.V:
    P=concat(Z[s],Z[s]);k=trace(P);assert k
    socle_traces[s]=k
    multiplet('s_'+s,4,s,s,0,[add((1/k,P))])
assert len(basis)==83 and len(mults)==55

blocks={}
for i,b in enumerate(basis):blocks.setdefault((b['degree'],b['source'],b['target']),[]).append(i)
pair_inverse={}
for key,ids in blocks.items():
    d,s,t=key;dual=blocks[(4-d,t,s)]
    M=[[trace(concat(basis[i]['poly'],basis[j]['poly'])) for i in ids] for j in dual]
    pair_inverse[key]=(ids,dual,inverse(M))

def coordinates(P,key):
    if key not in pair_inverse:return {}
    ids,dual,inv=pair_inverse[key]
    rhs=[trace(concat(P,basis[j]['poly'])) for j in dual]
    vals=[sum((a*b for a,b in zip(row,rhs)),F(0)) for row in inv]
    return {i:c for i,c in zip(ids,vals) if c}

table={}
for i,a in enumerate(basis):
    for j,b in enumerate(basis):
        if a['target']!=b['source'] or a['degree']+b['degree']>4:continue
        out=coordinates(concat(a['poly'],b['poly']),(a['degree']+b['degree'],a['source'],b['target']))
        if out:table[i,j]=out
def mul(X,Y):
    out={}
    for i,x in X.items():
        for j,y in Y.items():
            for k,c in table.get((i,j),{}).items():out[k]=out.get(k,F(0))+x*y*c
    return {k:c for k,c in out.items() if c}
def vsum(X,Y):
    out=dict(X)
    for k,c in Y.items():out[k]=out.get(k,F(0))+c
    return {k:c for k,c in out.items() if c}
def act(X,raising):
    out={}
    for i,x in X.items():
        b=basis[i];m=b['m'];r=b['r'];coef=m if raising else r-m
        if coef:
            j=mults[b['multiplet']]['basis_ids'][m-1 if raising else m+1]
            out[j]=out.get(j,F(0))+x*coef
    return {k:c for k,c in out.items() if c}

assoc=equiv=0
for i,a in enumerate(basis):
    for j,b in enumerate(basis):
        if a['target']!=b['source']:continue
        for raising in [True,False]:
            assert act(table.get((i,j),{}),raising)==vsum(mul(act({i:F(1)},raising),{j:F(1)}),mul({i:F(1)},act({j:F(1)},raising))),(i,j,raising)
            equiv+=1
        for k,c in enumerate(basis):
            if b['target']!=c['source'] or a['degree']+b['degree']+c['degree']>4:continue
            assert mul(table.get((i,j),{}),{k:F(1)})==mul({i:F(1)},table.get((j,k),{})),(i,j,k)
            assoc+=1

def fall(n,k):return factorial(n)//factorial(n-k) if 0<=k<=n else 0
def tau(s,r,t,u,v):
    if (s+r-t)%2 or not abs(s-r)<=t<=s+r:return F(0)
    h=(s+r-t)//2
    if not 0<=u+v-h<=t:return F(0)
    return F(sum((-1)**z*comb(h,z)*fall(s-u,h-z)*fall(u,z)*fall(r-v,z)*fall(v,h-z) for z in range(h+1)))
constants=[];channel_checks=0
for i,A in enumerate(mults):
    for j,B in enumerate(mults):
        if A['target']!=B['source'] or A['degree']+B['degree']>4:continue
        for k,C in enumerate(mults):
            if (C['degree'],C['source'],C['target'])!=(A['degree']+B['degree'],A['source'],B['target']):continue
            s,r,t=A['r'],B['r'],C['r']
            if (s+r-t)%2 or not abs(s-r)<=t<=s+r:continue
            h=(s+r-t)//2;scalar=None
            for u in range(s+1):
                for v in range(r+1):
                    coeff=tau(s,r,t,u,v);m=u+v-h
                    actual=table.get((A['basis_ids'][u],B['basis_ids'][v]),{}).get(C['basis_ids'][m],F(0)) if 0<=m<=t else F(0)
                    if coeff and scalar is None:scalar=actual/coeff
            assert scalar is not None
            for u in range(s+1):
                for v in range(r+1):
                    product=table.get((A['basis_ids'][u],B['basis_ids'][v]),{})
                    for m,bi in enumerate(C['basis_ids']):
                        expected=scalar*tau(s,r,t,u,v) if m==u+v-h else F(0)
                        assert product.get(bi,F(0))==expected,(i,j,k,u,v,m)
                        channel_checks+=1
            if scalar:constants.append({'first':i,'second':j,'output':k,'coefficient':str(scalar)})

report={'scope':'Fixed abstract candidate Q; actual VOA identification is a separate proof',
 'convention':'first then second is chronological concatenation; usual path-algebra product is second * first',
 'multiplets':mults,
 'basis':[{k:v for k,v in b.items() if k!='poly'}|{'representative':[{'arrows':list(p),'coefficient':str(c)} for p,c in b['poly'].items()]} for b in basis],
 'arrows':q.a,'central_Z_socle_traces':{s:str(v) for s,v in socle_traces.items()},
 'multiplication':[{'first':i,'second':j,'terms':{str(k):str(v) for k,v in out.items()}} for (i,j),out in table.items()],
 'equivariant_constants':constants,'transvectant':'T^t_{s,r}, first input U_s then second input U_r, unnormalized determinant differential',
 'checks':{'basis_dimension':83,'multiplets':55,'nonzero_products':len(table),'nonzero_channels':len(constants),'associativity_composable_triples_degree_le_4':assoc,'equivariance_equations':equiv,'channel_coefficient_equations':channel_checks},
 'runtime':sys.version,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report['checks'],indent=2))
print('All other products and channels are zero. All arithmetic exact.')
