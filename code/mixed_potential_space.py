from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
V='TABCD'
a=[]
for x,y,m in [('T','A',1),('T','B',1),('A','C',2),('A','D',2),('B','C',2),('B','D',2)]:
 for s,t in [(x,y),(y,x)]:
  for c in ([0] if m==1 else [1,-1]): a.append((s,t,c))
idx={x:i for i,x in enumerate(a)}
out={v:[i for i,x in enumerate(a) if x[0]==v] for v in V}
def paths(s,n):
 if n==0:return [()]
 return [(i,)+p for i in out[s] for p in paths(a[i][1],n-1)]
def end(p):return a[p[-1]][1]
def canon(p):return min(p[k:]+p[:k] for k in range(len(p)))
cycles=[p for s in V for p in paths(s,4) if end(p)==s]
cs=sorted({canon(p) for p in cycles if sum(a[i][2] for i in p)==0}); ci={x:i for i,x in enumerate(cs)}
def val(p):return ci.get(canon(p))
rows=[]
def rel(ps):
 r={}
 for p,c in ps:
  j=val(p)
  if j is not None:r[j]=r.get(j,F(0))+c
 r={j:c for j,c in r.items() if c}
 if r:rows.append(r)
for p in cycles:
 if sum(a[i][2] for i in p)==-2:
  rel([(p[:k]+(idx[(a[i][0],a[i][1],1)],)+p[k+1:],F(1)) for k,i in enumerate(p) if a[i][2]==-1])
 # every return triplet zero (at any vertex, cyclicity takes care)
 if a[p[0]][0]==a[p[1]][1] and a[p[0]][2] and a[p[1]][2]:
  u,v=a[p[0]][2],a[p[1]][2]
  if u==v:rel([(p,F(1))])
  else:
   q=(idx[(a[p[0]][0],a[p[0]][1],v)],idx[(a[p[1]][0],a[p[1]][1],u)])+p[2:]
   rel([(p,F(1)),(q,F(1))])
piv={}
for r in rows:
 while r:
  j=min(r)
  if j in piv:
   c=r[j]
   for k,v in piv[j].items():
    r[k]=r.get(k,F(0))-c*v
    if not r[k]:del r[k]
  else:
   c=r[j];r={k:v/c for k,v in r.items()};piv[j]=r;break
free=[j for j in range(len(cs)) if j not in piv]
basis=[]
for t in free:
 v={t:F(1)}
 for j in sorted(piv,reverse=True):
  c=-sum(w*v.get(k,0) for k,w in piv[j].items() if k!=j)
  if c:v[j]=c
 basis.append(v)
print('arrows',len(a),'closed words',len(cycles),'weight0 cyclic orbits',len(cs),'constraints rank',len(piv),'potential dim',len(basis))
def coeff(p):
 j=val(p)
 return [b.get(j,F(0)) for b in basis] if j is not None else [F(0)]*len(basis)
def fmt(c):
 return ' '.join(f'{"+" if v >= 0 else ""}{v}z{i}' for i,v in enumerate(c) if v) or '0'
def name(p):return ''.join(a[i][0]+('+' if a[i][2]==1 else '-' if a[i][2]==-1 else '') for i in p)+end(p)
for i,b in enumerate(basis):
 print('\nz'+str(i),[(name(cs[j]),str(v)) for j,v in b.items()])
# Loop scalar path bases choose (+,-) for U1 wedge representatives. rank invariant under scaling.
for s in V:
 pp=[p for p in paths(s,2) if end(p)==s and (a[p[0]][2],a[p[1]][2]) in [(0,0),(1,-1)]]
 print('\nLOOP',s,[name(p) for p in pp])
 for p in pp:print([fmt(coeff(p+q)) for q in pp])
# cross blocks weight+- for U1 irreps, ++ U2, wedge U0
for s,t in [('A','B'),('C','D'),('T','C'),('T','D')]:
 pp=[p for p in paths(s,2) if end(p)==t];qq=[p for p in paths(t,2) if end(p)==s]
 print('\nCROSS',s,t)
 for p in pp:
  print(name(p),[fmt(coeff(p+q)) for q in qq])
 print('columns',[name(q) for q in qq])
Path(__file__).with_suffix('.json').write_text(json.dumps({'arrows':a,'cyclic_orbits':cs,'basis':[{str(j):str(v) for j,v in b.items()} for b in basis]},indent=2))
